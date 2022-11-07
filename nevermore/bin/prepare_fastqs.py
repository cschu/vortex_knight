#!/usr/bin/env python3

import argparse
import itertools
import os
import pathlib
import re
import shutil
import subprocess
import sys

from collections import Counter


def check_pairwise(r1, r2):
	""" Checks if two sets of read files contain the same prefixes.

	Input:
	 - r1/r2: lists of tuples (prefix, filename) of R1/R2 paired-end read files

	Raises error if one prefix is not found in either r1 or r2.

	"""
	r1_r2 = tuple(prefix[:-1] for prefix, _ in itertools.chain(r1, r2))
	for prefix in set(r1_r2):
		if r1_r2.count(prefix) != 2:
			raise ValueError(f"Missing mates for prefix {prefix}.")


def transfer_file(source, dest, remote_input=False):
	""" Transfers a file depending on its location and compressed state.

	Input:
	 - path to source file
	 - path to destination file
	 - whether source file is considered to be located on a remote file system

	"""
	resolved_src = pathlib.Path(source).resolve()
	if source.endswith(".gz"):
		if remote_input:
			# if file is on remote file system, copy it to destination
			shutil.copyfile(resolved_src, dest)
		else:
			# if file is gzipped and on local fs, just symlink it
			pathlib.Path(dest).symlink_to(resolved_src)
	elif source.endswith(".bz2"):
		bz2_pr = subprocess.Popen(("bzip2", "-dc", resolved_src), stdout=subprocess.PIPE)
		with open(dest, "wt") as _out:
			subprocess.run(("gzip", "-c", "-"), stdin=bz2_pr.stdout, stdout=_out)
	else:
		# if file is not compressed, gzip it to destination
		with open(dest, "wt") as _out:
			subprocess.run(("gzip", "-c", resolved_src), stdout=_out)


def transfer_multifiles(files, dest, remote_input=False, compression=None):
	""" Transfers a set of files depending on their location and compressed state.

	Input:
	 - list of source file paths
	 - path to destination file
	 - whether source files are considered to be located on a remote file system
	 - the compression type of the files (supported: None, gz, bz2)
	"""
	print(f"transfer_multifiles: remote={remote_input} compression={compression} dest={dest}", file=sys.stderr, flush=True)

	if len(files) > 1:
		src_files = tuple(os.path.abspath(f) for f in files)  # tuple(f.resolve() for f in files)
		cat_cmd = ("cat", ) + src_files

		if compression == ".gz":
			# multiple gzip-compressed files can just be concatenated
			with open(dest, "wt") as _out:
				subprocess.run(cat_cmd, stdout=_out)
		elif compression == ".bz2":
			cat_pr = subprocess.Popen(cat_cmd, stdout=subprocess.PIPE)
			bz2_pr = subprocess.Popen(("bzip2", "-dc", "-"), stdin=cat_pr.stdout, stdout=subprocess.PIPE)
			with open(dest, "wt") as _out:
				subprocess.run(("gzip", "-c", "-"), stdin=bz2_pr.stdout, stdout=_out)
		else:
			# multiple uncompressed files will be cat | gzipped
			cat_pr = subprocess.Popen(cat_cmd, stdout=subprocess.PIPE)
			with open(dest, "wt") as _out:
				subprocess.run(("gzip", "-c", "-"), stdin=cat_pr.stdout, stdout=_out)
			
	else:
		transfer_file(files[0], dest, remote_input=remote_input)


def process_sample(
	sample, fastqs, output_dir,
	fastq_suffix_pattern,
	remove_suffix=None, remote_input=False,
):
	""" Checks if a set of fastq files in a directory is a valid collection
	and transfers files to a destination dir upon success.

	Input:
	 - sample_id
	 - list of fastq files
	 - path to output directory
	 - suffix to strip off from filenames (e.g. _001)
	 - whether fastq files are located on remote file system
	"""

	if len(fastqs) == 1:
		# remove potential "single(s)" string from single fastq file name prefix
		sample_sub = re.sub(r"[._]singles?", "", sample)
		# 20221018: and attach it at the end of the sample name
		# - this might be a temporary fix, but @93a73d0
		# single-end samples without .singles-suffix cause problems 
		# with fastqc results in the collate step
		sample = sample_sub + ".singles"
		sample_dir = os.path.join(output_dir, sample)
		pathlib.Path(sample_dir).mkdir(parents=True, exist_ok=True)

		dest = os.path.join(sample_dir, f"{sample}_R1.fastq.gz")
		transfer_file(fastqs[0], dest, remote_input=remote_input)

	else:

		# check if all fastq files are compressed the same way
		suffixes = Counter(
			f[f.rfind("."):] in (".gz", ".bz2") for f in fastqs
		)

		if len(suffixes) > 1:
			raise ValueError(f"sample: {sample} has mixed compressed and uncompressed input files. Please check.")

		if suffixes.most_common()[0][0]:
			# all compressed
			suffixes = Counter(
				f[f.rfind("."):] for f in fastqs
			)
			if len(suffixes) > 1:
				raise ValueError(f"sample: {sample} has mixed gzip and bzip2 files. Please check.")
			compression = suffixes.most_common()[0][0]
		else:
			compression = None

		# extract the file name prefixes
		prefixes = [re.sub(fastq_suffix_pattern, "", os.path.basename(f)) for f in fastqs]
		if remove_suffix:
			# remove suffix pattern if requested
			prefixes = [re.sub(remove_suffix + r"$", "", p) for p in prefixes]

		print("PRE", prefixes, file=sys.stderr)

		# partition fastqs into R1, R2, and 'other' sets
		r1 = [(p, f) for p, f in zip(prefixes, fastqs) if re.search(r"[._R]1$", p)]
		r2 = [(p, f) for p, f in zip(prefixes, fastqs) if re.search(r"[._R]2$", p)]
		others = sorted(list(set(fastqs).difference({f for _, f in r1}).difference({f for _, f in r2})))

		# check if R1/R2 sets have equal sizes or are empty
		# R1 empty: potential scRNAseq (or any protocol with barcode reads in R1)
		# R2 empty: typical single end reads with (R?)1 suffix
		assert len(r2) == 0 or len(r1) == 0 or (r1 and len(r1) == len(r2)), "R1/R2 sets are not of the same length"

		# if R1 and R2 are of equal size, check if the prefixes match
		if len(r1) == len(r2) and r1:
			check_pairwise(r1, r2)

		# sort R1/R2 for concatenation, get rid off prefixes
		r1 = sorted(f for _, f in r1)
		r2 = sorted(f for _, f in r2)

		print("R1", r1, file=sys.stderr)
		print("R2", r2, file=sys.stderr)
		print("others", others, file=sys.stderr, flush=True)

		sample_dir = os.path.join(output_dir, sample)

		if r1 or r2:

			pathlib.Path(sample_dir).mkdir(parents=True, exist_ok=True)

			if r1:
				# if R1 is not empty, transfer R1-files
				dest = os.path.join(sample_dir, f"{sample}_R1.fastq.gz")
				transfer_multifiles(r1, dest, remote_input=remote_input, compression=compression)
			if r2:
				# if R2 is not empty, transfer R2-files,
				# if R1 is empty, rename R2 to R1 so that files can be processed as normal single-end
				target_r = "R2" if r1 else "R1"
				dest = os.path.join(sample_dir, f"{sample}_{target_r}.fastq.gz")
				transfer_multifiles(r2, dest, remote_input=remote_input, compression=compression)

		if others:
			# if single-end reads exist,
			# transfer them to <sample>.singles
			# these will be processed independently and merged with the paired-end reads
			# at a later stage
			sample_dir = sample_dir + ".singles"
			pathlib.Path(sample_dir).mkdir(parents=True, exist_ok=True)
			dest = os.path.join(sample_dir, f"{sample}.singles_R1.fastq.gz")
			transfer_multifiles(others, dest, remote_input=remote_input, compression=compression)
		

def is_fastq(f, valid_fastq_suffixes, valid_compression_suffixes):
	""" Checks if a file is a fastq file (compressed or uncompressed.)

	Input:
	 - filename

	Output:
	 - true if file is fastq else false

	"""
	# valid_fastq_suffixes = (".fastq", ".fq", ".txt")
	# valid_compression_suffixes = (".gz", ".bz2")

	prefix, suffix = os.path.splitext(f)
	if suffix in valid_fastq_suffixes:
		return True
	if suffix in valid_compression_suffixes:
		_, suffix = os.path.splitext(prefix)
		return suffix in valid_fastq_suffixes
	return False


def main():
	ap = argparse.ArgumentParser()
	ap.add_argument("-i", "--input_dir", type=str, default=".")
	ap.add_argument("-o", "--output_dir", type=str, default="prepared_samples")
	ap.add_argument("-p", "--prefix", type=str, required=True)
	ap.add_argument("--remote-input", action="store_true")
	ap.add_argument("--remove-suffix", type=str, default=None)
	ap.add_argument("--valid-fastq-suffixes", type=str, default="fastq,fq")
	ap.add_argument("--valid-compression-suffixes", type=str, default="gz,bz2")

	args = ap.parse_args()

	valid_fastq_suffixes = tuple(f".{suffix}" for suffix in args.valid_fastq_suffixes.split(","))
	print(valid_fastq_suffixes)
	valid_compression_suffixes = tuple(f".{suffix}" for suffix in args.valid_compression_suffixes.split(","))
	print(valid_compression_suffixes)

	fastq_file_suffix_pattern = r"[._](" + \
		args.valid_fastq_suffixes.replace(",", "|") + \
		")([._](" + \
		args.valid_compression_suffixes.replace(",", "|") + \
		"))?$"

	fastq_mate_pattern = r"([._R][12])$"

	def collect_fastq_files(input_dir, valid_fastq_suffixes, valid_compression_suffixes):
		return sorted(
				os.path.join(input_dir, f)
				for f in os.listdir(input_dir)
				if is_fastq(f, valid_fastq_suffixes, valid_compression_suffixes)
			)

	try:
		pwd, dirs, _ = next(os.walk(args.input_dir))
	except StopIteration:
		raise ValueError(f"Could not find input directory {args.input_dir} ({os.path.abspath(args.input_dir)})")

	samples = {}
	pathlib.Path(args.output_dir).mkdir(parents=True, exist_ok=True)

	for sample_dir in dirs:
		sample, sample_dir = sample_dir, os.path.join(pwd, sample_dir)

		samples.setdefault(sample, []).extend(
			collect_fastq_files(sample_dir, valid_fastq_suffixes, valid_compression_suffixes)			
		)

	root_fastqs = collect_fastq_files(args.input_dir, valid_fastq_suffixes, valid_compression_suffixes)

	if samples and root_fastqs:
		raise ValueError("Found {len(root_fastqs)} fastq files in input directory together with {len(samples)} sample directories. Please check input data.")
	elif root_fastqs:
		for f in root_fastqs:
			# sample = re.sub(r"[._](fastq|fq|txt)([._](gz|bz2))?$", "", os.path.basename(f))
			sample = re.sub(fastq_file_suffix_pattern, "", os.path.basename(f))
			sample = re.sub(fastq_mate_pattern, "", sample)
			samples.setdefault(sample, []).append(f)

	# # collect all fastq files from input directory
	# # assumption: fastq files are sym-linked into input_dir from a sample/files directory tree
	# # i.e., from a nextflow.Channel()
	# fastqs = sorted(
	# 	f 
	# 	for f in os.listdir(args.input_dir)
	# 	if is_fastq(f)
	# )
	# assert fastqs, f"Could not find any fastq files in '{args.input_dir}'."

	# samples = {}

	# # resolve the symlinks and group by sample ids
	# for f in fastqs:
	# 	full_f = pathlib.Path(os.path.join(args.input_dir, f))
	# 	if full_f.is_symlink():
	# 		link_target = full_f.resolve()
	# 		sample, *_, _ = str(link_target).replace(args.prefix, "").lstrip("/").split("/")
	# 		samples.setdefault(sample, []).append(full_f)

	# # if there is only one sample, this points at a flat sample hierarchy (input_dir/{file1,file2,file3,...})
	# # in this case, we derive the sample names from the files: same prefix -> same sample
	# # in a flat sample hierarchy, we cannot have mixed file prefixes, i.e. no lanes, no preprocessed singles, etc.
	# if len(samples) == 1:
	# 	tmp_samples = {}
	# 	for f in samples[list(samples.keys())[0]]:
	# 		sample = re.sub(r"\.(fastq|fq|txt)(.(gz|bz2))?$", "", os.path.basename(f.name))
	# 		sample = re.sub(r"([._R][12])$?", "", sample)
	# 		tmp_samples.setdefault(sample, []).append(f)
	# 	samples.clear()
	# 	samples = tmp_samples		

	# check and transfer the files
	for sample, fastqs in samples.items():
		try:
			process_sample(
				sample, fastqs, args.output_dir,
				fastq_file_suffix_pattern,
				remove_suffix=args.remove_suffix, remote_input=args.remote_input
			)
		except Exception as e:
			raise ValueError(f"Encountered problems processing sample '{sample}': {e}.\nPlease check your file names.")


if __name__ == "__main__":
	main()