#!/usr/bin/env python3

import argparse
import itertools
import os
import pathlib
import re
import sys


def check_pairwise(r1, r2):
	pairs = {}
	for p1, p2 in itertools.product(
		tuple(p[:-1] for p, _ in r1),
		tuple(p[:-1] for p, _ in r2)
	):
		pairs.setdefault(p1, set()).add(1)
		pairs.setdefault(p2, set()).add(2)

	for p, counts in pairs.items():
		if len(counts) < 2:
			raise ValueError(f"Missing mates for prefix {p}, files={str(counts)}")
		elif len(counts) > 2:
			raise ValueError(f"Too many files for prefix {p}, files={str(counts)}")
		else:
			...



def process_sample(input_dir, sample_id, fastqs, output_dir, remove_suffix=None):

	print("#!/usr/bin/bash")

	print("RESUF", remove_suffix, file=sys.stderr)

	if len(fastqs) == 1:
		sample_id = re.sub(r"[._]singles?", "", sample_id) + ".singles"
		pathlib.Path(os.path.join(output_dir, sample_id)).mkdir(parents=True, exist_ok=True)
		print(f"cp {os.path.join(input_dir, fastqs[0])} {os.path.join(output_dir, sample_id, sample_id)}_R1.fastq.gz")




	else:
		prefixes = [re.sub(r"\.(fastq|fq).gz$", "", f) for f in fastqs]
		if remove_suffix:
			prefixes = [re.sub(remove_suffix + r"$", "", p) for p in prefixes]

		print("PRE", prefixes, file=sys.stderr)

		r1 = [(p, f) for p, f in zip(prefixes, fastqs) if p.endswith("1")]
		r2 = [(p, f) for p, f in zip(prefixes, fastqs) if p.endswith("2")]
		others = set(fastqs).difference({f for _, f in r1}).difference({f for _, f in r2})

		assert len(r2) == 0 or len(r1) == len(r2), "R1/R2 sets are not of the same length"
		check_pairwise(r1, r2)

		r1 = [os.path.join(input_dir, f) for f in sorted(f for _, f in r1)]
		r2 = [os.path.join(input_dir, f) for f in sorted(f for _, f in r2)]


		print("R1", r1, file=sys.stderr)
		print("R2", r2, file=sys.stderr)

		r1_cmd = " ".join(["cat"] + r1) + f" > {os.path.join(output_dir, sample_id, sample_id)}_R1.fastq.gz"
		print(r1_cmd)

		if r2:
			r2_cmd = " ".join(["cat"] + r2) + f" > {os.path.join(output_dir, sample_id, sample_id)}_R2.fastq.gz"
			print(r2_cmd)

		if others:
			sample_id = sample_id + ".singles"
			pathlib.Path(os.path.join(output_dir, sample_id)).mkdir(parents=True, exist_ok=True)
			other_cmd = " ".join(["cat"] + list(others)) + f" > {os.path.join(output_dir, sample_id, sample_id)}_R1.fastq.gz"
			print(other_cmd)


		# ['HD_Path_4_S81_L001_R1_001', 'HD_Path_4_S81_L001_R2_001', 'HD_Path_4_S81_L002_R1_001', 'HD_Path_4_S81_L002_R2_001']
		# ['HD_Path_4_S81_L001_R1', 'HD_Path_4_S81_L001_R2', 'HD_Path_4_S81_L002_R1', 'HD_Path_4_S81_L002_R2']







def main():

	ap = argparse.ArgumentParser()
	ap.add_argument("-i", "--input_dir", type=str, default=".")
	ap.add_argument("-o", "--output_dir", type=str, default="prepared_samples")
	ap.add_argument("-s", "--sample_id", type=str, required=True)
	ap.add_argument("--remove-suffix", type=str, default=None)

	args = ap.parse_args()

	pathlib.Path(os.path.join(args.output_dir, args.sample_id)).mkdir(parents=True, exist_ok=True)


	walk = os.walk(args.input_dir)
	try:
		cwd, dirs, files = next(walk)
	except StopIteration:
		raise StopIteration(f"Invalid directory: {args.input_dir}")

	fastqs = sorted(f for f in files if f.endswith(".fq.gz") or f.endswith(".fastq.gz"))
	assert fastqs, "Could not find any fastq files."

	print(cwd, fastqs, file=sys.stderr)
	try:
		process_sample(cwd, args.sample_id, fastqs, args.output_dir, remove_suffix=args.remove_suffix)
	except Exception as e:
		raise ValueError(f"Encountered problems processing sample '{args.sample_id}': {e}.\nPlease check your file names.")





if __name__ == "__main__":
	main()