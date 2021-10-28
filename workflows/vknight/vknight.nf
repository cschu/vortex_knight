#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { kraken2 } from  "../../modules/vknight/profilers/kraken2"
include { mtags_extract; mtags_annotate; mtags_merge } from "../../modules/vknight/profilers/mtags"
include { motus2 } from  "../../modules/vknight/profilers/motus2"
include { mapseq; collate_mapseq_tables } from "../../modules/vknight/profilers/mapseq"
include { pathseq } from "../../modules/vknight/profilers/pathseq"


if (!params.publish_mode) {
	params.publish_mode = "symlink"
}

if (!params.output_dir) {
	params.output_dir = "vknight_out"
}

output_dir = "${params.output_dir}"

/*if (params.motus_database) {
	motus_database = "-db ${params.motus_database}"
} else {
	motus_database = ""
}*/

if (!params.motus2_min_length) {
	params.motus2_min_length = 30
}

if (!params.motus2_n_marker_genes) {
	params.motus2_n_marker_genes = 1
}

if (!params.pathseq_min_clipped_read_length) {
	params.pathseq_min_clipped_read_length = 31
}

def run_kraken2 = (!params.skip_kraken2 || params.run_kraken2);
def run_mtags = (!params.skip_mtags || params.run_mtags);
def run_mapseq = (run_mtags && (!params.skip_mapseq || params.run_mapseq) && params.mapseq_bin)
def run_motus2 = (!params.skip_motus2 || params.run_motus2);
def run_pathseq = (!params.skip_pathseq || params.run_pathseq);
def run_count_reads = (!params.skip_counts || params.run_counts);
def convert_fastq2bam = (run_pathseq || run_count_reads);


process collate_data {
	publishDir "$output_dir", mode: params.publish_mode

	input:
	path(outputs)

	output:
	path("combined_profiles_rds/*.rds"), emit: collated_data

	script:
	def kraken2_output = run_kraken2 ? "--kraken2_res_path kraken2/" : ""
	def motus_output = run_motus2 ? "--mOTUs_res_path motus/" : ""
	def mtags_output = run_mtags ? "--mTAGs_res_path mtags_tables/" : ""
	def counts_output = run_count_reads ? "--libsize_res_path libsize/ --lib_layout_res_path lib_layout/" : ""
	def pathseq_output = run_pathseq ? "--PathSeq_res_path pathseq/" : ""
	"""
	mkdir -p combined_profiles_rds/
	mkdir -p {kraken2,pathseq,motus,libsize,lib_layout,otu_tables,mtags_tables}

	find \$(pwd) -maxdepth 1 -name '*kraken2_report.txt' -exec ln -sf {} kraken2/ \\;
	find \$(pwd) -maxdepth 1 -name '*pathseq.scores' -exec ln -sf {} pathseq/ \\;
	find \$(pwd) -maxdepth 1 -name '*motus.txt' -exec ln -sf {} motus/ \\;
	find \$(pwd) -maxdepth 1 -name '*libsize.txt' -exec ln -sf {} libsize/ \\;
	find \$(pwd) -maxdepth 1 -name '*is_paired.txt' -exec ln -sf {} lib_layout/ \\;
	find \$(pwd) -maxdepth 1 -name '*bac_ssu.tsv' -exec ln -sf {} otu_tables/ \\;
	find \$(pwd) -maxdepth 1 -name 'merged_profile.genus.tsv' -exec ln -sf {} mtags_tables/ \\;

	Rscript ${params.collate_script} \\
	${kraken2_output} \\
	${motus_output} \\
	${mtags_output} \\
	${pathseq_output} \\
	${counts_output} \\
	--out_folder combined_profiles_rds/
	"""
}


workflow bam_analysis {
	take:
		bam_ch

	main:
		out_ch = Channel.empty()
    	if (run_pathseq) {
        	pathseq(bam_ch)
			out_ch = out_ch.concat(pathseq.out.scores)
    	}

		out_ch = out_ch
			.map { sample, files -> return files }

	emit:
		results = out_ch
}


workflow fastq_analysis {
	take:
		fastq_ch

	main:
		out_ch = Channel.empty()

		if (run_kraken2) {
			kraken2(fastq_ch)
			out_ch = out_ch.concat(kraken2.out.kraken2_out)
		}
	
		if (run_motus2) {
			motus2(fastq_ch, params.motus_database)
			out_ch = out_ch.concat(motus2.out.motus_out)
		}

		out_ch = out_ch
			.map { sample, files -> return files }

		if (run_mtags) {
			mtags_extract(fastq_ch)
	
			mtags_annotate(mtags_extract.out.mtags_out)
	
			mtags_merge(mtags_annotate.out.mtags_bins.collect())

			out_ch = out_ch.concat(mtags_merge.out.mtags_tables)

			if (run_mapseq) {
				mapseq(mtags_extract.out.mtags_out)
	
				collate_mapseq_tables(mapseq.out.bac_ssu.collect())
			}
		}

	emit:
		results = out_ch
}


workflow {

	/* perform fastq-based analyses */

	if (run_kraken2) {
		kraken2(combined_fastq_ch)
	}

	if (run_motus2) {
		motus2(combined_fastq_ch)
	}

	if (run_mtags) {
		mtags_extract(combined_fastq_ch)

		mtags_annotate(mtags_extract.out.mtags_out)

		mtags_merge(mtags_annotate.out.mtags_bins.collect())

		if (run_mapseq) {
			mapseq(mtags_extract.out.mtags_out)

			collate_mapseq_tables(mapseq.out.bac_ssu.collect())
		}
	}

	/* collate data */

	if (params.collate_script != null && params.collate_script != "") {
		data_to_collate_ch = Channel.empty()

		if (run_kraken2) {
			data_to_collate_ch = data_to_collate_ch.concat(kraken2.out.kraken2_out)
		}

		if (run_count_reads) {
			data_to_collate_ch = data_to_collate_ch.concat(count_reads.out.counts)
				.concat(count_reads.out.is_paired)
		}

		if (run_motus2) {
			data_to_collate_ch = data_to_collate_ch.concat(motus2.out.motus_out)
		}

		if (run_pathseq) {
			data_to_collate_ch = data_to_collate_ch.concat(pathseq.out.scores)
		}

		data_to_collate_ch = data_to_collate_ch
			.map { sample, files -> return files }

		if (run_mtags) {
			data_to_collate_ch = data_to_collate_ch.concat(mtags_merge.out.mtags_tables)
		}

		collate_data(data_to_collate_ch.collect())
	}
}
