#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { fastq_input; bam_input } from "./nevermore/workflows/input"
include { vknight_main } from "./vknight/workflows/vknight"


if (params.input_dir && params.remote_input_dir) {
	log.info """
		Cannot process both --input_dir and --remote_input_dir. Please check input parameters.
	""".stripIndent()
	exit 1
} else if (!params.input_dir && !params.remote_input_dir) {
	log.info """
		Neither --input_dir nor --remote_input_dir set.
	""".stripIndent()
	exit 1
}

def input_dir = (params.input_dir) ? params.input_dir : params.remote_input_dir
def fastq_input_pattern = input_dir + "/" + "**[._]{fastq.gz,fq.gz}"
def bam_input_pattern = input_dir + "/" + "**.bam"

def run_kraken2 = (!params.skip_kraken2 || params.run_kraken2) && !params.amplicon_seq;
def run_mtags = (!params.skip_mtags || params.run_mtags);
def run_mapseq = (run_mtags && (!params.skip_mapseq || params.run_mapseq) && params.mapseq_bin)
def run_motus = (!params.skip_motus || params.run_motus) && !params.amplicon_seq;
def run_pathseq = (!params.skip_pathseq || params.run_pathseq) && !params.amplicon_seq;
def run_read_counter = (!params.skip_read_counter || params.run_read_counter)

def run_bam_analysis = run_pathseq && !params.amplicon_seq
def run_fastq_analysis = (run_kraken2 || run_mtags || run_mapseq || run_motus || run_read_counter) && !params.amplicon_seq
def run_amplicon_analysis = params.amplicon_seq


workflow {

	fastq_input(
		Channel.fromPath(fastq_input_pattern)
	)

	bam_input(
		Channel.fromPath(bam_input_pattern)
	)

	fastq_ch = fastq_input.out.fastqs.concat(bam_input.out.bamfiles)
	fastq_ch.view()

	vknight_main(fastq_ch)

}
