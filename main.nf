#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { fastq_input; bam_input } from "./nevermore/workflows/input"
include { vknight_main } from "./vknight/workflows/vknight"
include { nevermore_main } from "./nevermore/workflows/nevermore"


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

params.fastq_input_pattern = "**[._]{fastq.gz,fq.gz,fastq.bz2,fq.bz2}"
def fastq_input_pattern = input_dir + "/" + params.fastq_input_pattern
def bam_input_pattern = input_dir + "/" + "**.bam"

print "PARAMS: ${params}"

workflow {

	if (params.bam_input_pattern) {
		bam_input(
			Channel.fromPath(bam_input_pattern)
		)
		fastq_ch = bam_input.out.fastqs
	} else {
		fastq_input(
			// Channel.fromPath(input_dir + "/*", type: "dir")
			Channel.fromPath(fastq_input_pattern),
			Channel.of(null)
		)
		fastq_ch = fastq_input.out.fastqs
	}

	fastq_ch.dump(pretty: true, tag: "fastq_ch_dump")

	nevermore_main(fastq_ch)
	
	vknight_main(nevermore_main.out.fastqs)

}
