#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { fastqc } from "../modules/qc/fastqc"
include { multiqc } from "../modules/qc/multiqc"
include { bwa_mem_align } from "../modules/align/bwa"
include { merge_and_sort } from "../modules/align/helpers"
include { merge_single_fastqs } from "../modules/converters/merge_fastqs"

def config_dir = (projectDir.endsWith("nevermore")) ? "${projectDir}/config" : "${projectDir}/nevermore/config"


workflow nevermore_align {

	take:

		fastq_ch

	main:

		/*	route all single-read files into a common channel */

		single_ch = fastq_ch
			.filter { it[0].is_paired == false }
			.map { sample, fastq ->
				def meta = [:]
				meta.id = fastq.name.replaceAll(/_R1.fastq.gz$/, "")
				meta.is_paired = false
				meta.merged = false
				return tuple(meta, fastq)
			}

		/*	route all paired-end read files into a common channel */

		paired_ch = fastq_ch
			.filter { it[0].is_paired == true }
			.map { sample, fastq ->
				def meta = [:]
				meta.id = sample.id
				meta.is_paired = true
				meta.merged = true
				return tuple(meta, fastq)
			}

		/*	group all single-read files by sample and route into merge-channel */

		merged_single_ch = single_ch
			.map { sample, fastq ->
				return tuple(
					sample.id.replaceAll(/.(orphans|singles|chimeras)$/, ".singles"),
					fastq
				)
			}
			.groupTuple(sort: true)
			.map { sample_id, files ->
				def meta = [:]
				meta.id = sample_id
				meta.is_paired = false
				meta.merged = true
				return tuple(meta, files)
			}

		/*	then merge single-read file groups into single files */

		merge_single_fastqs(merged_single_ch)

		/* 	take all single-read files except for the qc-survivors,
			concat with merged single-read files (takes care of single-end qc-survivors),
			concat with paired-end files,
			and route them into a channel for post-qc fastqc analysis
		*/

		fastqc_in_ch = single_ch
			.filter { ! it[0].id.endsWith(".singles") }
			.map { sample, fastq ->
				def meta = [:]
				meta.id = fastq.name.replaceAll(/_R1.fastq.gz$/, "")
				meta.is_paired = false
				meta.merged = false
				return tuple(meta, fastq)
			}
			.concat(merge_single_fastqs.out.fastq)
			.concat(paired_ch)

		/*	perform post-qc fastqc analysis and generate multiqc report on merged single-read and paired-end sets */

		fastqc(fastqc_in_ch, "qc")

		multiqc(
			fastqc.out.stats
				.filter { it[0].merged == true || it[0].is_paired == true }
				.map { sample, report -> return report }
				.collect(),
			"${config_dir}/multiqc.config",
			"qc"
		)

		/*	align merged single-read and paired-end sets against reference */

		bwa_mem_align(
			paired_ch.concat(merge_single_fastqs.out.fastq),
			params.reference,
			true
		)

		/*	merge paired-end and single-read alignments into single per-sample bamfiles */

		aligned_ch = bwa_mem_align.out.bam
			.map { sample, bam ->
				sample_id = sample.id.replaceAll(/.(orphans|singles|chimeras)$/, "")
				return tuple(sample_id, bam)
			}
			.groupTuple(sort: true)

		merge_and_sort(aligned_ch, true)


	emit:

		alignments = merge_and_sort.out.bam
		read_counts = fastqc.out.counts
		aln_counts = merge_and_sort.out.flagstats

}

