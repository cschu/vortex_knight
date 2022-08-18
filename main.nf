#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { fastq_input } from "./nevermore/workflows/input"
include { bam2fq } from "./nevermore/modules/converters/bam2fq"
include { fq2bam } from "./nevermore/modules/converters/fq2bam"
include { prepare_fastqs } from "./nevermore/modules/converters/prepare_fastqs"
include { nevermore_simple_preprocessing } from "./nevermore/workflows/nevermore"
include { amplicon_analysis; bam_analysis; fastq_analysis } from "./vknight/workflows/vknight"
include { classify_sample } from "./nevermore/modules/functions"
include { remove_host_kraken2; remove_host_kraken2_individual } from "./nevermore/modules/decon/kraken2"
include { flagstats } from "./nevermore/modules/stats"
include { collate_results } from "./vknight/modules/collate"

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

def get_basecounts = (!params.skip_basecounts || params.run_basecounts);
def convert_fastq2bam = (run_pathseq || get_basecounts);

def do_preprocessing = (!params.skip_preprocessing || params.run_preprocessing)

def run_bam_analysis = run_pathseq && !params.amplicon_seq
def run_fastq_analysis = (run_kraken2 || run_mtags || run_mapseq || run_motus || run_read_counter) && !params.amplicon_seq
def run_amplicon_analysis = params.amplicon_seq


process transfer_bams {
	input:
		path(bamfiles)
	output:
		path("bam/*.bam"), emit: bamfiles
	script:
		"""
		mkdir -p bam/
		find . -maxdepth 1 -type l -name '*.bam' | xargs -I {} readlink {} | xargs -I {} rsync -avP {} bam/
		"""
}

workflow bam_input {
	take:
		bam_ch
	main:
		transfer_bams(bam_ch.collect())
		bam_ch = transfer_bams.out.bamfiles
			.map { file ->
				def sample = file.name.replaceAll(/.bam$/, "")
				return tuple(sample, file)
			}
			.groupTuple(sort: true)
			.map { classify_sample(it[0], it[1]) }

		bam2fq(bam_ch)
		bam_ch = bam2fq.out.reads
			.map { classify_sample(it[0].id, it[1]) }
	emit:
		bamfiles = bam_ch
}

workflow vknight_main {
	take:
		fastq_ch
	main:
		results_ch = Channel.empty()

		if (do_preprocessing) {

			raw_fastq_ch = fastq_ch.out.fastqs.concat(bfastq_ch)

			nevermore_simple_preprocessing(raw_fastq_ch)

			preprocessed_ch = nevermore_simple_preprocessing.out.main_reads_out
			results_ch = results_ch
				.concat(nevermore_simple_preprocessing.out.raw_counts)
				.map { sample, files -> files }

			if (params.remove_host) {

				if (params.remove_host == "individual") {

					remove_host_kraken2_individual(nevermore_simple_preprocessing.out.main_reads_out, params.remove_host_kraken2_db)

					preprocessed_ch = remove_host_kraken2_individual.out.reads


				} else if (params.remove_host == "pair") {

					remove_host_kraken2(nevermore_simple_preprocessing.out.main_reads_out, params.remove_host_kraken2_db)

					preprocessed_ch = remove_host_kraken2.out.reads

				}

			}

		} else {

			preprocessed_ch = fastq_ch.out.fastqs
				.concat(bfastq_ch)

		}


		if (get_basecounts || run_bam_analysis) {

			fq2bam(preprocessed_ch)

			if (get_basecounts) {

				flagstats(fq2bam.out.reads)

				flagstat_results_ch = flagstats.out.flagstats
					.concat(flagstats.out.counts)
					.concat(flagstats.out.is_paired)
					.map { sample, files -> files }
				results_ch = results_ch.concat(flagstat_results_ch)

			}

			if (run_bam_analysis) {

				bam_analysis(fq2bam.out.reads)
				results_ch = results_ch.concat(bam_analysis.out.results)

			}

		}

		if (run_fastq_analysis) {

			fastq_analysis(preprocessed_ch)
			results_ch = results_ch.concat(fastq_analysis.out.results)

		}

		if (run_amplicon_analysis) {

			amplicon_analysis(preprocessed_ch)
			results_ch = results_ch.concat(amplicon_analysis.out.results)

		}

		if (!params.skip_collate) {
			collate_results(
				results_ch.collect(),
				"${projectDir}/scripts/ExtractProfiledCounts_210823.R",
				params.GTDB_markers
			)
		}
	
	emit:
		results = results_ch


}


workflow {

	fastq_input(
		Channel.fromPath(fastq_input_pattern)
	)

	bam_input(
		Channel.fromPath(bam_input_pattern)
	)

	fastq_ch = fastq_input.out.fastqs.concat(bam_input.out.bamfiles)
	fastq_ch.view()

	if (params.run_vknight) {
		vknight_main(fastq_ch)
	}


}

// 		/*bam2fq(bam_ch.out.bamfiles)

// 	bfastq_ch = bam2fq.out.reads
// 		.map { classify_sample(it[0].id, it[1]) } */

// 	/*fastq_ch = Channel
// 		//.fromPath(params.input_dir + "/" + "**.{fastq,fq,fastq.gz,fq.gz}")
// 		.fromPath(params.input_dir + "/" + "**[._]{fastq.gz,fq.gz}")
// 		.map { file ->
// 				def sample = file.name.replaceAll(/.?(fastq|fq)(.gz)?$/, "")
// 				sample = sample.replaceAll(/_R?[12]$/, "")
// 				return tuple(sample, file)
// 		}
// 		.groupTuple(sort: true)
//         .map { classify_sample(it[0], it[1]) }

// 	bam_ch = Channel
// 		.fromPath(params.input_dir + "/" + "**.bam")
// 		.map { file ->
// 			def sample = file.name.replaceAll(/.bam$/, "")
// 			return tuple(sample, file)
// 		}
// 		.groupTuple(sort: true)
// 		.map { classify_sample(it[0], it[1]) }

// 	bam2fq(bam_ch)

// 	bfastq_ch = bam2fq.out.reads
// 		.map { classify_sample(it[0].id, it[1]) }

// 	prepare_fastqs(fastq_ch)  */

// 	results_ch = Channel.empty()

// 	if (do_preprocessing) {

// 		// raw_fastq_ch = prepare_fastqs.out.reads.concat(bfastq_ch)
// 		raw_fastq_ch = fastq_ch.out.fastqs.concat(bfastq_ch)

// 		nevermore_simple_preprocessing(raw_fastq_ch)

// 		preprocessed_ch = nevermore_simple_preprocessing.out.main_reads_out
// 		results_ch = results_ch
// 			.concat(nevermore_simple_preprocessing.out.raw_counts)
// 			.map { sample, files -> files }

// 		if (params.remove_host) {

// 			if (params.remove_host == "individual") {

// 				remove_host_kraken2_individual(nevermore_simple_preprocessing.out.main_reads_out, params.remove_host_kraken2_db)

// 				preprocessed_ch = remove_host_kraken2_individual.out.reads


// 			} else if (params.remove_host == "pair") {

// 				remove_host_kraken2(nevermore_simple_preprocessing.out.main_reads_out, params.remove_host_kraken2_db)

// 				preprocessed_ch = remove_host_kraken2.out.reads

// 			}

// 		}

// 	} else {

// 		//preprocessed_ch = prepare_fastqs.out.reads
// 		preprocessed_ch = fastq_ch.out.fastqs
// 			.concat(bfastq_ch)

// 	}


// 	if (get_basecounts || run_bam_analysis) {

// 		fq2bam(preprocessed_ch)

// 		if (get_basecounts) {

// 	        flagstats(fq2bam.out.reads)

// 			flagstat_results_ch = flagstats.out.flagstats
// 				.concat(flagstats.out.counts)
// 				.concat(flagstats.out.is_paired)
// 				.map { sample, files -> files }
// 			results_ch = results_ch.concat(flagstat_results_ch)

// 		}

// 		if (run_bam_analysis) {

// 			bam_analysis(fq2bam.out.reads)
// 			results_ch = results_ch.concat(bam_analysis.out.results)

// 		}

//     }

// 	if (run_fastq_analysis) {

// 		fastq_analysis(preprocessed_ch)
// 		results_ch = results_ch.concat(fastq_analysis.out.results)

// 	}

// 	if (run_amplicon_analysis) {

// 		amplicon_analysis(preprocessed_ch)
// 		results_ch = results_ch.concat(amplicon_analysis.out.results)

// 	}

// 	if (!params.skip_collate) {
// 		collate_results(
// 			results_ch.collect(),
// 			"${projectDir}/scripts/ExtractProfiledCounts_210823.R",
// 			params.GTDB_markers
// 		)
// 	}

// }
