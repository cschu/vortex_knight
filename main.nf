#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { bam2fq; fq2bam; prepare_fastqs } from "./modules/vknight/convert"
include { nevermore_simple_preprocessing } from "./workflows/nevermore/nevermore"
include { amplicon_analysis; bam_analysis; fastq_analysis } from "./workflows/vknight/vknight"
include { classify_sample } from "./modules/nevermore/functions"
include { remove_host_kraken2 } from "./modules/nevermore/decon/kraken2"
include { flagstats; count_reads_flagstats } from "./modules/vknight/stats"


def run_kraken2 = (!params.skip_kraken2 || params.run_kraken2) && !params.amplicon_seq;
def run_mtags = (!params.skip_mtags || params.run_mtags);
def run_mapseq = (run_mtags && (!params.skip_mapseq || params.run_mapseq) && params.mapseq_bin)
def run_motus2 = (!params.skip_motus2 || params.run_motus2) && !params.amplicon_seq;
def run_pathseq = (!params.skip_pathseq || params.run_pathseq) && !params.amplicon_seq;
def run_read_counter = (!params.skip_read_counter || params.run_read_counter)

def get_basecounts = (!params.skip_basecounts || params.run_basecounts);
def convert_fastq2bam = (run_pathseq || get_basecounts);

def do_preprocessing = (!params.skip_preprocessing || params.run_preprocessing)

def run_bam_analysis = run_pathseq && !params.amplicon_seq
def run_fastq_analysis = (run_kraken2 || run_mtags || run_mapseq || run_motus2 || run_read_counter) && !params.amplicon_seq
def run_amplicon_analysis = params.amplicon_seq


workflow {

	fastq_ch = Channel
		.fromPath(params.input_dir + "/" + "**.{fastq,fq,fastq.gz,fq.gz}")
		.map { file ->
				def sample = file.name.replaceAll(/.(fastq|fq)(.gz)?$/, "")
				sample = sample.replaceAll(/_R?[12]$/, "")
				return tuple(sample, file)
		}
		.groupTuple(sort: true)
        .map { classify_sample(it[0], it[1]) }

	bam_ch = Channel
		.fromPath(params.input_dir + "/" + "**.bam")
		.map { file ->
			def sample = file.name.replaceAll(/.bam$/, "")
			return tuple(sample, file)
		}
		.groupTuple(sort: true)

	bam2fq(bam_ch)

	bfastq_ch = bam2fq.out.reads
		.map { classify_sample(it[0], it[1]) }

	prepare_fastqs(fastq_ch)


	if (do_preprocessing) {

		raw_fastq_ch = prepare_fastqs.out.reads.concat(bfastq_ch)

		nevermore_simple_preprocessing(raw_fastq_ch)


		if (params.remove_host) {

			remove_host_kraken2(nevermore_simple_preprocessing.out.main_reads_out, params.remove_host_kraken2_db)

			preprocessed_ch = remove_host_kraken2.out.reads

		} else {

			preprocessed_ch = nevermore_simple_preprocessing.out.main_reads_out

		}

	} else {

		preprocessed_ch = prepare_fastqs.out.reads
			.concat(bfastq_ch)

	}


	if (get_basecounts || run_bam_analysis) {

		fq2bam(preprocessed_ch)

		if (get_basecounts) {

	        flagstats(fq2bam.out.reads)

    	    count_reads_flagstats(flagstats.out.flagstats)

		}

		if (run_bam_analysis) {

			bam_analysis(fq2bam.out.reads)

		}

    }


	if (run_fastq_analysis) {

		fastq_analysis(preprocessed_ch)

	}

	if (run_amplicon_analysis) {

		amplicon_analysis(preprocessed_ch)

	}

}
