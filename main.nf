#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { bam2fq; fq2bam; prepare_fastqs } from "./modules/vknight/convert"
include { nevermore_simple_preprocessing } from "./workflows/nevermore/nevermore"
include { amplicon_analysis; bam_analysis; fastq_analysis; collate_data } from "./workflows/vknight/vknight"
include { classify_sample } from "./modules/nevermore/functions"
include { remove_host_kraken2 } from "./modules/nevermore/decon/kraken2"
include { flagstats; count_reads } from "./modules/vknight/stats"


def run_kraken2 = (!params.skip_kraken2 || params.run_kraken2) && !params.amplicon_seq;
def run_mtags = (!params.skip_mtags || params.run_mtags);
def run_mapseq = (run_mtags && (!params.skip_mapseq || params.run_mapseq) && params.mapseq_bin)
def run_motus2 = (!params.skip_motus2 || params.run_motus2)  && !params.amplicon_seq;
def run_pathseq = (!params.skip_pathseq || params.run_pathseq) && !params.amplicon_seq;
def run_count_reads = (!params.skip_counts || params.run_counts) &&!params.amplicon_seq;
def convert_fastq2bam = (run_pathseq || run_count_reads);

def do_preprocessing = (!params.skip_preprocessing || params.run_preprocessing)

def run_bam_analysis = run_pathseq && !params.amplicon_seq
def run_fastq_analysis = (run_kraken2 || run_mtags || run_mapseq || run_motus2) && !params.amplicon_seq
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
        .map { classify_sample(it[0], it[1]) }

	bam2fq(bam_ch)

	if (do_preprocessing) {

		prepare_fastqs(fastq_ch)

		raw_fastq_ch = prepare_fastqs.out.reads.concat(bam2fq.out.reads)

		nevermore_simple_preprocessing(raw_fastq_ch)


		if (params.remove_host) {

			remove_host_kraken2(nevermore_simple_preprocessing.out.main_reads_out, params.remove_host_kraken2_db)

			preprocessed_ch = remove_host_kraken2.out.reads

		} else {

			preprocessed_ch = nevermore_simple_preprocessing.out.main_reads_out

		}

	} else {

		preprocessed_ch = fastq_ch

	}


	if (run_count_reads || run_bam_analysis) {

		fq2bam(preprocessed_ch)

		if (run_count_reads) {

	        flagstats(fq2bam.out.reads)

    	    count_reads(flagstats.out.flagstats)

		}

		if (run_bam_analysis) {

			bam_analysis(fq2bam.out.reads)

		}

    }


	if (run_fastq_analysis) {

		preprocessed_ch.view()
		fastq_analysis(preprocessed_ch)

	}

	if (run_amplicon_analysis) {

		amplicon_analysis(preprocessed_ch)

	}


//	/* collate data */
//
//	if (params.collate_script != null && params.collate_script != "") {
//		data_to_collate_ch = Channel.empty()
//
//		if (run_kraken2) {
//			data_to_collate_ch = data_to_collate_ch.concat(kraken2.out.kraken2_out)
//		}
//
//		if (run_count_reads) {
//			data_to_collate_ch = data_to_collate_ch.concat(count_reads.out.counts)
//				.concat(count_reads.out.is_paired)
//		}
//
//		if (run_motus2) {
//			data_to_collate_ch = data_to_collate_ch.concat(motus2.out.motus_out)
//		}
//
//		if (run_pathseq) {
//			data_to_collate_ch = data_to_collate_ch.concat(pathseq.out.scores)
//		}
//
//		data_to_collate_ch = data_to_collate_ch
//			.map { sample, files -> return files }
//
//		if (run_mtags) {
//			data_to_collate_ch = data_to_collate_ch.concat(mtags_merge.out.mtags_tables)
//		}
//
//		collate_data(data_to_collate_ch.collect())
//	}
}
