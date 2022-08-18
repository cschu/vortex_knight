#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { kraken2 } from  "../modules/profilers/kraken2"
include { mtags_extract; mtags_annotate; mtags_merge } from "../modules/profilers/mtags"
include { motus } from  "../modules/profilers/motus"
include { mapseq; mapseq_with_customdb; collate_mapseq_tables } from "../modules/profilers/mapseq"
include { pathseq } from "../modules/profilers/pathseq"
include { read_counter } from "../modules/profilers/read_counter"
include { fq2fa } from "../../nevermore/modules/converters/fq2fa"


if (!params.publish_mode) {
	params.publish_mode = "symlink"
}

if (!params.output_dir) {
	params.output_dir = "vknight_out"
}

output_dir = "${params.output_dir}"

if (!params.motus_min_length) {
	params.motus_min_length = 30
}

if (!params.motus_n_marker_genes) {
	params.motus_n_marker_genes = 1
}

if (!params.pathseq_min_clipped_read_length) {
	params.pathseq_min_clipped_read_length = 31
}

def run_kraken2 = (!params.skip_kraken2 || params.run_kraken2);
def run_mtags = (!params.skip_mtags || params.run_mtags);
def run_mapseq = (run_mtags && (!params.skip_mapseq || params.run_mapseq) && params.mapseq_bin)
def run_motus = (!params.skip_motus || params.run_motus);
def run_pathseq = (!params.skip_pathseq || params.run_pathseq);
def run_read_counter = (!params.skip_read_counter || params.run_read_counter)


workflow bam_analysis {
	take:
		bam_ch

	main:
		out_ch = Channel.empty()
    	if (run_pathseq) {
			pathseq(bam_ch, params.pathseq_database)
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
			kraken2(fastq_ch, params.kraken_database)
			out_ch = out_ch.concat(kraken2.out.kraken2_out)
		}

		if (run_motus) {
			motus(fastq_ch, params.motus_database)
			out_ch = out_ch.concat(motus.out.motus_out)
		}

		if (run_read_counter) {
			read_counter(fastq_ch, params.read_counter_database)
			out_ch = out_ch.concat(read_counter.out.read_counter_out)
		}

		out_ch = out_ch
			.map { sample, files -> return files }

		if (run_mtags) {
			mtags_extract(fastq_ch)
	
			mtags_annotate(mtags_extract.out.mtags_out)
	
			mtags_merge(mtags_annotate.out.mtags_bins.collect())

			out_ch = out_ch.concat(mtags_merge.out.mtags_tables)

			if (run_mapseq) {

				mapseq_ch = Channel.empty()

				if (params.mapseq_db) {

					mapseq_with_customdb(mtags_extract.out.mtags_out, params.mapseq_db)
					mapseq_ch = mapseq_with_customdb.out.bac_ssu.collect()
					out_ch = out_ch.concat(mapseq_with_customdb.out.bac_ssu)

				} else {

					mapseq(mtags_extract.out.mtags_out)
					mapseq_ch = mapseq.out.bac_ssu.collect()
					out_ch = out_ch.concat(mapseq.out.bac_ssu)

				}
	
				collate_mapseq_tables(mapseq_ch)
				out_ch = out_ch.concat(collate_mapseq_tables.out.ssu_tables)
			}
		}

	emit:
		results = out_ch
}


workflow amplicon_analysis {
	take:
		fastq_ch

	main:
		out_ch = Channel.empty()

		fq2fa(fastq_ch)

		mapseq_ch = Channel.empty()

		if (params.mapseq_db) {

			mapseq_with_customdb(fq2fa.out.reads, params.mapseq_db)
			mapseq_ch = mapseq_with_customdb.out.bac_ssu.collect() 
			out_ch = out_ch.concat(mapseq_with_customdb.out.bac_ssu)

		} else {

			mapseq(fq2fa.out.reads)
			mapseq_ch = mapseq.out.bac_ssu.collect()
			out_ch = out_ch.concat(mapseq.out.bac_ssu)

		}

		collate_mapseq_tables(mapseq_ch)
		out_ch = out_ch.concat(collate_mapseq_tables.out.ssu_tables)


	emit:
		results = out_ch
}


workflow vknight_main {
	take:
		fastq_ch
	main:
		results_ch = Channel.empty()

		if (do_preprocessing) {

			nevermore_simple_preprocessing(fastq_ch)

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

