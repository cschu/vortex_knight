#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { kraken2 } from  "../modules/profilers/kraken2"
include { mtags_extract; mtags_annotate; mtags_merge } from "../modules/profilers/mtags"
include { motus } from  "../modules/profilers/motus"
include { mapseq; mapseq_with_customdb; collate_mapseq_tables } from "../modules/profilers/mapseq"
include { pathseq } from "../modules/profilers/pathseq"
include { read_counter } from "../modules/profilers/read_counter"
include { idtaxa } from "../modules/profilers/idtaxa"
include { run_metaphlan4 as metaphlan4; collate_metaphlan4_tables } from "../../nevermore/modules/profilers/metaphlan4"

include { fq2fa } from "../../nevermore/modules/converters/fq2fa"
include { fastqc } from "../../nevermore/modules/qc/fastqc"
include { multiqc } from "../../nevermore/modules/qc/multiqc"
include { fq2bam } from "../../nevermore/modules/converters/fq2bam"
include { nevermore_simple_preprocessing } from "../../nevermore/workflows/nevermore"
include { remove_host_kraken2; remove_host_kraken2_individual } from "../../nevermore/modules/decon/kraken2"
include { flagstats } from "../../nevermore/modules/stats"
include { collate_results } from "../modules/collate"
include { collate_stats } from "../../nevermore/modules/collate"
include { starmap } from "../modules/decon/starmap"

include { vk_decon } from "./decon"


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

params.kraken2_db = params.kraken_database
params.motus_db = params.motus_database
params.read_counter_db = params.read_counter_database
params.pathseq_db = params.pathseq_database



def run_kraken2 = (!params.skip_kraken2 || params.run_kraken2);
def run_mtags = (!params.skip_mtags || params.run_mtags || params.run_downstream_mtags);
def run_mapseq = (run_mtags && (!params.skip_mapseq || params.run_mapseq) && params.mapseq_bin)
def run_motus = (!params.skip_motus || params.run_motus);
def run_pathseq = (!params.skip_pathseq || params.run_pathseq);
def run_read_counter = (!params.skip_read_counter || params.run_read_counter)
def run_idtaxa = (!params.skip_idtaxa || params.run_idtaxa)
def run_metaphlan4 = (!params.skip_metaphlan4 || params.run_metaphlan4)

def get_basecounts = (!params.skip_basecounts || params.run_basecounts);

def run_bam_analysis = run_pathseq && !params.amplicon_seq
def run_fastq_analysis = (run_kraken2 || run_mtags || run_mapseq || run_motus || run_read_counter) && !params.amplicon_seq
def run_amplicon_analysis = params.amplicon_seq



workflow bam_analysis {
	take:
		bam_ch

	main:
		out_ch = Channel.empty()
    	if (run_pathseq) {
			pathseq(bam_ch, params.pathseq_db)
			out_ch = out_ch.mix(pathseq.out.scores)
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

		if (run_metaphlan4) {
			metaphlan4(fastq_ch, params.metaphlan4_db)
			collate_metaphlan4_tables(metaphlan4.out.mp4_table)
		}

		if (run_kraken2) {
			kraken2(fastq_ch, params.kraken2_db)
			out_ch = out_ch.mix(kraken2.out.kraken2_out)
		}

		if (run_motus) {
			motus(fastq_ch, params.motus_db)
			out_ch = out_ch.mix(motus.out.motus_out)
		}

		if (run_read_counter) {
			read_counter(fastq_ch, params.read_counter_db)
			out_ch = out_ch.mix(read_counter.out.read_counter_out)
		}

		out_ch = out_ch
			.map { sample, files -> return files }

		if (run_mtags || run_mapseq) {
			mtags_extract(fastq_ch)
	
			if (run_mtags) {
				mtags_annotate(mtags_extract.out.mtags_out)
	
				mtags_merge(mtags_annotate.out.mtags_bins.collect())

				out_ch = out_ch.mix(mtags_merge.out.mtags_tables)
			}

			if (run_mapseq) {

				mapseq_ch = Channel.empty()

				if (params.mapseq_db) {

					mapseq_with_customdb(mtags_extract.out.mtags_out, params.mapseq_db)
					mapseq_ch = mapseq_with_customdb.out.bac_ssu.collect()
					out_ch = out_ch.mix(mapseq_with_customdb.out.bac_ssu)

				} else {

					mapseq(mtags_extract.out.mtags_out)
					mapseq_ch = mapseq.out.bac_ssu.collect()
					out_ch = out_ch.mix(mapseq.out.bac_ssu)

				}
	
				collate_mapseq_tables(mapseq_ch)
				out_ch = out_ch.mix(collate_mapseq_tables.out.ssu_tables)
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

		if (run_mapseq && params.mapseq_db) {

			mapseq_with_customdb(fq2fa.out.reads, params.mapseq_db)
			mapseq_ch = mapseq_with_customdb.out.bac_ssu.collect() 
			out_ch = out_ch.mix(mapseq_with_customdb.out.bac_ssu)

		} else if (run_mapseq) {

			mapseq(fq2fa.out.reads)
			mapseq_ch = mapseq.out.bac_ssu.collect()
			out_ch = out_ch.mix(mapseq.out.bac_ssu)

		} 
		
		// out_ch = Channel.empty()
		if (run_idtaxa) {
			idtaxa(fq2fa.out.reads, params.idtaxa_classifier_db)
			out_ch = out_ch.mix(idtaxa.out.count_table)
		}

		if (run_mapseq) {
			collate_mapseq_tables(mapseq_ch)
			out_ch = out_ch.mix(collate_mapseq_tables.out.ssu_tables)
		} 


	emit:
		results = out_ch
}


workflow vknight_main {
	take:
		fastq_ch
	main:
		results_ch = Channel.empty()
		preprocessed_ch = fastq_ch
		
		if (params.remove_host == "vk_decon") {

			vk_decon(preprocessed_ch)
			preprocessed_ch = vk_decon.out.reads

		}

		if (get_basecounts || run_bam_analysis) {

			fq2bam(preprocessed_ch)

			if (get_basecounts) {

				flagstats(fq2bam.out.reads, "basecounts")

				flagstat_results_ch = flagstats.out.flagstats
					.mix(flagstats.out.counts)
					.mix(flagstats.out.is_paired)
					.map { sample, files -> files }
				results_ch = results_ch.mix(flagstat_results_ch)

				flagstats_libtype(
					flagstats.out.is_paired
						.map { id, file -> file }
						.collate()
				)

			}

			if (run_bam_analysis) {

				bam_analysis(fq2bam.out.reads)
				results_ch = results_ch.mix(bam_analysis.out.results)

			}

		}

		if (run_fastq_analysis) {

			fastq_analysis(preprocessed_ch)
			results_ch = results_ch.mix(fastq_analysis.out.results)

		}

		if (run_amplicon_analysis) {

			amplicon_analysis(preprocessed_ch)
			results_ch = results_ch.mix(amplicon_analysis.out.results)

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
