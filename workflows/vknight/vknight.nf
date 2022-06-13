#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { kraken2 } from  "../../modules/vknight/profilers/kraken2"
include { mtags_extract; mtags_annotate; mtags_merge } from "../../modules/vknight/profilers/mtags"
include { motus2 } from  "../../modules/vknight/profilers/motus2"
include { mapseq; mapseq_with_customdb; collate_mapseq_tables } from "../../modules/vknight/profilers/mapseq"
include { pathseq } from "../../modules/vknight/profilers/pathseq"
include { read_counter } from "../../modules/vknight/profilers/read_counter"


if (!params.publish_mode) {
	params.publish_mode = "symlink"
}

if (!params.output_dir) {
	params.output_dir = "vknight_out"
}

output_dir = "${params.output_dir}"

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

		if (run_motus2) {
			motus2(fastq_ch, params.motus_database)
			out_ch = out_ch.concat(motus2.out.motus_out)
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

		mtags_extract(fastq_ch)

		mapseq_ch = Channel.empty()

		if (params.mapseq_db) {

			mapseq_with_customdb(mtags_extract.out.mtags_out, params.mapseq_db)
			mapseq_ch = mapseq_with_customdb.out.bac_ssu.collect() 

		} else {

			mapseq(mtags_extract.out.mtags_out)
			mapseq_ch = mapseq.out.bac_ssu.collect() 

		}

		collate_mapseq_tables(mapseq_ch)
		out_ch = out_ch.concat(collate_mapseq_tables.out.ssu_tables)


	emit:
		results = out_ch
}


