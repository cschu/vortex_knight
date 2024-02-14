process collate_results {
	publishDir params.output_dir, mode: params.publish_mode

	input:
	path(results)
	path(collate_script)
	path(gtdb_markers)

	output:
	path("collated/*.rds"), emit: collated, optional: true

	script:
	"""
	mkdir -p collated/

	mkdir -p kraken2/
	(mv *kraken2_report.txt kraken2/) || :

	mkdir -p motus/
	(mv *motus.txt motus/) || :

	mkdir -p pathseq/
	(mv *pathseq.scores pathseq/) || :

	mkdir -p libsize/
	(mv *libsize.txt libsize/) || :

	mkdir -p liblayout/
	(mv *is_paired.txt liblayout/) || :

	mkdir -p flagstats/
	(mv *flagstats.txt flagstats/) || :

	mkdir -p mapseq/
	(mv *mseq mapseq/) || :

	mkdir -p mtags_tables/
	(mv merged_profile.genus.tsv mtags_tables/) || :

	mkdir -p read_counter/
	(mv *read_counter.txt read_counter/) || :
	
	mkdir -p idtaxa/
	(mv *IDTaxa.tsv idtaxa/) || :

	mkdir -p raw_counts/
	(mv *.txt raw_counts/) || :



	Rscript --vanilla ${collate_script} \
		--libdir \$(dirname \$(readlink ${collate_script})) \
		--gtdb_markers ${gtdb_markers} \
		--kraken2_res_path kraken2/ \
		--mOTUs_res_path motus/ \
		--PathSeq_res_path pathseq/ \
		--mTAGs_res_path mtags_tables/ \
		--mapseq_res_path mapseq/ \
		--libsize_res_path libsize/ \
		--lib_layout_res_path liblayout/ \
		--N_raw_counts_path raw_counts/ \
		--read_counter_res_path read_counter/ \
		--IDtaxa_res_path idtaxa/ \
		--out_folder collated/ 
	"""
}