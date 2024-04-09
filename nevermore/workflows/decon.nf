include { remove_host_kraken2_individual; remove_host_kraken2 } from "../modules/decon/kraken2"
include { sortmerna } from "../modules/decon/sortmerna"


workflow nevermore_decon {

		take:
			fastq_ch
		
		main:
			preprocessed_ch = Channel.empty()

			if (params.run_sortmerna) {

				fastq_ch
					.branch {
						metaT: it[0].containsKey("library_source") && it[0].library_source == "metaT"
						metaG: true
					}
					.set { for_sortmerna_ch }

				sortmerna(for_sortmerna_ch.metaT, params.sortmerna_db)
				preprocessed_ch = for_sortmerna_ch.metaG
					.concat(sortmerna.out.fastqs)

			} else {

				preprocessed_ch = fastq_ch

			}
	
			if (params.remove_host) {
	
				remove_host_kraken2_individual(preprocessed_ch, params.remove_host_kraken2_db)
	
				preprocessed_ch = remove_host_kraken2_individual.out.reads				
	
			}

		emit:
			reads = preprocessed_ch

}
