process kraken2 {
    publishDir params.output_dir, mode: "copy"
    container "registry.git.embl.org/schudoma/kraken2-docker:latest"
	label 'kraken2'
	label "large"

    input:
    tuple val(sample), path(fastqs)
	path(kraken_db)

    output:
    tuple val(sample), path("${sample.id}/${sample.id}.kraken2_report.txt"), emit: kraken2_out

    script:

    def r1_files = fastqs.findAll( { it.name.endsWith("_R1.fastq.gz") && !it.name.matches("(.*)(singles|orphans|chimeras)(.*)") } )
	def r2_files = fastqs.findAll( { it.name.endsWith("_R2.fastq.gz") } )
	def orphans = fastqs.findAll( { it.name.matches("(.*)(singles|orphans|chimeras)(.*)") } )

    def is_paired = (r1_files.size() != 0 && r2_files.size() != 0) ? "--paired" : ""
	// if (r1_files.size() != 0 && r2_files.size() != 0) {
	// 	// input_files += "-f ${r1_files.join(' ')} -r ${r2_files.join(' ')}"
    //     is_paired += "--paired"
	// } 
    // else if (r1_files.size() != 0) {
	// 	input_files += "-s ${r1_files.join(' ')}"
	// } else if (r2_files.size() != 0) {
	// 	input_files += "-s ${r2_files.join(' ')}"
	// }
    
    // if (orphans.size() != 0) {
	// 	input_files += " -s ${orphans.join(' ')}"
	// }



    // def is_paired = (sample.is_paired) ? "--paired" : "";
    def kraken_params = "--use-mpa-style --gzip-compressed ${is_paired}"



    """
    mkdir -p ${sample.id}
    kraken2 --db ${kraken_db} --threads ${task.cpus} --minimum-hit-groups ${params.kraken2_min_hit_groups} --report ${sample.id}/${sample.id}.kraken2_report.txt ${kraken_params} ${fastqs}
    """
}
