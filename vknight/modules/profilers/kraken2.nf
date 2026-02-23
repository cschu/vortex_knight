process kraken2 {
    publishDir params.output_dir, mode: "copy", pattern: "**.kraken2_report.txt"
    container "quay.io/biocontainers/kraken2:2.17.1--pl5321h077b44d_0"
	label 'kraken2'
	label "highmemMedium"
    tag "${sample.id}"

    input:
    tuple val(sample), path(fastqs)
	path(kraken_db)

    output:
    tuple val(sample), path("kraken2/${sample.id}/${sample.id}.kraken2_report.txt"), emit: kraken2_out

    script:
    def r1_files = fastqs.findAll( { it.name.endsWith("_R1.fastq.gz") && !it.name.matches("(.*)(singles|orphans|chimeras)(.*)") } )
	def r2_files = fastqs.findAll( { it.name.endsWith("_R2.fastq.gz") } )
	def orphans = fastqs.findAll( { it.name.matches("(.*)(singles|orphans|chimeras)(.*)") } )
    def is_paired = (r1_files.size() != 0 && r2_files.size() != 0) ? "--paired" : ""
    def kraken_params = "--use-mpa-style --gzip-compressed ${is_paired}"

    """
    mkdir -p kraken2/${sample.id}/
    kraken2 --db ${kraken_db} --threads ${task.cpus} --minimum-hit-groups ${params.kraken2_min_hit_groups} --report kraken2/${sample.id}/${sample.id}.kraken2_report.txt ${kraken_params} ${fastqs}
    """
}
