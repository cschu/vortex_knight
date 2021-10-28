process kraken2 {
    publishDir params.output_dir, mode: params.publish_mode

    input:
    tuple val(sample), path(reads)
	path(kraken_db)

    output:
    tuple val(sample), path("${sample.id}/${sample.id}.kraken2_report.txt"), emit: kraken2_out

    script:
    def is_paired = (sample.is_paired) ? "--paired" : "";
    """
    mkdir -p ${sample.id}
    kraken2 --db ${kraken_db} --threads $task.cpus --minimum-hit-groups {params.kraken2_min_hit_groups} --gzip-compressed --report ${sample.id}/${sample.id}.kraken2_report.txt ${is_paired} \$(ls $reads)
    """
}
