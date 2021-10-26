process kraken2 {
    publishDir params.output_dir, mode: params.publish_mode

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("${sample.id}/${sample.id}.kraken2_report.txt"), emit: kraken2_out

    script:
    def is_paired = (sample.is_paired) ? "--paired" : "";
    """
    mkdir -p ${sample.id}
    kraken2 --db ${params.kraken_database} --threads $task.cpus --gzip-compressed --report ${sample.id}/${sample.id}.kraken2_report.txt ${is_paired} \$(ls $reads)
    """
}
