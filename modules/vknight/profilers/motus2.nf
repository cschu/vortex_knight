process motus2 {
    publishDir params.output_dir, mode: params.publish_mode

    input:
    tuple val(sample), path(reads)
	path(motus_db)

    output:
    tuple val(sample), path("${sample.id}/${sample.id}.motus.txt"), emit: motus_out

    script:
    def motus_input = (sample.is_paired) ? "-f ${sample.id}_R1.fastq.gz -r ${sample.id}_R2.fastq.gz" : "-s ${sample.id}_R1.fastq.gz";
    """
    mkdir -p ${sample.id}
    motus profile -t $task.cpus -k ${params.motus2_tax_level} -c -v 7 -l ${params.motus2_min_length} -g ${params.motus2_n_marker_genes} -db ${motus_db} ${motus_input} > ${sample.id}/${sample.id}.motus.txt
    """
}
