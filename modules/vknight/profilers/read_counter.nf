process read_counter {
    publishDir params.output_dir, mode: params.publish_mode

    input:
    tuple val(sample), path(reads)
	path(marker_db)

    output:
    tuple val(sample), path("${sample.id}/${sample.id}.read_counter.txt"), emit: read_counter_out

    script:
    def read_counter_input = (sample.is_paired) ? "-f ${sample.id}_R1.fastq.gz -r ${sample.id}_R2.fastq.gz" : "-s ${sample.id}_R1.fastq.gz";
    """
    mkdir -p ${sample.id}
	read_counter map -t $task.cpus -db {marker_db} -l {params.read_counter_min_length} ${read_counter_input} > ${sample.id}/${sample.id}.read_counter.txt
    """
}

