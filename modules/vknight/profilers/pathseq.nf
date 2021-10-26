process pathseq {
    publishDir params.output_dir, mode: params.publish_mode

    input:
    tuple val(sample), path(bam)

    output:
    tuple val(sample), path("${sample.id}/${sample.id}.pathseq.score*"), emit: scores
    tuple val(sample), path("${sample.id}/${sample.id}.pathseq.bam*"), emit: bam

    script:
    """
    mkdir -p ${sample.id}
    maxmem=\$(echo \"$task.memory\"| sed 's/ GB/g/g')
    gatk --java-options \"-Xmx\$maxmem\" PathSeqPipelineSpark \\
        --input $bam \\
        --filter-bwa-image ${params.pathseq_database}/reference.fasta.img \\
        --kmer-file ${params.pathseq_database}/host.hss \\
        --min-clipped-read-length ${params.pathseq_min_clipped_read_length} \\
        --microbe-fasta ${params.pathseq_database}/microbe.fasta \\
        --microbe-bwa-image ${params.pathseq_database}/microbe.fasta.img \\
        --taxonomy-file ${params.pathseq_database}/microbe.db \\
        --output ${sample.id}/${sample.id}.pathseq.bam \\
        --scores-output ${sample.id}/${sample.id}.pathseq.scores \\
        --score-metrics ${sample.id}/${sample.id}.pathseq.score_metrics
    """
}
