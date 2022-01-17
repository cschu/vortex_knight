process pathseq {
    publishDir params.output_dir, mode: params.publish_mode

    input:
    tuple val(sample), path(bam)
	path(pathseq_db)

    output:
    tuple val(sample), path("${sample.id}/${sample.id}.pathseq.score*"), emit: scores
    tuple val(sample), path("${sample.id}/${sample.id}.pathseq.bam*"), emit: bam

    script:
    def maxmem = task.memory.toGiga()
    """
    mkdir -p ${sample.id}
    gatk --java-options \"-Xmx${maxmem}g\" PathSeqPipelineSpark \\
        --input $bam \\
        --filter-bwa-image ${pathseq_db}/reference.fasta.img \\
        --kmer-file ${pathseq_db}/host.hss \\
        --min-clipped-read-length ${params.pathseq_min_clipped_read_length} \\
        --microbe-fasta ${pathseq_db}/microbe.fasta \\
        --microbe-bwa-image ${pathseq_db}/microbe.fasta.img \\
        --taxonomy-file ${pathseq_db}/microbe.db \\
        --output ${sample.id}/${sample.id}.pathseq.bam \\
        --scores-output ${sample.id}/${sample.id}.pathseq.scores \\
        --score-metrics ${sample.id}/${sample.id}.pathseq.score_metrics
    """
}
