params.pathseq_skip_host_alignment = false
params.pathseq_skip_quality_filters = false


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

	def microbe_seq = ""
	if (params.pathseq_db_microbe_fasta) {
		microbe_seq = "--microbe-fasta ${params.pathseq_db_microbe_fasta}"
	} else {
		microbe_seq = "--microbe-dict ${params.pathseq_db_microbe_dict}"
	}

    """
    mkdir -p ${sample.id}

	gatk --java-options \"-Xmx${maxmem}g\" PathSeqPipelineSpark \\
		--input ${bam} \\
		--filter-bwa-image ${params.pathseq_db_filter_bwa_image} \\
		--kmer-file ${params.pathseq_db_kmer_file} \\
		--min-clipped-read-length ${params.pathseq_min_clipped_read_length} \\
		${microbe_seq} \\
		--microbe-bwa-image ${params.pathseq_db_microbe_bwa_image} \\
		--taxonomy-file ${params.pathseq_db_taxonomy_file} \\
		--output ${sample.id}/${sample.id}.pathseq.bam \\
		--scores-output ${sample.id}/${sample.id}.pathseq.scores \\
		--score-metrics ${sample.id}/${sample.id}.pathseq.score_metrics \\
		--is-host-aligned ${params.pathseq_skip_host_alignment} \\
		--skip-quality-filters ${params.pathseq_skip_quality_filters}
    """
}
