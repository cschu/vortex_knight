params.pathseq_skip_host_alignment = false
params.pathseq_skip_quality_filters = false
params.pathseq_filter_duplicates = true


process pathseq {
    publishDir params.output_dir, mode: "copy"

    input:
    tuple val(sample), path(bam)
	path(pathseq_db)

    output:
    tuple val(sample), path("${sample.id}/${sample.id}.pathseq.score*"), emit: scores
    tuple val(sample), path("${sample.id}/${sample.id}.pathseq.bam*"), emit: bam
	tuple val(sample), path("${sample.id}/${sample.id}.pathseq.filter_metrics"), emit: filter_metrics

    script:
    def maxmem = task.memory.toGiga()

	def microbe_seq = ""
	if (params.pathseq_db_microbe_fasta) {
		microbe_seq = "--microbe-fasta ${params.pathseq_db_microbe_fasta}"
	} else {
		microbe_seq = "--microbe-dict ${params.pathseq_db_microbe_dict}"
	}

	def filter_duplicates = params.pathseq_filter_duplicates == true
	// if (params.pathseq_filter_duplicates) {
	// 	filter_duplicates += "--filter-duplicates true"
	// } else {
	// 	filter_duplicates += "--filter-duplicates false"
	// }

    """
    mkdir -p ${sample.id}

	gatk --java-options \"-Xmx${maxmem}g\" PathSeqPipelineSpark \\
		--input ${bam} \\
		--filter-bwa-image ${pathseq_db}/pathseq_host.fa.img \\
		--kmer-file ${pathseq_db}/pathseq_host.bfi \\
		--min-clipped-read-length ${params.pathseq_min_clipped_read_length} \\
		${microbe_seq} \\
		--microbe-bwa-image ${params.pathseq_db}/pathseq_microbe.fa.img \\
		--taxonomy-file ${params.pathseq_db_taxonomy_file} \\
		--output ${sample.id}/${sample.id}.pathseq.bam \\
		--scores-output ${sample.id}/${sample.id}.pathseq.scores \\
		--score-metrics ${sample.id}/${sample.id}.pathseq.score_metrics \\
		--is-host-aligned ${params.pathseq_skip_host_alignment} \\
		--skip-quality-filters ${params.pathseq_skip_quality_filters} \\
		--filter-metrics ${sample.id}/${sample.id}.pathseq.filter_metrics \\
		--filter-duplicates ${filter_duplicates}
    """
}
