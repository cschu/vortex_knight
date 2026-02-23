params.pathseq_skip_host_alignment = false
params.pathseq_skip_quality_filters = false
params.pathseq_filter_duplicates = true


process pathseq {
    publishDir params.output_dir, mode: "copy", pattern: "**.scores"
	container "broadinstitute/gatk:4.1.6.0"
	// container "quay.io/biocontainers/gatk4:4.6.2.0--py310hdfd78af_1"
	label "highmemXLarge"
	tag "${sample.id}"

    input:
    tuple val(sample), path(bam)
	path(pathseq_db)

    output:
    tuple val(sample), path("pathseq/${sample.id}/${sample.id}.pathseq.score*"), emit: scores
    tuple val(sample), path("pathseq/${sample.id}/${sample.id}.pathseq.bam*"), emit: bam
	tuple val(sample), path("pathseq/${sample.id}/${sample.id}.pathseq.filter_metrics"), emit: filter_metrics

    script:
    def maxmem = task.memory.toGiga()
	def filter_duplicates = params.pathseq_filter_duplicates == true
	def microbe_seq_input = (params.new_pathseq) 
		? "--microbe-dict ${params.pathseq_db}/pathseq_microbe.dict"
		: "--microbe-fasta ${pathseq_db}/pathseq_microbe.fa"

	"""
    mkdir -p ${sample.id}

	gatk --java-options \"-Xmx${maxmem}g\" PathSeqPipelineSpark \\
		--conf spark.hadoop.hadoop.security.authentication=simple \\
		--input ${bam} \\
		--filter-bwa-image ${pathseq_db}/pathseq_host.fa.img \\
		--kmer-file ${pathseq_db}/pathseq_host.bfi \\
		${microbe_seq_input} \\
		--microbe-bwa-image ${pathseq_db}/pathseq_microbe.fa.img \\
		--taxonomy-file ${pathseq_db}/pathseq_taxonomy.db \\
		--min-clipped-read-length ${params.pathseq_min_clipped_read_length} \\
		--is-host-aligned ${params.pathseq_skip_host_alignment} \\
		--skip-quality-filters ${params.pathseq_skip_quality_filters} \\
		--filter-duplicates ${filter_duplicates} \\
		--output pathseq/${sample.id}/${sample.id}.pathseq.bam \\
		--scores-output pathseq/${sample.id}/${sample.id}.pathseq.scores \\
		--score-metrics pathseq/${sample.id}/${sample.id}.pathseq.score_metrics \\
		--filter-metrics pathseq/${sample.id}/${sample.id}.pathseq.filter_metrics
    """
}
