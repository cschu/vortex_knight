params.pathseq_skip_host_alignment = false
params.pathseq_skip_quality_filters = false
params.pathseq_filter_duplicates = true


process pathseq {
    publishDir params.output_dir, mode: "copy"
	// biocontainers docker image causes issues on EMBL HPC
	// https://gatk.broadinstitute.org/hc/en-us/community/posts/360078378372--Could-not-determine-local-host-name
	// needs to be gatk4!
	// container "quay.io/biocontainers/gatk:3.8--py36_4"
	// container "registry.git.embl.org/schudoma/pathseq-docker:latest"
	container "registry.git.embl.org/schudoma/pathseq-docker:with_user"
	label "highmemXLarge"

    input:
    tuple val(sample), path(bam)
	path(pathseq_db)

    output:
    tuple val(sample), path("${sample.id}/${sample.id}.pathseq.score*"), emit: scores
    tuple val(sample), path("${sample.id}/${sample.id}.pathseq.bam*"), emit: bam
	tuple val(sample), path("${sample.id}/${sample.id}.pathseq.filter_metrics"), emit: filter_metrics

    script:
    def maxmem = task.memory.toGiga()

	// def microbe_seq = ""
	// if (params.pathseq_db_microbe_fasta) {
	// 	microbe_seq = "--microbe-fasta ${params.pathseq_db_microbe_fasta}"
	// } else {
	// 	microbe_seq = "--microbe-dict ${params.pathseq_db_microbe_dict}"
	// }

	def filter_duplicates = params.pathseq_filter_duplicates == true
	// if (params.pathseq_filter_duplicates) {
	// 	filter_duplicates += "--filter-duplicates true"
	// } else {
	// 	filter_duplicates += "--filter-duplicates false"
	// }

	// if [[ -f "${params.pathseq_db}/pathseq_microbe.fa" ]]; then
	// 	microbe_seq=""
	// else
	// 	microbe_seq="--microbe-dict ${params.pathseq_db}/*.dict"
	// fi
	"""
    mkdir -p ${sample.id}

	gatk --java-options \"-Xmx${maxmem}g\" PathSeqPipelineSpark \\
		--input ${bam} \\
		--filter-bwa-image ${pathseq_db}/pathseq_host.fa.img \\
		--kmer-file ${pathseq_db}/pathseq_host.bfi \\
		--min-clipped-read-length ${params.pathseq_min_clipped_read_length} \\
		--microbe-fasta ${pathseq_db}/pathseq_microbe.fa \\
		--microbe-bwa-image ${pathseq_db}/pathseq_microbe.fa.img \\
		--taxonomy-file ${pathseq_db}/pathseq_taxonomy.db \\
		--output ${sample.id}/${sample.id}.pathseq.bam \\
		--scores-output ${sample.id}/${sample.id}.pathseq.scores \\
		--score-metrics ${sample.id}/${sample.id}.pathseq.score_metrics \\
		--is-host-aligned ${params.pathseq_skip_host_alignment} \\
		--skip-quality-filters ${params.pathseq_skip_quality_filters} \\
		--filter-metrics ${sample.id}/${sample.id}.pathseq.filter_metrics \\
		--filter-duplicates ${filter_duplicates}
    """
}
