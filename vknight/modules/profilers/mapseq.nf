process mapseq {
	container "ghcr.io/cschu/vknight_profilers:main"
	label "mapseq"
	publishDir params.output_dir, mode: "copy"

    input:
    tuple val(sample), path(seqs)

    output:
    path("${sample.id}/${sample.id}_R*bac_ssu.mseq"), emit: bac_ssu

    script:
	// def r2_cmd = (sample.is_paired) ? "${params.mapseq_bin} -nthreads ${task.cpus} --outfmt simple ${sample.id}_R2.fastq.gz_bac_ssu.fasta > ${sample.id}/${sample.id}_R2_bac_ssu.mseq" : ""
	def r2_cmd = (sample.is_paired) ? "${params.mapseq_bin} -nthreads ${task.cpus} --outfmt simple ${sample.id}_R2*.fasta > ${sample.id}/${sample.id}_R2_bac_ssu.mseq" : ""

	"""
    mkdir -p ${sample.id}
    ${params.mapseq_bin} -nthreads ${task.cpus} --outfmt simple ${sample.id}_R1*.fasta > ${sample.id}/${sample.id}_R1_bac_ssu.mseq
	${r2_cmd}
	"""
    // ${params.mapseq_bin} -nthreads ${task.cpus} --outfmt simple ${sample.id}_R1.fastq.gz_bac_ssu.fasta > ${sample.id}/${sample.id}_R1_bac_ssu.mseq
}


process mapseq_with_customdb {
	container "ghcr.io/cschu/vknight_profilers:main"
	label "mapseq"
	publishDir params.output_dir, mode: "copy"

	input:
	tuple val(sample), path(seqs)
	path(db)

	output:
	path("${sample.id}/${sample.id}_R*bac_ssu.mseq"), emit: bac_ssu

	script:
	def db_string = "${db}/*.fna ${db}/*.tax"
	// def r2_cmd = (sample.is_paired) ? "${params.mapseq_bin} -nthreads ${task.cpus} --outfmt simple ${sample.id}_R2.fastq.gz_bac_ssu.fasta ${db_string} > ${sample.id}/${sample.id}_R2_bac_ssu.mseq" : ""
	def r2_cmd = (sample.is_paired) ? "${params.mapseq_bin} -nthreads ${task.cpus} --outfmt simple ${sample.id}_R2*.fasta ${db_string} > ${sample.id}/${sample.id}_R2_bac_ssu.mseq" : ""

	"""
    mkdir -p ${sample.id}                                                                                               	
    ${params.mapseq_bin} -nthreads ${task.cpus} --outfmt simple ${sample.id}_R1*.fasta ${db_string} > ${sample.id}/${sample.id}_R1_bac_ssu.mseq
	${r2_cmd}
	"""
    //${params.mapseq_bin} -nthreads ${task.cpus} --outfmt simple ${sample.id}_R1.fastq.gz_bac_ssu.fasta ${db_string} > ${sample.id}/${sample.id}_R1_bac_ssu.mseq
}


process collate_mapseq_tables {
	container "ghcr.io/cschu/vknight_profilers:main"
    publishDir params.output_dir, mode: "copy"

    input:
    path(mapped_seqs)

    output:
    path("mapseq_tables/mapseq_counts_genus_*_bac_ssu.tsv"), emit: ssu_tables

    script:
	"""
	mkdir -p mapseq_tables
	${params.mapseq_bin} -otutable -tl 5 \$(ls *_R1_bac_ssu.mseq) | sed 's/_R1_bac_ssu.mseq//g' > mapseq_tables/mapseq_counts_genus_fwd_bac_ssu.tsv	
	r2=\$(ls *_R2_bac_ssu.mseq | wc -l)
    if [[ ! "\$r2" -eq 0 ]]; then
        ${params.mapseq_bin} -otutable -tl 5 \$(ls *_R2_bac_ssu.mseq) | sed 's/_R2_bac_ssu.mseq//g' > mapseq_tables/mapseq_counts_genus_rev_bac_ssu.tsv
	fi
	"""
}
