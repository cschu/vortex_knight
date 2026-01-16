params.motus_tax_level = "mOTU"
params.motus_min_length = 75
params.motus_n_marker_genes = 3
params.motus_readcount_type = "insert.scaled_counts"
params.motus_run_mapsnv = false
params.motus_full_rank_taxonomy = false
params.motus_print_counts = false

process motus {
    container "quay.io/biocontainers/motus:3.1.0--pyhdfd78af_0"

    input:
    tuple val(sample), path(fastqs)
	path(motus_db)

    output:
    tuple val(sample), path("${sample.id}/${sample.id}.motus.txt"), emit: motus_profile
    tuple val(sample), path("${sample.id}/${sample.id}.motus.bam"), emit: motus_bam, optional: true
    tuple val(sample), path("MOTUS_DONE_SENTINEL")

    script:

    def input_files = ""
	def r1_files = fastqs.findAll( { it.name.endsWith("_R1.fastq.gz") && !it.name.matches("(.*)(singles|orphans|chimeras)(.*)") } )
	def r2_files = fastqs.findAll( { it.name.endsWith("_R2.fastq.gz") } )
	def orphans = fastqs.findAll( { it.name.matches("(.*)(singles|orphans|chimeras)(.*)") } )

	if (r1_files.size() != 0 && r2_files.size() != 0) {
		input_files += "-f ${r1_files.join(' ')} -r ${r2_files.join(' ')}"
	} else if (r1_files.size() != 0) {
		input_files += "-s ${r1_files.join(' ')}"
	} else if (r2_files.size() != 0) {
		input_files += "-s ${r2_files.join(' ')}"
	}
    
    if (orphans.size() != 0) {
		input_files += " -s ${orphans.join(' ')}"
	}

    def mapsnv_cmd = ""
    if (params.motus_run_mapsnv) {
        mapsnv_cmd += "motus map_snv -t ${task.cpus} -db ${motus_db} ${input_files} > ${sample.id}/${sample.id}.motus.bam"
    } // else {
    //     mapsnv_cmd += "touch ${sample.id}/${sample.id}.motus.bam"
    // }

    def additional_options = []
    if (params.motus_full_rank_taxonomy == true) {
        additional_options.add("-q")
    }
    
    if (params.motus_print_counts) {
        additional_options.add("-c")
    }

    """
    set -e -o pipefail
    mkdir -p ${sample.id}
    motus profile -v 7 -db ${motus_db} ${additional_options.join(" ")} \
	-t ${task.cpus} \
	-n ${sample.id} \
	-k ${params.motus_tax_level} \
	-l ${params.motus_min_length} \
	-g ${params.motus_n_marker_genes} \
	-y ${params.motus_readcount_type} \
	-db ${motus_db} \
	${input_files} > ${sample.id}/${sample.id}.motus.txt

    ${mapsnv_cmd}
    touch MOTUS_DONE_SENTINEL
    """

    // # generate profile
    // motus profile \
    // -f QC/${ID}.R1.fastq.gz -r QC/${ID}.R2.fastq.gz \
    // -g 1 \
    // -t 30 \
    // -y insert.raw_counts \
    // -o out_align/${ID}.profile.txt 
}


process motus4 {
    container "registry.git.embl.org/schudoma/motus4-docker:latest"

    input:
    tuple val(sample), path(fastqs)
	path(motus_db)

    output:
    tuple val(sample), path("${sample.id}/${sample.id}.motus4.txt"), emit: motus_profile
    tuple val(sample), path("MOTUS4_DONE_SENTINEL")

    script:

    def input_files = ""
	def r1_files = fastqs.findAll( { it.name.endsWith("_R1.fastq.gz") && !it.name.matches("(.*)(singles|orphans|chimeras)(.*)") } )
	def r2_files = fastqs.findAll( { it.name.endsWith("_R2.fastq.gz") } )
	def orphans = fastqs.findAll( { it.name.matches("(.*)(singles|orphans|chimeras)(.*)") } )

	if (r1_files.size() != 0 && r2_files.size() != 0) {
		input_files += "-f ${r1_files.join(' ')} -r ${r2_files.join(' ')}"
	} else if (r1_files.size() != 0) {
		input_files += "-s ${r1_files.join(' ')}"
	} else if (r2_files.size() != 0) {
		input_files += "-s ${r2_files.join(' ')}"
	}
    
    if (orphans.size() != 0) {
		input_files += " -s ${orphans.join(' ')}"
	}

    def additional_options = []
    if (params.motus_full_rank_taxonomy == true) {
        additional_options.add("-q")
    }
    
    if (params.motus_print_counts) {
        additional_options.add("-c")
    }

    """
    set -e -o pipefail
    mkdir -p ${sample.id}
    motus profile -v 7 -db ${motus_db} ${additional_options.join(" ")} \
	-t ${task.cpus} \
	-n ${sample.id} \
	-k ${params.motus_tax_level} \
	-l ${params.motus_min_length} \
	-g ${params.motus_n_marker_genes} \
	-y ${params.motus_readcount_type} \
	-db ${motus_db} \
	${input_files} > ${sample.id}/${sample.id}.motus4.txt

    touch MOTUS4_DONE_SENTINEL
    """
}



process motus_merge {
    container "quay.io/biocontainers/motus:3.1.0--pyhdfd78af_0"

    input:
    path(profiles)
    path(motus_db)

    output:
    path("motus_profiles/motus_merged.txt")

    script:
    """
    mkdir -p motus_profiles/ input/

    for f in ${profiles}; do ln -sf ../\$f input/; done

    motus merge -db ${motus_db} -d input/ -o motus_profiles/motus_merged.txt
    """

}