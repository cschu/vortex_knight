params.motus_tax_level = "mOTU"  // -k
params.motus_min_length = 75  // -l
params.motus_n_marker_genes = 3  // -g 
params.motus_run_mapsnv = false
params.motus_full_rank_taxonomy = false  // -q 
params.motus_print_counts = false // -c 
params.motus_readcount_type = "insert.scaled_counts"  // -y

// deal with motus readcount types
// motus < 4: base.coverage, insert.raw_counts, insert.scaled_counts
// motus 4: INSERT_RAW, INSERT_NORM, INSERT_SCALED, BASE_RAW, BASE_NORM
def motus_readcount_type = "insert.scaled_counts"
def motus4_readcount_type = "INSERT_SCALED"
if (params.motus_readcount_type in [ "base.coverage", "BASE_RAW" ]) {
    motus_readcount_type = "base.coverage"
    motus4_readcount_type = "BASE_RAW"
} else if (params.motus_readcount_type in [ "insert.raw_counts", "INSERT_RAW" ]) {
    motus_readcount_type = "insert.raw_counts"
    motus4_readcount_type = "INSERT_RAW"
} else if (params.motus_readcount_type in [ "insert.scaled_counts", "INSERT_SCALED" ]) {
    motus_readcount_type = "insert.scaled_counts"
    motus4_readcount_type = "INSERT_SCALED"
} else if (params.motus_readcount_type in [ "INSERT_NORM", "BASE_NORM" ]) {
    motus_readcount_type = null
    motus4_readcount_type = param.motus_readcount_type
}







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
    motus profile -v 7 ${additional_options.join(" ")} \
	-t ${task.cpus} \
	-n ${sample.id} \
	-k ${params.motus_tax_level} \
	-l ${params.motus_min_length} \
	-g ${params.motus_n_marker_genes} \
	-y ${motus_readcount_type} \
	-db ${motus_db} \
    -o ${sample.id}/${sample.id}.motus.txt \
	${input_files}

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

    """
    set -e -o pipefail
    mkdir -p ${sample.id}

    export MOTUS_DB=\$(readlink ${motus_db})

    motus profile -v 7 \
	-t ${task.cpus} \
	-n ${sample.id} \
	-l ${params.motus_min_length} \
	-g ${params.motus_n_marker_genes} \
	-y ${motus4_readcount_type} \
	-o ${sample.id}/${sample.id}.motus4.txt \
	${input_files}

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