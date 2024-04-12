process run_metaphlan4 {
	container "docker://quay.io/biocontainers/metaphlan:4.1.0--pyhca03a8a_0"
	
	input:
	tuple val(sample), path(fastqs)
	path(mp4_db)

	output:
	tuple val(sample), path("${sample.id}.mp4.txt"), emit: mp4_table
	tuple val(sample), path("${sample.id}.mp4.sam.bz2"), emit: mp4_sam, optional: (params.run_samestr || params.samestr_compatible_output)
	// tuple val(sample), path("${sample.id}.bowtie2.bz2"), emit: mp4_bt2
	
	script:
	def mp4_params = "--bowtie2db ${mp4_db} --input_type fastq --nproc ${task.cpus} --tmp_dir tmp/"
	def mp4_input = ""
	def bt2_out = "--bowtie2out ${sample.id}.bowtie2.bz2"

	def samestr_params = ""
	if (params.run_samestr || params.samestr_compatible_output) {
		samestr_params = "--samout ${sample.id}.mp4.sam.bz2"
	}

	
	if (fastqs instanceof Collection && fastqs.size() == 2) {
		mp4_input = "${sample.id}_R1.fastq.gz,${sample.id}_R2.fastq.gz"
	} else if (fastqs instanceof Collection && fastqs.size() == 3) {
		mp4_input = "${sample.id}_R1.fastq.gz,${sample.id}_R2.fastq.gz,${sample.id}.singles_R1.fastq.gz"
	} else {
		mp4_input = "${fastqs}"
	}

	"""
	mkdir -p tmp/

	metaphlan ${mp4_input} ${mp4_params} ${bt2_out} -o ${sample.id}.mp4.txt ${samestr_params}
	touch ${sample.id}.mp4.sam.bz2
	"""
}


// process convert_mp4_to_gtdb {

// 	input:
// 	tuple val(sample), path(mp4_table)
// 	path(mp4_db)

// 	output:
// 	tuple val(sample), path("${sample.id}.mp4.gtdb.txt"), emit: mp4_table_gtdb

// 	script:
// 	"""
// 	sgb_to_gtdb_profile.py -i ${mp4_table}
// 	"""

// }


process combine_metaphlan4 {
	container "docker://quay.io/biocontainers/metaphlan:4.1.0--pyhca03a8a_0"

	input:
	tuple val(sample), path(bt2)

	output:
	tuple val(sample), path("metaphlan4/${sample.id}/${sample.id}.mp4.txt"), emit: mp4_table

	script:
	def mp4_params = "--input_type bowtie2out --nproc ${task.cpus} --tmp_dir tmp/"
	def bt2_out = "--bowtie2out ${sample.id}.bowtie2.bz2"
	def mp4_input = "${sample.id}.bowtie2.bz2,${sample.id}.singles.bowtie2.bz2"
	"""
	mkdir -p metaphlan4/${sample.id}/

	metaphlan ${mp4_input} ${mp4_params} -o metaphlan4/${sample.id}/${sample.id}.mp4.txt
	"""
}


process collate_metaphlan4_tables {
	container "docker://quay.io/biocontainers/metaphlan:4.1.0--pyhca03a8a_0"

	input:
	path(tables)

	output:
	path("metaphlan4_abundance_table.txt"), emit: mp4_abundance_table

	script:
	"""
	mkdir -p metaphlan4/

	merge_metaphlan_tables.py ${tables} > metaphlan4_abundance_table.txt
	"""

}