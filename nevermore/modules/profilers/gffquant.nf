process run_gffquant {
	publishDir "${params.output_dir}", mode: params.publish_mode

	input:
	tuple val(sample), path(bam)
	path(db)

	output:
	tuple val(sample), path("profiles/${sample}/*.txt.gz"), emit: results

	script:
	def emapper_version = (params.emapper_version) ? "--emapper_version ${params.emapper_version}" : ""
	"""
	echo $sample $bam
	mkdir -p logs/
	mkdir -p profiles/
	gffquant ${db} ${bam} -o profiles/${sample}/${sample} ${params.gffquant_params} > logs/${sample}.o 2> logs/${sample}.e
	"""
}


process collate_feature_counts {
	publishDir "${params.output_dir}", mode: params.publish_mode

	input:
	tuple val(sample), path(count_tables)

	output:
	path("profiles/collated/*.txt.gz"), emit: collated, optional: true

	script:
	"""
	mkdir -p profiles/collated/
	collate_counts . -o profiles/collated/collated -c uniq_scaled
	collate_counts . -o profiles/collated/collated -c combined_scaled
	"""
}
