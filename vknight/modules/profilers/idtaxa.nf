params.idtaxa_error_rate_threshold = 40
params.idtaxa_strand = "both"

process idtaxa {

	input:
	tuple val(sample), path(fasta)
	path(idtaxa_classifier_db)

	script:

	def r1_cmd = "Idtaxa_function.R ${fasta[0]} ${idtaxa_classifier_db} ${params.idtaxa_error_rate_threshold} ${params.idtaxa_strand}"
	def r2_cmd = (sample.is_paired) ? "Idtaxa_function.R ${fasta[0]} ${idtaxa_classifier_db} ${params.idtaxa_error_rate_threshold} ${params.idtaxa_strand}" : ""

	"""
	mkdir -p ${sample.id}/

	${r1_cmd}
	${r2_cmd}
	"""

}