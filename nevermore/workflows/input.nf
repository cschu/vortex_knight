nextflow.enable.dsl=2

include { classify_sample } from "../modules/functions"


process transfer_fastqs {
	input:
		path(fastqs)
	output:
		path("fastq/*/**"), emit: fastqs
	script:
		"""
		transfer.py -i . -o fastq/
		"""
}

process prepare_fastqs {
	publishDir params.output_dir, mode: "${params.publish_mode}"

	input:
		tuple val(sample_id), path(files)
	output:
		tuple val(sample_id), path("${sample_id}/*.fastq.gz"), emit: paired, optional: true
		tuple val("${sample_id}.singles"), path("${sample_id}.singles/*.fastq.gz"), emit: single, optional: true
    script:
		"""
		prepare_fastqs.py -i . -o . -s ${sample_id} --remove-suffix ${params.suffix_pattern} > run.sh
     	bash run.sh
        """
}

workflow remote_fastq_input {
	take:
		fastq_ch

	main:

		transfer_fastqs(fastq_ch.collect())
		res_ch = transfer_fastqs.out.fastqs.flatten()

	emit:
		fastqs = res_ch
}


workflow fastq_input {
	take:
		fastq_ch
	
	main:
		if (params.remote_input_dir) {
			fastq_ch = remote_fastq_input(fastq_ch)
		}

		fastq_ch = fastq_ch
			.map { file ->
    			def sample = file.getParent().getName()
				return tuple(sample, file)
			}
			.groupTuple(sort: true)

		prepare_fastqs(fastq_ch)


	emit:
		fastqs = prepare_fastqs.out.paired
			.concat(prepare_fastqs.out.single)
			.map { classify_sample(it[0], it[1]) }
}


/*
		


    179     fastq_ch = Channel
    180         .fromPath(params.input_dir + "/**[._]{fastq.gz,fq.gz}").collect()
    181     // fastq_ch.view()
    182
    183     copy_or_gzip(fastq_ch, "/scratch/schudoma/copy_or_zip.py")
    184     copy_or_gzip.out.fastqs.view()
    185
    186     fastq_ch = copy_or_gzip.out.fastqs.flatten()
    187
    188
    189         .map { file ->
    190             def sample = file.getParent().getName()
    191             fname = file.name.replaceAll(/.(fastq|fq)(.gz)?$/, "").replaceAll(params.suffix_pattern, "")
    192             return tuple(sample, file)
    193         }
    194         .groupTuple(sort: true)
    195     // fastq_ch.view()
    196
    197     rename_or_merge_fastqs(fastq_ch, "/scratch/schudoma/get_samples.py")
    198
    199
    200 }
*/	