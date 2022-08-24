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

process transfer_bams {
	input:
		path(bamfiles)
	output:
		path("bam/*.bam"), emit: bamfiles
	script:
		"""
		mkdir -p bam/
		find . -maxdepth 1 -type l -name '*.bam' | xargs -I {} readlink {} | xargs -I {} rsync -avP {} bam/
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

		fastq_ch = prepare_fastqs.out.paired
			.concat(prepare_fastqs.out.single)
			.map { classify_sample(it[0], it[1]) } 


	emit:
		fastqs = fastq_ch
}

workflow bam_input {
	take:
		bam_ch
	main:
		transfer_bams(bam_ch.collect())
		transfer_bams.out.bamfiles.view()
		bam_ch = transfer_bams.out.bamfiles.flatten()
			.map { file ->
				def sample = file.name.replaceAll(/.bam$/, "")
				return tuple(sample, file)
			}
			.groupTuple(sort: true)
			.map { classify_sample(it[0], it[1]) }

		bam2fq(bam_ch)
		bam_ch = bam2fq.out.reads
			.map { classify_sample(it[0].id, it[1]) }
	emit:
		bamfiles = bam_ch
}