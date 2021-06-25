#!/usr/bin/env nextflow

nextflow.enable.dsl=2

if (!params.file_pattern) {
	params.file_pattern = "**[12].fastq.gz"
}

if (!params.single_file_pattern) {
	params.single_file_pattern = "**singles.fastq.gz"
}

if (!params.file_suffix) {
	params.file_suffix = ".fastq.gz"
}

if (!params.preprocessed_singles) {
	run_singles = false
} else {
	run_singles = true
}

if (!params.publish_mode) {
	params.publish_mode = "link"
}

if (!params.output_dir) {
	params.output_dir = "nevermore_out"
}

output_dir = "${params.output_dir}"

suffix_pattern = params.file_suffix


process qc_preprocess {
	conda "bioconda::bbmap"

	input:
	tuple val(sample), path(reads)

	output:
	stdout
	val "${sample}", emit: sample_id
	path "${sample}/${sample}.qc_R1.fastq.gz", emit: r1_fq
    path "${sample}/${sample}.qc_R2.fastq.gz", optional: true, emit: r2_fq
	path "${sample}/${sample}.qc_S1.fastq.gz", optional: true, emit: u_fq

	script:
	def qc_params = "qtrim=rl trimq=25 maq=25 minlen=45"
	def r1_index = reads[0].name.endsWith("1${suffix_pattern}") ? 0 : 1
	def r1 = "in=" + reads[r1_index] + " out=${sample}/${sample}.qc_R1.fastq.gz"
	def r2 = (reads[1] && file(reads[1])) ? "in2=" + ( reads[ r1_index == 0 ? 1 : 0 ] + " out2=${sample}/${sample}.qc_R2.fastq.gz outs=${sample}/${sample}.qc_S1.fastq.gz" ) : ""

	"""
	maxmem=\$(echo \"$task.memory\"| sed 's/ GB/g/g')
	bbduk.sh -Xmx\$maxmem t=$task.cpus ${qc_params} ${r1} ${r2}
	"""
}

process qc_preprocess_singles {
	conda "bioconda::bbmap"

	input:
	tuple val(sample), path(reads)

	output:
	stdout
	val "$sample", emit: sample_id
	path "${sample}/${sample}.qc_U.fastq.gz", emit: u_fq

	script:
	def qc_params = "qtrim=rl trimq=25 maq=25 minlen=45"

	"""
	maxmem=\$(echo \"$task.memory\"| sed 's/ GB/g/g')
	bbduk.sh -Xmx\$maxmem t=$task.cpus ${qc_params} in=${reads[0]} out=${sample}/${sample}.qc_U.fastq.gz
	"""
}

process merge_singles {

	input:
	val(sample)
	path(reads)
	path(singles)

	output:
	stdout
	val "$sample", emit: sample_id
    path "${sample}/${sample}.qc_U.fastq.gz", emit: u_fq

	script:
	"""
	mkdir -p $sample
	cat ${reads} ${singles} > ${sample}/${sample}.qc_U.fastq.gz
	"""

}

process decontaminate {
	conda "bioconda::bwa bioconda::'samtools>=1.11'"

	input:
	val(sample)
	path(reads1)
	path(reads2)

	output:
	stdout
	val "$sample", emit: sample_id
	path "${sample}/${sample}.no_human_R1.fastq.gz", emit: r1_fq
	path "${sample}/${sample}.no_human_R2.fastq.gz", optional: true, emit: r2_fq

	script:
	def r1 = reads1
	def r2 = file(reads2) ? reads2 : ""

	def r1_out = "-1 ${sample}/${sample}.no_human_R1.fastq.gz"
	def r2_out = file(reads2) ? "-2 ${sample}/${sample}.no_human_R2.fastq.gz" : ""

	"""
	cpus=\$(expr \"$task.cpus\" - 4)
	mkdir -p $sample
	bwa mem -t \$cpus ${params.human_ref} ${r1} ${r2} | samtools collate --threads 2 -f -O - | samtools fastq -@ 2 -f 4 ${r1_out} ${r2_out} -s singletons.fastq.gz -
	"""
}

process decontaminate_singles {
	conda "bioconda::bwa bioconda::'samtools>=1.11'"

	input:
	val(sample)
	path(reads)

	output:
	stdout
	val "$sample", emit: sample_id
	path "${sample}/${sample}.no_human_U.fastq.gz", emit: u_fq

	script:

	"""
	cpus=\$(expr \"$task.cpus\" - 4)
	mkdir -p $sample
	bwa mem -t \$cpus ${params.human_ref} ${reads} | samtools collate --threads 2 -f -O - | samtools fastq -@ 2 -f 4 -s ${sample}/${sample}.no_human_U.fastq.gz -
	"""

}

process align {
	conda "bioconda::bwa bioconda::'samtools>=1.11'"

	input:
	val(sample)
	path(reads1)
	path(reads2)

	output:
	stdout
	val "$sample", emit: sample_id
	path "${sample}.main.bam", emit: bam

	script:
	def r1 = reads1
	def r2 =  file(reads2) ? reads2 : ""

	"""
	cpus=\$(expr \"$task.cpus\" - 4)
	bwa mem -a -t \$cpus ${params.reference} ${r1} ${r2} | samtools view -F 4 -buSh - | samtools sort -@ 2 -o ${sample}.main.bam -
	"""
}

process align_singles {
	conda "bioconda::bwa bioconda::'samtools>=1.11'"

	input:
	val(sample)
	path(reads)

	output:
	stdout
	val "$sample", emit: sample_id
	path "${sample}.main.singles.bam", emit: bam

	script:

	"""
	cpus=\$(expr \"$task.cpus\" - 4)
	bwa mem -a -t \$cpus ${params.reference} ${reads} | samtools view -F 4 -buSh - | samtools sort -@ 2 -o ${sample}.main.singles.bam -
	"""
}

process rename_bam {
	publishDir "$output_dir", mode: params.publish_mode

	input:
	val(sample)
	path(bam)

	output:
	stdout
	val "$sample", emit: sample_id
	path "${sample}/${sample}.bam", emit: bam

	script:

	"""
	mkdir -p $sample
	cp $bam ${sample}/${sample}.bam
	"""
}

process merge_and_sort {
	conda "bioconda::bwa bioconda::'samtools>=1.11'"
	publishDir "$output_dir", mode: params.publish_mode

	input:
	val(sample)
	path(main_bam)
	path(singles_bam)

	output:
	stdout
	val "$sample", emit: sample_id
	path "${sample}/${sample}.bam", emit: bam

	script:
	def s_bam = file(singles_bam) ? singles_bam : ""

	if (file(singles_bam)) {
		"""
		mkdir -p $sample
		samtools merge -@ $task.cpus "${sample}/${sample}.bam" ${main_bam} ${s_bam}
		"""
	} else {
		"""
		mkdir -p $sample
		cp ${main_bam} "${sample}/${sample}.bam"
		"""
	}
}

process gffquant {
	conda params.gffquant_env
	publishDir "$output_dir", mode: params.publish_mode

	input:
	val(sample)
	path(bam)

	output:
	stdout
	val "$sample", emit: sample_id
	path "${sample}/${sample}.seqname.uniq.txt", emit: uniq_seq
	path "${sample}/${sample}.seqname.dist1.txt", emit: dist1_seq
	path "${sample}/${sample}.feature_counts.txt", emit: feat_counts
	path "${sample}/${sample}.gene_counts.txt", emit: gene_counts
	path "${sample}/${sample}.covsum.txt", emit: covsum

	script:

	"""
	mkdir -p ${sample}
	gffquant ${params.gffquant_db} ${bam} -o ${sample}/${sample} -m ${params.gffquant_mode} --ambig_mode ${params.gffquant_amode}
	"""
}



workflow {
	reads_ch = Channel
	    .fromPath(params.input_dir + "/" + params.file_pattern)
    	.map { file ->
        	def sample = file.name.replaceAll(suffix_pattern, "")
        	sample = sample.replaceAll(/_[12]$/, "")
	        return tuple(sample, file)
    	}
	    .groupTuple()
	reads_ch.view()

	aux_reads_ch = Channel
		.fromPath(params.input_dir + "/" + params.single_file_pattern)
		.map { file ->
			def sample = file.name.replaceAll(suffix_pattern, "")
			sample = sample.replaceAll(/.singles/, "")
			return tuple(sample, file)
		}
		.groupTuple()
	aux_reads_ch.view()

	if (params.r1_only) {
		qc_preprocess_singles(reads_ch)
		decontaminate_singles(qc_preprocess_singles.out.sample_id, qc_preprocess_singles.out.u_fq)
		align_singles(decontaminate_singles.out.sample_id, decontaminate_singles.out.u_fq)
		rename_bam(align_singles.out.sample_id, align_singles.out.bam)
		gffquant(rename_bam.out.sample_id, rename_bam.out.bam)
	} else {
		qc_preprocess(reads_ch)
		decontaminate(qc_preprocess.out.sample_id, qc_preprocess.out.r1_fq, qc_preprocess.out.r2_fq)
		align(decontaminate.out.sample_id, decontaminate.out.r1_fq, decontaminate.out.r2_fq)
		if (run_singles) {
			qc_preprocess_singles(aux_reads_ch)
			merge_singles(qc_preprocess.out.sample_id, qc_preprocess.out.u_fq, qc_preprocess_singles.out.u_fq)
			decontaminate_singles(merge_singles.out.sample_id, merge_singles.out.u_fq)
			align_singles(decontaminate_singles.out.sample_id, decontaminate_singles.out.u_fq)
			merge_and_sort(align.out.sample_id, align.out.bam, align_singles.out.bam)
		} else {
			decontaminate_singles(qc_preprocess.out.sample_id, qc_preprocess.out.u_fq)
			align_singles(decontaminate_singles.out.sample_id, decontaminate_singles.out.u_fq)
			merge_and_sort(align.out.sample_id, align.out.bam, align_singles.out.bam)
		}
		gffquant(merge_and_sort.out.sample_id, merge_and_sort.out.bam)
	}

}
