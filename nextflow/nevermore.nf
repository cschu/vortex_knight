#!/usr/bin/env nextflow

nextflow.enable.dsl=2

if (!params.file_pattern) {
	params.file_pattern = "**[12].fastq.gz"
}

if (!params.single_file_pattern) {
	params.single_file_pattern = "**singles.fastq.gz"
}

if (!params.publish_mode) {
	params.publish_mode = "link"
}

if (!params.output_dir) {
	params.output_dir = "nevermore_out"
}

output_dir = "${params.output_dir}"

suffix_pattern = params.file_pattern.replaceAll(/\*\*/, "").substring(4)
single_suffix_pattern = params.single_file_pattern.replaceAll(/\*\*/, "").substring(7)

//qc_params = "qtrim=rl trimq=25 maq=25 minlen=45"

process qc_preprocess {
	conda "bioconda::bbmap"
    publishDir "$output_dir", mode: params.publish_mode

	input:
	tuple val(sample), path(reads)

	output:
	stdout
	val "${sample}", emit: sample_id
	path "${sample}/${sample}.qc_R1.fastq.gz", emit: r1_fq
    path "${sample}/${sample}.qc_R2.fastq.gz", optional: true, emit: r2_fq
	path "${sample}/${sample}.qc_S1.fastq.gz", emit: singles_fq

	script:
	def qc_params = "qtrim=rl trimq=25 maq=25 minlen=45"
	def r1_index = reads[0].name.endsWith("1${suffix_pattern}") ? 0 : 1
	def r1 = "in=" + reads[r1_index] + " out=${sample}/${sample}.qc_R1.fastq.gz"
	def r2 = file(reads[1]) ? "in2=" + ( reads[ r1_index == 0 ? 1 : 0 ] + " out2=${sample}/${sample}.qc_R2.fastq.gz outs=${sample}/${sample}.qc_S1.fastq.gz" ) : ""

	"""
	maxmem=\$(echo \"$task.memory\"| sed 's/ GB/g/g')
	bbduk.sh -Xmx\$maxmem t=$task.cpus ${qc_params} ${r1} ${r2}
	touch ${sample}/${sample}.qc_S1.fastq.gz
	"""
}

process qc_preprocess_singles {
	conda "bioconda::bbmap"
	publishDir "$output_dir", mode: params.publish_mode

	input:
	tuple val(sample), path(reads)
	path(singles)

	output:
	stdout
	val "$sample", emit: sample_id
	path "${sample}/${sample}.qc_U.fastq.gz", emit: u_fq

	script:
	def qc_params = "qtrim=rl trimq=25 maq=25 minlen=45"

	"""
	maxmem=\$(echo \"$task.memory\"| sed 's/ GB/g/g')
	bbduk.sh -Xmx\$maxmem t=$task.cpus ${qc_params} in=${reads[0]} out=${sample}/${sample}.qc_S.fastq.gz
	cat ${singles} ${sample}/${sample}.qc_S.fastq.gz > ${sample}/${sample}.qc_U.fastq.gz
	"""
}

process decontaminate {
	conda "bioconda::bwa bioconda::'samtools>=1.11'"
	publishDir "$output_dir", mode: params.publish_mode

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
	publishDir "$output_dir", mode: params.publish_mode

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
	publishDir "$output_dir", mode: params.publish_mode

	input:
	val(sample)
	path(reads1)
	path(reads2)

	output:
	stdout
	val "$sample", emit: sample_id
	path "${sample}/${sample}.main.bam", emit: bam

	script:
	def r1 = reads1
	def r2 =  file(reads2) ? reads2 : ""

	"""
	cpus=\$(expr \"$task.cpus\" - 4)
	mkdir -p $sample
	bwa mem -a -t \$cpus ${params.reference} ${r1} ${r2} | samtools view -F 4 -buSh - | samtools sort -@ 2 -o ${sample}/${sample}.main.bam -
	"""
}

process align_singles {
	conda "bioconda::bwa bioconda::'samtools>=1.11'"
	publishDir "$output_dir", mode: params.publish_mode

	input:
	val(sample)
	path(reads)

	output:
	stdout
	val "$sample", emit: sample_id
	path "${sample}/${sample}.main.singles.bam", emit: bam

	script:

	"""
	cpus=\$(expr \"$task.cpus\" - 4)
	mkdir -p $sample
	bwa mem -a -t \$cpus ${params.reference} ${reads} | samtools view -F 4 -buSh - | samtools sort -@ 2 -o ${sample}/${sample}.main.singles.bam -
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

	"""
	samtools merge -@ $task.cpus "${sample}/${sample}.bam" ${main_bam} ${singles_bam}
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
			def sample = file.name.replaceAll(single_suffix_pattern, "")
			sample = sample.replaceAll(/.singles/, "")
			return tuple(sample, file)
		}
		.groupTuple()
	aux_reads_ch.view()

	qc_preprocess(reads_ch)
	qc_preprocess_singles(aux_reads_ch, qc_preprocess.out.singles_fq)

	decontaminate(qc_preprocess.out.sample_id, qc_preprocess.out.r1_fq, qc_preprocess.out.r2_fq)
	decontaminate_singles(qc_preprocess_singles.out.sample_id, qc_preprocess_singles.out.u_fq)

	align(decontaminate.out.sample_id, decontaminate.out.r1_fq, decontaminate.out.r2_fq)
	align_singles(decontaminate_singles.out.sample_id, decontaminate_singles.out.u_fq)

	merge_and_sort(align.out.sample_id, align.out.bam, align_singles.out.bam)
}


/*in=<file>           Main input. in=stdin.fq will pipe from stdin.
in2=<file>          Input for 2nd read of pairs in a different file.
ref=<file,file>     Comma-delimited list of reference files.
                    In addition to filenames, you may also use the keywords:
                    adapters, artifacts, phix, lambda, pjet, mtst, kapa
literal=<seq,seq>   Comma-delimited list of literal reference sequences.
touppercase=f       (tuc) Change all bases upper-case.
interleaved=auto    (int) t/f overrides interleaved autodetection.
qin=auto            Input quality offset: 33 (Sanger), 64, or auto.
reads=-1            If positive, quit after processing X reads or pairs.
copyundefined=f     (cu) Process non-AGCT IUPAC reference bases by making all
                    possible unambiguous copies.  Intended for short motifs
                    or adapter barcodes, as time/memory use is exponential.
samplerate=1        Set lower to only process a fraction of input reads.
samref=<file>       Optional reference fasta for processing sam files.

Output parameters:
out=<file>          (outnonmatch) Write reads here that do not contain
                    kmers matching the database.  'out=stdout.fq' will pipe
                    to standard out.
out2=<file>         (outnonmatch2) Use this to write 2nd read of pairs to a
                    different file.
*/
