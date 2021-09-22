#!/usr/bin/env nextflow

nextflow.enable.dsl=2

merge_paired_end = (!params.skip_merge_pairs) ? true : false;


process qc_preprocess {
	input:
	tuple val(sample), path(reads)

	output:
	tuple val("${sample}"), path("${sample}/${sample}.qc_[RO]*.fastq.gz"), optional: true, emit: qc_reads_p
	tuple val("${sample}"), path("${sample}/${sample}.qc_U.fastq.gz"), optional: true, emit: qc_reads_s

	script:
	//def qc_params = "qtrim=rl trimq=25 maq=25 minlen=45 ref=${params.adapters} ktrim=r k=23 mink=11 hdist=1 tpe tbo"


	if (reads.size() == 2) {
		"""
		maxmem=\$(echo \"$task.memory\"| sed 's/ GB/g/g')
		mkdir -p ${sample}

		bbduk.sh -Xmx\$maxmem t=$task.cpus ${params.qc_params} in=${sample}_R1.fastq.gz in2=${sample}_R2.fastq.gz out=${sample}/${sample}.qc_R1.fastq.gz out2=${sample}/${sample}.qc_R2.fastq.gz outs=${sample}/${sample}.qc_O.fastq.gz
		"""
	} else {
		"""
		maxmem=\$(echo \"$task.memory\"| sed 's/ GB/g/g')
		mkdir -p ${sample}

		bbduk.sh -Xmx\$maxmem t=$task.cpus ${params.qc_params} in=${sample}_R1.fastq.gz out=${sample}/${sample}.qc_U.fastq.gz
		"""
	}
}


process qc_bbmerge {
	input:
	tuple val(sample), path(reads)

	output:
	tuple val(sample), path("${sample}/${sample}*.fastq.gz"), emit: merged_reads

	script:
	def merge_params = "rsem=t extend2=20 iterations=5 ecct vstrict"
	"""
	maxmem=\$(echo \"$task.memory\"| sed 's/ GB/g/g')
	mkdir -p ${sample}

	bbmerge.sh -Xmx\$maxmem t=$task.cpus ${merge_params} in=${sample}.qc_R1.fastq.gz in2=${sample}.qc_R2.fastq.gz out=${sample}/${sample}.merged_M.fastq.gz outu1=${sample}/${sample}.merged_R1.fastq.gz outu2=${sample}/${sample}.merged_R2.fastq.gz
	"""
}


process concat_singles {
	input:
	tuple val(sample), path(reads)

	output:
	tuple val(sample), path("${sample}/${sample}_concat_singles.fastq.gz"), emit: concat_reads

	script:
	"""
	mkdir -p $sample
	cat ${reads} > ${sample}/${sample}_concat_singles.fastq.gz
	"""
}


process fastqc {
    publishDir "${params.output_dir}", mode: params.publish_mode, pattern: "raw_counts/*.txt"

	input:
	tuple val(sample), path(reads)

	output:
	tuple val(sample), path("fastqc/*/*fastqc_data.txt"), emit: reports
	tuple val(sample), path("raw_counts/${sample}.txt"), emit: counts

	script:
	def process_r2 = (reads.size() == 2) ? "fastqc -t $task.cpus --extract --outdir=fastqc ${sample}_R2.fastq.gz && mv fastqc/${sample}_R2_fastqc/fastqc_data.txt fastqc/${sample}_R2_fastqc/${sample}_R2_fastqc_data.txt" : "";

	"""
	mkdir -p fastqc
	mkdir -p raw_counts
	fastqc -t $task.cpus --extract --outdir=fastqc ${sample}_R1.fastq.gz && mv fastqc/${sample}_R1_fastqc/fastqc_data.txt fastqc/${sample}_R1_fastqc/${sample}_R1_fastqc_data.txt
	${process_r2}	

	grep "Total Sequences" fastqc/*/*data.txt > seqcount.txt
	echo \$(wc -l seqcount.txt)\$'\t'\$(head -n1 seqcount.txt | cut -f 2) > raw_counts/${sample}.txt
	"""
}


process multiqc {
    publishDir "${params.output_dir}", mode: params.publish_mode
	
	input:
	path(reports)

	output:
	path("multiqc_report.html")

	script:
	def send_report = (params.email) ? "echo . | mailx -s 'multiqc_report' -a multiqc_report.html ${params.email}" : ""
	"""
	multiqc -c ${params.multiqc_config} .
	${send_report}
	"""
}


workflow nevermore_preprocess {
	take:
		reads_ch

	main: 
		
		/*
			Normalise fastq file names by symlinking to R1/R2.
		*/

		//prepare_fastqs(reads_ch)
		fastqc(reads_ch)

		multiqc(
			fastqc.out.reports.map { sample, report -> report }.collect()
		)


		/*
			Preprocess (QC) reads, then split preprocessed pairs from single-end / orphaned reads.
		*/

		qc_preprocess(reads_ch) //prepare_fastqs.out.reads)

		/*
			Merge paired-end reads, then split merged reads from those that failed to merge.
		*/

		if (merge_paired_end) {

			qc_bbmerge(qc_preprocess.out.qc_reads_p)

			qc_bbmerge.out.merged_reads
				.multiMap { sample, reads ->
					merged: (reads[0] != null && reads[0].name.endsWith(".fastq.gz")) ? [sample, reads[0]] : 0
					paired: (reads[1] != null && reads[1].name.endsWith(".fastq.gz") && reads[2] != null && reads[2].name.endsWith(".fastq.gz")) ? [sample, [reads[1], reads[2]]] : 0
					other: 0
				}
				.set { merged_reads_ch }

			merged_merged_ch = merged_reads_ch.merged.filter({ it != 0 })
			paired_ch = merged_reads_ch.paired.filter({ it != 0 })

			/*
				Redirect all unpaired reads into a common channel, then concatenate them into a single unpaired fastq file.
			*/
	
			single_reads_ch = merged_merged_ch.concat(qc_preprocess.out.qc_reads_s)
				.map { sample, reads ->
					return tuple(sample.replaceAll(/.singles$/, ""), reads)
				}
	        	.groupTuple(sort: true)
	
			concat_singles(single_reads_ch)

			singles_ch = concat_singles.out.concat_reads

		} else {

			merged_merged_ch = Channel.empty()
			paired_ch = qc_preprocess.out.qc_reads_p 
			singles_ch = qc_preprocess.out.qc_reads_s
		
		}

		singles_ch = singles_ch
			.map { sample, reads -> return tuple("${sample}.singles", reads) }

	emit:
		preprocessed = paired_ch.concat(singles_ch)
}
