#!/usr/bin/env nextflow

nextflow.enable.dsl=2


if (!params.publish_mode) {
	params.publish_mode = "symlink"
}

if (!params.output_dir) {
	params.output_dir = "vknight_out"
}

output_dir = "${params.output_dir}"

if (params.motus_database) {
	motus_database = "-db ${params.motus_database}"
} else {
	motus_database = ""
}


process bam2fq {
	input:
	tuple val(sample), path(bam)

	output:
	stdout
	tuple val(sample), path("out/${sample}*.fastq.gz"), emit: reads

	script:
	"""
	mkdir -p out
	samtools collate -@ $task.cpus -u -O $bam | samtools fastq -F 3840 -0 ${sample}_other.fastq.gz -1 ${sample}_R1.fastq.gz -2 ${sample}_R2.fastq.gz

	if [[ -z "\$(gzip -dc ${sample}_R1.fastq.gz | head -n 1)" ]];
	then
		mv ${sample}_other.fastq.gz out/${sample}_R1.fastq.gz;
	else
		mv ${sample}_R1.fastq.gz out/;
		if [[ ! -z "\$(gzip -dc ${sample}_R2.fastq.gz | head -n 1)" ]];
		then
			mv ${sample}_R2.fastq.gz out/;
		fi;
	fi;
	rm -rf *.fastq.gz
	"""
}


process prepare_fastqs {
	input:
	tuple val(sample), path(fq)

	output:
	tuple val(sample), path("out/${sample}_R*.fastq.gz") ,emit: reads

	script:
	if (fq.size() == 2) {
		"""
		mkdir -p out
		ln -sf ../${fq[0]} out/${sample}_R1.fastq.gz
		ln -sf ../${fq[1]} out/${sample}_R2.fastq.gz
		"""
	} else {
		"""
		mkdir -p out
		ln -sf ../${fq[0]} out/${sample}_R1.fastq.gz
		"""
	}
}


process fq2bam {
	input:
	tuple val(sample), path(fq)

	output:
	stdout
	tuple val(sample), path("out/${sample}.bam"), emit: reads

	script:

	if (fq.size() == 2) {
		"""
		mkdir -p out
		gatk FastqToSam -F1 ${fq[0]} -F2 ${fq[1]} -O out/${sample}.bam -SM $sample
		"""
	} else {
		"""
		mkdir -p out
		gatk FastqToSam -F1 ${fq[0]} -O out/${sample}.bam -SM $sample
		"""
	}
}


process count_reads {
	publishDir "$output_dir", mode: params.publish_mode

	input:
	tuple val(sample), path(bam)

	output:
	//tuple val(sample),
	path("${sample}/${sample}.libsize.txt"), emit: counts
	path("${sample}/${sample}.is_paired.txt"), emit: is_paired
	path("${sample}/${sample}.flagstats.txt"), emit: flagstats

	script:
	"""
	mkdir -p ${sample}
	samtools flagstat $bam > "${sample}/${sample}.flagstats.txt"
	head -n 1 "${sample}/${sample}.flagstats.txt" | awk '{print \$1 + \$3}' > "${sample}/${sample}.libsize.txt"
	grep -m 1 "paired in sequencing" "${sample}/${sample}.flagstats.txt" | awk '{npaired = \$1 + \$3; if (npaired==0) {print "unpaired"} else {print "paired"};}' > "${sample}/${sample}.is_paired.txt"
	"""
}


process kraken2 {
	publishDir "$output_dir", mode: params.publish_mode

	input:
	tuple val(sample), path(reads)

	output:
	tuple val(sample), path("${sample}/${sample}.kraken2_report.txt"), emit: kraken2_out

	script:
	def is_paired = (reads.size() == 2) ? "--paired" : "";
	"""
	mkdir -p ${sample}
	kraken2 --db ${params.kraken_database} --threads $task.cpus --gzip-compressed --report ${sample}/${sample}.kraken2_report.txt ${is_paired} \$(ls $reads)
	"""
}


process mtag_extraction {
	input:
	tuple val(sample), path(reads)

	output:
	tuple val(sample), path("*_bac_ssu.fasta"), emit: mtag_out

	script:
	def mtag_input = (reads.size() == 2) ? "-i1 ${sample}_R1.fastq.gz -i2 ${sample}_R2.fastq.gz" : "-is ${sample}_R1.fastq.gz";

	"""
	mtags extract -t $task.cpus -o . ${mtag_input}
	"""
}


process mapseq {
	input:
	tuple val(sample), path(seqs)

	output:
	path("${sample}/${sample}_R*bac_ssu.mseq"), emit: bac_ssu

	script:
	if (seqs.size() == 2) {
		"""
		mkdir -p ${sample}
		${params.mapseq_bin} ${sample}_R1.fastq.gz_bac_ssu.fasta > ${sample}/${sample}_R1_bac_ssu.mseq
		${params.mapseq_bin} ${sample}_R2.fastq.gz_bac_ssu.fasta > ${sample}/${sample}_R2_bac_ssu.mseq
		"""
	} else {
		"""
		mkdir -p ${sample}
		${params.mapseq_bin} ${sample}_R1.fastq.gz_bac_ssu.fasta > ${sample}/${sample}_R1_bac_ssu.mseq
		"""
	}
}


process collate_mapseq_tables {
	publishDir "$output_dir", mode: params.publish_mode

	input:
	path(mapped_seqs)

	output:
	path("otu_tables/mapseq_counts_genus_*_bac_ssu.tsv"), emit: ssu_tables

	script:
	if (mapped_seqs.size() % 2 == 0) {
		"""
		mkdir -p otu_tables
		${params.mapseq_bin} -otutable -tl 5 \$(ls *_R1_bac_ssu.mseq) | sed 's/_R1_bac_ssu.mseq//g' > otu_tables/mapseq_counts_genus_fwd_bac_ssu.tsv
		${params.mapseq_bin} -otutable -tl 5 \$(ls *_R2_bac_ssu.mseq) | sed 's/_R2_bac_ssu.mseq//g' > otu_tables/mapseq_counts_genus_rev_bac_ssu.tsv
		"""
	} else {
		"""
		mkdir -p otu_tables
		${params.mapseq_bin} -otutable -tl 5 \$(ls *_R1_bac_ssu.mseq) | sed 's/_R1_bac_ssu.mseq//g' > otu_tables/mapseq_counts_genus_fwd_bac_ssu.tsv
		"""
	}
}


process motus {
	publishDir "$output_dir", mode: params.publish_mode

	input:
	tuple val(sample), path(reads)

	output:
	tuple val(sample), path("${sample}/${sample}.motus.txt"), emit: motus_out

	script:
	def motus_input = (reads.size() == 2) ? "-f ${sample}_R1.fastq.gz -r ${sample}_R2.fastq.gz" : "-s ${sample}_R1.fastq.gz";
	"""
	mkdir -p ${sample}
	motus profile -t $task.cpus -g 1 -k genus -c -v 7 ${motus_database} ${motus_input} > ${sample}/${sample}.motus.txt
	"""
}


process pathseq {
	publishDir "$output_dir", mode: params.publish_mode

	input:
	tuple val(sample), path(bam)

	output:
	val(sample), emit: sample
	path("${sample}/${sample}.pathseq.bam"), emit: bam
	path("${sample}/${sample}.pathseq.bam.sbi"), emit: sbi
	path("${sample}/${sample}.pathseq.score_metrics"), emit: score_metrics
	path("${sample}/${sample}.pathseq.scores"), emit: txt

	script:
	"""
	mkdir -p ${sample}
	maxmem=\$(echo \"$task.memory\"| sed 's/ GB/g/g')
	gatk --java-options \"-Xmx\$maxmem\" PathSeqPipelineSpark \\
		--input $bam \\
		--filter-bwa-image ${params.pathseq_database}/reference.fasta.img \\
		--kmer-file ${params.pathseq_database}/host.hss \\
		--min-clipped-read-length 31 \\
		--microbe-fasta ${params.pathseq_database}/microbe.fasta \\
		--microbe-bwa-image ${params.pathseq_database}/microbe.fasta.img \\
		--taxonomy-file ${params.pathseq_database}/microbe.db \\
		--output ${sample}/${sample}.pathseq.bam \\
		--scores-output ${sample}/${sample}.pathseq.scores \\
		--score-metrics ${sample}/${sample}.pathseq.score_metrics
	"""
}


workflow {

	/*
		Collect all fastq files from the input directory tree.
		It is recommended to use one subdirectory per sample.
	*/

	fastq_ch = Channel
		.fromPath(params.input_dir + "/" + "**.{fastq,fq,fastq.gz,fq.gz}")
		.map { file ->
				def sample = file.name.replaceAll(/.(fastq|fq)(.gz)?$/, "")
				sample = sample.replaceAll(/_R?[12]$/, "")
				return tuple(sample, file)
		}
		.groupTuple(sort: true)

	/*
		Normalise the input fastq naming scheme
		r1: <sample>_R1.fastq.gz
		r2: <sample>_R2.fastq.gz
	*/

	prepare_fastqs(fastq_ch)

	/*
		Collect all bam files from the input directory tree.
	*/

	bam_ch = Channel
		.fromPath(params.input_dir + "/" + "**.bam")
		.map { file ->
			def sample = file.name.replaceAll(/.bam$/, "")
			return tuple(sample, file)
		}
		.groupTuple(sort: true)

	/*
		Convert input bam to fastq.
		This gets rid of:
		- alignment information (e.g. when reads were prealigned against host)
		- optical duplicates, secondary/suppl alignments, reads that don't pass qual filters
	*/

	bam2fq(bam_ch)

	/*
		Combine the normalised and bam-extracted fastqs
	*/

	combined_fastq_ch = prepare_fastqs.out.reads.concat(bam2fq.out.reads)

	/*
		Convert all fastqs to bam files
	*/

	fq2bam(combined_fastq_ch)

	combined_bam_ch = fq2bam.out.reads

	/* perform bam-based analyses */

	count_reads(combined_bam_ch)

	pathseq(combined_bam_ch)

	/* perform fastq-based analyses */

	kraken2(combined_fastq_ch)

	motus(combined_fastq_ch)

	mtag_extraction(combined_fastq_ch)

	mapseq(mtag_extraction.out.mtag_out)

	collate_mapseq_tables(mapseq.out.bac_ssu.collect())
}
