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
	params.output_dir = "vknight_out"
}

output_dir = "${params.output_dir}"

suffix_pattern = params.file_suffix


/*
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
*/

process bam2fq {
	input:
    tuple val(sample), path(bam)

	output:
	stdout
	val(sample), emit: sample
	path("out/${sample}_R1.fastq.gz"), emit: r1
	path("out/${sample}_R2.fastq.gz"), optional: true, emit: r2
	path("out/${sample}_singles.fastq.gz"), optional: true, emit: singles

	script:
	"""
	mkdir -p out
	samtools fastq -@ $task.cpus -1 out/${sample}_R1.fastq.gz -2 out/${sample}_R2.fastq.gz -s out/${sample}_singles.fastq.gz $bam
	"""	
}


process fq2bam {
	input:
	tuple val(sample), path(fq)

	output:
	stdout
	val(sample), emit: sample
	path("out/${sample}.bam"), emit: bam

	script:
	//def r1 = "in=" + reads[r1_index] + " out=${sample}/${sample}.qc_R1.fastq.gz"
    //def r2 = (reads[1] && file(reads[1])) ? "in2=" + ( reads[ r1_index == 0 ? 1 : 0 ] + " out2=${sample}/${sample}.qc_R2.fastq.gz outs=${sample}/${sample}.qc_S1.fastq.gz" ) : ""

	//def 

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


process make_dummy_fastqs {
	input:
	tuple val(sample), path(fq)

	output:
	val(sample), emit: sample
	path("out/${sample}_R1.fastq.gz"), emit: r1
	path("out/${sample}_R2.fastq.gz"), optional: true, emit: r2

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


process make_dummy_bam {
	input:
	tuple val(sample), path(bam)

	output:
	val(sample), emit: sample
	path("out/${sample}.bam"), emit: bam

	script:
	"""
	mkdir -p out
	ln -sf ../${bam} out/${sample}.bam
	"""
}


process count_reads {
	publishDir "$output_dir", mode: params.publish_mode

	input:
	val(sample)
	path(bam)

	output:
	val(sample), emit: sample
	path("${sample}/${sample}.libsize.txt")

	script:
	"""
	mkdir -p ${sample}
	samtools flagstat $bam | head -n 1 | awk '{print \$1 + \$3}' > "${sample}/${sample}.libsize.txt"
	"""

}


process kraken2_single {
	publishDir "$output_dir", mode: params.publish_mode

	input:
	val(sample)
	path(r1)

	output:
	val(sample), emit: sample
	path("${sample}/${sample}.kraken2_report.txt"), emit: report

	script:
	"""
	mkdir -p ${sample}
	/g/scb/zeller/fspringe/Software/kraken2/kraken2 --db /g/scb/zeller/jawirbel/total_RNAseq/databases/kraken2_standard --threads $task.cpus --gzip-compressed --report ${sample}/${sample}.kraken2_report.txt $r1
	""" 
	
}

process kraken2_paired {
	publishDir "$output_dir", mode: params.publish_mode

	input:
	val(sample)
	path(r1)
	path(r2)

	output:
	val(sample), emit: sample
	path("${sample}/${sample}.kraken2_report.txt"), emit: report

	script:
	"""
	mkdir -p ${sample}
	/g/scb/zeller/fspringe/Software/kraken2/kraken2 --db /g/scb/zeller/jawirbel/total_RNAseq/databases/kraken2_standard --threads $task.cpus --gzip-compressed --report ${sample}/${sample}.kraken2_report.txt --paired $r1 $r2
	"""

}

process mtag_extraction_single {
	input:
	val(sample)
	path(r1)

	output:
	val(sample), emit: sample
	//path("${r1}_arc_lsu.fasta"), emit: arc_lsu_r1
	//path("${r1}_arc_ssu.fasta"), emit: arc_ssu_r1
	path("${r1}_bac_lsu.fasta"), emit: bac_lsu_r1
	path("${r1}_bac_ssu.fasta"), emit: bac_ssu_r1
	//path("${r1}_euk_lsu.fasta"), emit: euk_lsu_r1
	//path("${r1}_euk_ssu.fasta"), emit: euk_ssu_r1
	//path("${r1}_read.map"), emit: readmap_r1

	script:
	"""
	mtags extract -is $r1 -o . -t $task.cpus
	"""
}


process mtag_extraction_paired {
	input:
	val(sample)
	path(r1)
	path(r2)

	output:
	val(sample), emit: sample
	//path("${r1}_arc_lsu.fasta"), emit: arc_lsu_r1
	//path("${r1}_arc_ssu.fasta"), emit: arc_ssu_r1
	path("${r1}_bac_lsu.fasta"), emit: bac_lsu_r1
	path("${r1}_bac_ssu.fasta"), emit: bac_ssu_r1
	//path("${r1}_euk_lsu.fasta"), emit: euk_lsu_r1
	//path("${r1}_euk_ssu.fasta"), emit: euk_ssu_r1
	//path("${r1}_read.map"), emit: readmap_r1
	//path("${r2}_arc_lsu.fasta"), emit: arc_lsu_r2
	//path("${r2}_arc_ssu.fasta"), emit: arc_ssu_r2
	path("${r2}_bac_lsu.fasta"), emit: bac_lsu_r2
	path("${r2}_bac_ssu.fasta"), emit: bac_ssu_r2
	//path("${r2}_euk_lsu.fasta"), emit: euk_lsu_r2
	//path("${r2}_euk_ssu.fasta"), emit: euk_ssu_r2
	//path("${r2}_read.map"), emit: readmap_r2
	script:
	"""
	mtags extract -i1 $r1 -i2 $r2 -o . -t $task.cpus      
	"""
}


process mapseq_single {
	input:
	val(sample)
    path(bac_lsu_r1)
    path(bac_ssu_r1)

	output:
	val(sample), emit: sample
	path("${sample}/${sample}_R1_bac_lsu.mseq"), emit: bac_lsu_r1
	path("${sample}/${sample}_R1_bac_ssu.mseq"), emit: bac_ssu_r1

	script:
	"""
	mkdir -p ${sample}
	/g/scb/zeller/fspringe/Software/mapseq-1.2.6-linux/mapseq $bac_lsu_r1 > ${sample}/${sample}_R1_bac_lsu.mseq
	/g/scb/zeller/fspringe/Software/mapseq-1.2.6-linux/mapseq $bac_ssu_r1 > ${sample}/${sample}_R1_bac_ssu.mseq
	"""	
}


process mapseq_paired {
	input:
	val(sample)
    path(bac_lsu_r1)
    path(bac_ssu_r1)
    path(bac_lsu_r2)
    path(bac_ssu_r2)

	output:
	val(sample), emit: sample
	path("${sample}/${sample}_R1_bac_lsu.mseq"), emit: bac_lsu_r1
	path("${sample}/${sample}_R1_bac_ssu.mseq"), emit: bac_ssu_r1
	path("${sample}/${sample}_R2_bac_lsu.mseq"), emit: bac_lsu_r2
	path("${sample}/${sample}_R2_bac_ssu.mseq"), emit: bac_ssu_r2

	script:
	"""
	mkdir -p ${sample}
	/g/scb/zeller/fspringe/Software/mapseq-1.2.6-linux/mapseq $bac_lsu_r1 > ${sample}/${sample}_R1_bac_lsu.mseq
	/g/scb/zeller/fspringe/Software/mapseq-1.2.6-linux/mapseq $bac_ssu_r1 > ${sample}/${sample}_R1_bac_ssu.mseq
	/g/scb/zeller/fspringe/Software/mapseq-1.2.6-linux/mapseq $bac_lsu_r2 > ${sample}/${sample}_R2_bac_lsu.mseq
	/g/scb/zeller/fspringe/Software/mapseq-1.2.6-linux/mapseq $bac_ssu_r2 > ${sample}/${sample}_R2_bac_ssu.mseq
	"""	
}


process collate_mapseq_single {
	publishDir "$output_dir", mode: params.publish_mode

	input:
    path(bac_lsu_r1)
    path(bac_ssu_r1)

	output:
	path("otu_tables/mapseq_counts_genus_fwd_bac_ssu.tsv")

	script:
	"""
	mkdir -p otu_tables
	/g/scb/zeller/fspringe/Software/mapseq-1.2.6-linux/mapseq -otutable -tl 5 $bac_ssu_r1 > otu_tables/mapseq_counts_genus_fwd_bac_ssu.tsv
	"""
}


process collate_mapseq_paired {
	publishDir "$output_dir", mode: params.publish_mode

	input:
    path(bac_lsu_r1)
    path(bac_ssu_r1)
    path(bac_lsu_r2)
    path(bac_ssu_r2)

	output:
	path("otu_tables/mapseq_counts_genus_fwd_bac_ssu.tsv")
	path("otu_tables/mapseq_counts_genus_rev_bac_ssu.tsv")

	script:
	"""
	mkdir -p otu_tables
	/g/scb/zeller/fspringe/Software/mapseq-1.2.6-linux/mapseq -otutable -tl 5 $bac_ssu_r1 > otu_tables/mapseq_counts_genus_fwd_bac_ssu.tsv
	/g/scb/zeller/fspringe/Software/mapseq-1.2.6-linux/mapseq -otutable -tl 5 $bac_ssu_r2 > otu_tables/mapseq_counts_genus_rev_bac_ssu.tsv
	"""
}



process motus_single {
	publishDir "$output_dir", mode: params.publish_mode

	input:
	val(sample)
	path(r1)

	output:
	val(sample), emit: sample
	path("${sample}/${sample}.motus.txt"), emit: motus_out

	script:
	"""
	mkdir -p ${sample}
    motus profile -s $r1 -l 50 -t $task.cpus -g 1 -k genus -c -v 7 > ${sample}/${sample}.motus.txt
	"""
}


process motus_paired {
	publishDir "$output_dir", mode: params.publish_mode

	input:
	val(sample)
	path(r1)
	path(r2)

	output:
	val(sample), emit: sample
	path("${sample}/${sample}.motus.txt"), emit: motus_out

	script:
	"""
	mkdir -p ${sample}
    motus profile -f $r1 -r $r2 -l 50 -t $task.cpus -g 1 -k genus -c -v 7 > ${sample}/${sample}.motus.txt
	"""
}


process pathseq {
	publishDir "$output_dir", mode: params.publish_mode

	input:
	val(sample)
	path(bam)

	output:
	val(sample), emit: sample
	path("${sample}/${sample}.pathseq.bam"), emit: bam
	path("${sample}/${sample}.pathseq.bam.sbi"), emit: sbi
	path("${sample}/${sample}.pathseq.score_metrics"), emit: score_metrics
	path("${sample}/${sample}.pathseq.scores"), emit: txt


	//maxmem="-Xmx"$(echo "256 GB"| sed 's/ GB/g/g')
	//echo $maxmem
	//gatk --java-options "-Xmx$maxmem" PathSeqPipelineSpark

	script:
	"""
	mkdir -p ${sample}
	maxmem=\$(echo \"$task.memory\"| sed 's/ GB/g/g')
	gatk --java-options \"-Xmx\$maxmem\" PathSeqPipelineSpark \\
		--input $bam \\
		--filter-bwa-image ${params.pathseq_database}/reference.fasta.img \\
		--kmer-file ${params.pathseq_database}/host.hss \\
		--min-clipped-read-length 70 \\
		--microbe-fasta ${params.pathseq_database}/microbe.fasta \\
		--microbe-bwa-image ${params.pathseq_database}/microbe.fasta.img \\
		--taxonomy-file ${params.pathseq_database}/microbe.db \\
		--output ${sample}/${sample}.pathseq.bam \\
		--scores-output ${sample}/${sample}.pathseq.scores \\
		--score-metrics ${sample}/${sample}.pathseq.score_metrics
	"""
}



workflow {
	fastq_ch = Channel
    	.fromPath(params.input_dir + "/" + "**.{fastq,fq,fastq.gz,fq.gz}")
        .map { file ->
				def sample = file.name.replaceAll(/.(fastq|fq)(.gz)?$/, "")
				sample = sample.replaceAll(/_R?[12]$/, "")
				return tuple(sample, file)
		}
		.groupTuple()
	fastq_ch.view()

	bam_ch = Channel
		.fromPath(params.input_dir + "/" + "**.bam")
		.map { file ->
			def sample = file.name.replaceAll(/.bam$/, "")
			return tuple(sample, file)
		}
		.groupTuple(sort: true)
	bam_ch.view()

	make_dummy_fastqs(fastq_ch)
	make_dummy_bam(bam_ch)
	
	bam2fq(bam_ch)
	fq2bam(fastq_ch)

	if (params.bam_input) {
		count_reads(make_dummy_bam.out.sample, make_dummy_bam.out.bam)
		pathseq(make_dummy_bam.out.sample, make_dummy_bam.out.bam)
		if (bam2fq.out.r2 != null) {
			kraken2_paired(bam2fq.out.sample, bam2fq.out.r1, bam2fq.out.r2)
			motus_paired(bam2fq.out.sample, bam2fq.out.r1, bam2fq.out.r2)
			mtag_extraction_paired(bam2fq.out.sample, bam2fq.out.r1, bam2fq.out.r2)
			mapseq_paired(mtag_extraction_paired.out.sample, mtag_extraction_paired.out.bac_lsu_r1, mtag_extraction_paired.out.bac_ssu_r1, mtag_extraction_paired.out.bac_lsu_r2, mtag_extraction_paired.out.bac_ssu_r2)
			collate_mapseq_paired(mapseq_paired.out.bac_lsu_r1.collect(), mapseq_paired.out.bac_ssu_r1.collect(), mapseq_paired.out.bac_lsu_r2.collect(), mapseq_paired.out.bac_ssu_r2.collect())
		} else {
			kraken2_single(bam2fq.out.sample, bam2fq.out.r1)
			motus_single(bam2fq.out.sample, bam2fq.out.r1)
			mtag_extraction_single(bam2fq.out.sample, bam2fq.out.r1)
			mapseq_single(mtag_extraction_single.out.sample, mtag_extraction_single.out.bac_lsu_r1, mtag_extraction_single.out.bac_ssu_r1)
			collate_mapseq_single(mapseq_single.out.bac_lsu_r1.collect(), mapseq_single.out.bac_ssu_r1.collect())
		}
	} else {
		count_reads(fq2bam.out.sample, fq2bam.out.bam)
		pathseq(fq2bam.out.sample, fq2bam.out.bam)
		if (make_dummy_fastqs.out.r2 != null) {
			//kraken2_paired(make_dummy_fastqs.out.sample, make_dummy_fastqs.out.r1, make_dummy_fastqs.out.r2)
			//motus_paired(make_dummy_fastqs.out.sample, make_dummy_fastqs.out.r1, make_dummy_fastqs.out.r2)
            mtag_extraction_paired(make_dummy_fastqs.out.sample, make_dummy_fastqs.out.r1, make_dummy_fastqs.out.r2)
			mapseq_paired(mtag_extraction_paired.out.sample, mtag_extraction_paired.out.bac_lsu_r1, mtag_extraction_paired.out.bac_ssu_r1, mtag_extraction_paired.out.bac_lsu_r2, mtag_extraction_paired.out.bac_ssu_r2)
			collate_mapseq_paired(mapseq_paired.out.bac_lsu_r1.collect(), mapseq_paired.out.bac_ssu_r1.collect(), mapseq_paired.out.bac_lsu_r2.collect(), mapseq_paired.out.bac_ssu_r2.collect())
		} else {
			kraken2_single(make_dummy_fastqs.out.sample, make_dummy_fastqs.out.r1)
			motus_single(make_dummy_fastqs.out.sample, make_dummy_fastqs.out.r1)
			mtag_extraction_single(make_dummy_fastqs.out.sample, make_dummy_fastqs.out.r1)
			mapseq_single(mtag_extraction_single.out.sample, mtag_extraction_single.out.bac_lsu_r1, mtag_extraction_single.out.bac_ssu_r1)
			collate_mapseq_single(mapseq_single.out.bac_lsu_r1.collect(), mapseq_single.out.bac_ssu_r1.collect())
		}
	}


}
