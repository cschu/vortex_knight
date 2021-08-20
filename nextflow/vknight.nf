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
	val(sample), emit: sample
	path("out/${sample}_R1.fastq.gz"), emit: r1
	path("out/${sample}_R2.fastq.gz"), optional: true, emit: r2

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


process fq2bam {
	input:
	tuple val(sample), path(fq)

	output:
	stdout
	val(sample), emit: sample
	path("out/${sample}.bam"), emit: bam

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
	kraken2 --db ${params.kraken_database} --threads $task.cpus --gzip-compressed --report ${sample}/${sample}.kraken2_report.txt $r1
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
	kraken2 --db ${params.kraken_database} --threads $task.cpus --gzip-compressed --report ${sample}/${sample}.kraken2_report.txt --paired $r1 $r2
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
	${params.mapseq_bin} $bac_lsu_r1 > ${sample}/${sample}_R1_bac_lsu.mseq
	${params.mapseq_bin} $bac_ssu_r1 > ${sample}/${sample}_R1_bac_ssu.mseq
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
	${params.mapseq_bin} $bac_lsu_r1 > ${sample}/${sample}_R1_bac_lsu.mseq
	${params.mapseq_bin} $bac_ssu_r1 > ${sample}/${sample}_R1_bac_ssu.mseq
	${params.mapseq_bin} $bac_lsu_r2 > ${sample}/${sample}_R2_bac_lsu.mseq
	${params.mapseq_bin} $bac_ssu_r2 > ${sample}/${sample}_R2_bac_ssu.mseq
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
	${params.mapseq_bin} -otutable -tl 5 \$(ls *_R1_bac_lsu.mseq) | sed 's/_R1_bac_ssu.mseq//g' > otu_tables/mapseq_counts_genus_fwd_bac_ssu.tsv
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
	${params.mapseq_bin} -otutable -tl 5 \$(ls *_R1_bac_ssu.mseq) | sed 's/_R1_bac_ssu.mseq//g' > otu_tables/mapseq_counts_genus_fwd_bac_ssu.tsv
	${params.mapseq_bin} -otutable -tl 5 \$(ls *_R2_bac_ssu.mseq) | sed 's/_R2_bac_ssu.mseq//g' > otu_tables/mapseq_counts_genus_rev_bac_ssu.tsv
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
    motus profile ${motus_database} -s $r1 -l 50 -t $task.cpus -g 1 -k genus -c -v 7 > ${sample}/${sample}.motus.txt
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
    motus profile ${motus_database} -f $r1 -r $r2 -l 50 -t $task.cpus -g 1 -k genus -c -v 7 > ${sample}/${sample}.motus.txt
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
	fastq_ch = Channel
    	.fromPath(params.input_dir + "/" + "**.{fastq,fq,fastq.gz,fq.gz}")
        .map { file ->
				def sample = file.name.replaceAll(/.(fastq|fq)(.gz)?$/, "")
				sample = sample.replaceAll(/_R?[12]$/, "")
				return tuple(sample, file)
		}
		.groupTuple(sort: true)
	//fastq_ch.view()

	bam_ch = Channel
		.fromPath(params.input_dir + "/" + "**.bam")
		.map { file ->
			def sample = file.name.replaceAll(/.bam$/, "")
			return tuple(sample, file)
		}
		.groupTuple(sort: true)
	//bam_ch.view()

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
			kraken2_paired(make_dummy_fastqs.out.sample, make_dummy_fastqs.out.r1, make_dummy_fastqs.out.r2)
			motus_paired(make_dummy_fastqs.out.sample, make_dummy_fastqs.out.r1, make_dummy_fastqs.out.r2)
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
