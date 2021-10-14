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

if (!params.motus2_min_length) {
	params.motus2_min_length = 30
}

if (!params.motus2_n_marker_genes) {
	params.motus2_n_marker_genes = 1
}

if (!params.pathseq_min_clipped_read_length) {
	params.pathseq_min_clipped_read_length = 31
}

if (!params.kraken2_min_hit_groups) {
	params.kraken2_min_hit_groups = 10
}

def run_kraken2 = (!params.skip_kraken2 || params.run_kraken2);
def run_mtags = (!params.skip_mtags || params.run_mtags);
def run_mapseq = (run_mtags && (!params.skip_mapseq || params.run_mapseq) && params.mapseq_bin)
def run_motus2 = (!params.skip_motus2 || params.run_motus2);
def run_pathseq = (!params.skip_pathseq || params.run_pathseq);
def run_count_reads = (!params.skip_counts || params.run_counts);
def convert_fastq2bam = (run_pathseq || run_count_reads);

def do_preprocessing = (!params.skip_preprocessing || params.run_preprocessing)

def run_bam_analysis = run_pathseq //|| run_count_reads
def run_fastq_analysis = run_kraken2 || run_mtags || run_mapseq || run_motus2

include { nevermore_preprocess } from "./workflows/nevermore/nevermore"
include { bam_analysis; fastq_analysis; collate_data } from "./workflows/vknight/vknight"


process bam2fq {
	publishDir "$output_dir", mode: params.publish_mode

	input:
	tuple val(sample), path(bam)

	output:
	stdout
	tuple val(sample), path("fastq/${sample}/${sample}*.fastq.gz"), emit: reads

	script:
	"""
	set -o pipefail
	mkdir -p fastq/${sample}
	samtools collate -@ $task.cpus -u -O $bam | samtools fastq -F 0x900 -0 ${sample}_other.fastq.gz -1 ${sample}_R1.fastq.gz -2 ${sample}_R2.fastq.gz

	if [[ "\$?" -eq 0 ]];
	then

		if [[ -z "\$(gzip -dc ${sample}_R1.fastq.gz | head -n 1)" ]];
		then
			if [[ ! -z "\$(gzip -dc ${sample}_other.fastq.gz | head -n 1)" ]];
			then
				mv -v ${sample}_other.fastq.gz fastq/${sample}/${sample}_R1.fastq.gz;
			fi;
		else
				mv -v ${sample}_R1.fastq.gz fastq/${sample}/;
				if [[ ! -z "\$(gzip -dc ${sample}_R2.fastq.gz | head -n 1)" ]];
				then
					mv -v ${sample}_R2.fastq.gz fastq/${sample}/;
				fi;
		fi;

		ls -l *.fastq.gz
		ls -l fastq/${sample}/*.fastq.gz
		rm -rf *.fastq.gz
	fi;
	"""
}


process prepare_fastqs {
	publishDir "$output_dir", mode: params.publish_mode

	input:
	tuple val(sample), path(fq)

	output:
	tuple val(sample), path("fastq/${sample}/${sample}_R*.fastq.gz"), emit: reads

	script:
	if (fq.size() == 2) {
		"""
		mkdir -p fastq/${sample}
		ln -sf ../../${fq[0]} fastq/${sample}/${sample}_R1.fastq.gz
		ln -sf ../../${fq[1]} fastq/${sample}/${sample}_R2.fastq.gz
		"""
	} else {
		"""
		mkdir -p fastq/${sample}
		ln -sf ../../${fq[0]} fastq/${sample}/${sample}_R1.fastq.gz
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


process remove_host_kraken2 {
	publishDir "$output_dir", mode: params.publish_mode

	input:
	tuple val(sample), path(fq)

	output:
	tuple val(sample), path("no_host/${sample}/${sample}_R*.fastq.gz"), emit: reads

	script:
	def out_options = (fq.size() == 2) ? "--paired --unclassified-out ${sample}#.fastq" : "--unclassified-out ${sample}_1.fastq"
	def move_r2 = (fq.size() == 2) ? "gzip -c ${sample}_2.fastq > no_host/${sample}/${sample}_R2.fastq.gz" : ""

	"""
	mkdir -p no_host/${sample}
	kraken2 --threads $task.cpus --db ${params.remove_host_kraken2_db} ${out_options} --output kraken_read_report.txt --report kraken_report.txt --report-minimizer-data --gzip-compressed --minimum-hit-groups ${params.kraken2_min_hit_groups} $fq

	gzip -c ${sample}_1.fastq > no_host/${sample}/${sample}_R1.fastq.gz
	${move_r2}
	"""
}


process flagstats {
    publishDir "$output_dir", mode: params.publish_mode

    input:
    tuple val(sample), path(bam)

    output:
    tuple val(sample), path("${sample}/${sample}.flagstats.txt"), emit: flagstats

    script:
    """
    mkdir -p ${sample}
    samtools flagstat $bam > "${sample}/${sample}.flagstats.txt"
    """
}


process count_reads {
    publishDir "$output_dir", mode: params.publish_mode

    input:
    //tuple val(sample), path(bam)
    tuple val(sample), path(flagstats)

    output:
    tuple val(sample), path("${sample}/${sample}.libsize.txt"), emit: counts
    tuple val(sample), path("${sample}/${sample}.is_paired.txt"), emit: is_paired
    //tuple val(sample), path("${sample}/${sample}.flagstats.txt"), emit: flagstats

    script:
    """
    mkdir -p ${sample}
    head -n 1 "${flagstats[0]}" | awk '{print \$1 + \$3}' > "${sample}/${sample}.libsize.txt"
    grep -m 1 "paired in sequencing" "${flagstats[0]}" | awk '{npaired = \$1 + \$3; if (npaired==0) {print "unpaired"} else {print "paired"};}' > "${sample}/${sample}.is_paired.txt"
    """
    //samtools flagstat $bam > "${sample}/${sample}.flagstats.txt"
    //head -n 1 "${sample}/${sample}.flagstats.txt" | awk '{print \$1 + \$3}' > "${sample}/${sample}.libsize.txt"
    //grep -m 1 "paired in sequencing" "${sample}/${sample}.flagstats.txt" | awk '{npaired = \$1 + \$3; if (npaired==0) {print "unpaired"} else {print "paired"};}' > "${sample}/${sample}.is_paired.txt"
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

	bam_ch = Channel
		.fromPath(params.input_dir + "/" + "**.bam")
		.map { file ->
			def sample = file.name.replaceAll(/.bam$/, "")
			return tuple(sample, file)
		}
		.groupTuple(sort: true)

	bam2fq(bam_ch)

	if (do_preprocessing) {

		prepare_fastqs(fastq_ch)

		raw_fastq_ch = prepare_fastqs.out.reads.concat(bam2fq.out.reads)

		nevermore_preprocess(raw_fastq_ch)


		if (params.remove_host) {

			remove_host_kraken2(nevermore_preprocess.out.preprocessed)

			preprocessed_ch = remove_host_kraken2.out.reads

		} else {

			preprocessed_ch = nevermore_preprocess.out.preprocessed

		}

	} else {

		preprocessed_ch = fastq_ch

	}


	if (run_count_reads || run_bam_analysis) {

		fq2bam(preprocessed_ch)

		if (run_count_reads) {

	        flagstats(fq2bam.out.reads)

    	    count_reads(flagstats.out.flagstats)

		}

		if (run_bam_analysis) {

			bam_analysis(fq2bam.out.reads)

		}

    }


	if (run_fastq_analysis) {

		fastq_analysis(preprocessed_ch)

	}


//	/* collate data */
//
//	if (params.collate_script != null && params.collate_script != "") {
//		data_to_collate_ch = Channel.empty()
//
//		if (run_kraken2) {
//			data_to_collate_ch = data_to_collate_ch.concat(kraken2.out.kraken2_out)
//		}
//
//		if (run_count_reads) {
//			data_to_collate_ch = data_to_collate_ch.concat(count_reads.out.counts)
//				.concat(count_reads.out.is_paired)
//		}
//
//		if (run_motus2) {
//			data_to_collate_ch = data_to_collate_ch.concat(motus2.out.motus_out)
//		}
//
//		if (run_pathseq) {
//			data_to_collate_ch = data_to_collate_ch.concat(pathseq.out.scores)
//		}
//
//		data_to_collate_ch = data_to_collate_ch
//			.map { sample, files -> return files }
//
//		if (run_mtags) {
//			data_to_collate_ch = data_to_collate_ch.concat(mtags_merge.out.mtags_tables)
//		}
//
//		collate_data(data_to_collate_ch.collect())
//	}
}
