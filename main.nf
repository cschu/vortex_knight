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

def run_kraken2 = (!params.skip_kraken2 || params.run_kraken2);
def run_mtags = (!params.skip_mtags || params.run_mtags);
def run_mapseq = (run_mtags && (!params.skip_mapseq || params.run_mapseq) && params.mapseq_bin)
def run_motus2 = (!params.skip_motus2 || params.run_motus2);
def run_pathseq = (!params.skip_pathseq || params.run_pathseq);
def run_count_reads = (!params.skip_counts || params.run_counts);
def convert_fastq2bam = (run_pathseq || run_count_reads);


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


process count_reads {
	publishDir "$output_dir", mode: params.publish_mode

	input:
	tuple val(sample), path(bam)

	output:
	tuple val(sample), path("${sample}/${sample}.libsize.txt"), emit: counts
	tuple val(sample), path("${sample}/${sample}.is_paired.txt"), emit: is_paired
	tuple val(sample), path("${sample}/${sample}.flagstats.txt"), emit: flagstats

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


process mtags_extract {
	input:
	tuple val(sample), path(reads)

	output:
	tuple val(sample), path("*_bac_ssu.fasta"), emit: mtags_out

	script:
	def mtags_input = (reads.size() == 2) ? "-f ${sample}_R1.fastq.gz -r ${sample}_R2.fastq.gz" : "-s ${sample}_R1.fastq.gz";

	"""
	mtags extract -t $task.cpus -o . ${mtags_input}
	"""
}


process mtags_annotate {
	input:
	tuple val(sample), path(seqs)

	output:
	path("${sample}.bins"), emit: mtags_bins

	script:
	def mtags_input = (seqs.size() == 2) ? "-f ${sample}_R1.fastq.gz_bac_ssu.fasta -r ${sample}_R2.fastq.gz_bac_ssu.fasta" : "-s ${sample}_R1.fastq.gz_bac_ssu.fasta";
	"""
	mtags annotate -t $task.cpus -o . -n ${sample} ${mtags_input}
	"""
}


process mtags_merge {
	publishDir "$output_dir", mode: params.publish_mode

	input:
	path(mtags_bins)

	output:
	path("mtags_tables/merged_profile.*"), emit: mtags_tables

	script:
	"""
	mkdir mtags_tables
	mtags merge -i *bins -o mtags_tables/merged_profile
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
	path("mapseq_tables/mapseq_counts_genus_*_bac_ssu.tsv"), emit: ssu_tables

	script:
	if (mapped_seqs.size() == 2) {
		"""
		mkdir -p mapseq_tables
		${params.mapseq_bin} -otutable -tl 5 \$(ls *_R1_bac_ssu.mseq) | sed 's/_R1_bac_ssu.mseq//g' > mapseq_tables/mapseq_counts_genus_fwd_bac_ssu.tsv
		${params.mapseq_bin} -otutable -tl 5 \$(ls *_R2_bac_ssu.mseq) | sed 's/_R2_bac_ssu.mseq//g' > mapseq_tables/mapseq_counts_genus_rev_bac_ssu.tsv
		"""
	} else {
		"""
		mkdir -p mapseq_tables
		${params.mapseq_bin} -otutable -tl 5 \$(ls *_R1_bac_ssu.mseq) | sed 's/_R1_bac_ssu.mseq//g' > mapseq_tables/mapseq_counts_genus_fwd_bac_ssu.tsv
		"""
	}
}



process motus2 {
	publishDir "$output_dir", mode: params.publish_mode

	input:
	tuple val(sample), path(reads)

	output:
	tuple val(sample), path("${sample}/${sample}.motus.txt"), emit: motus_out

	script:
	def motus_input = (reads.size() == 2) ? "-f ${sample}_R1.fastq.gz -r ${sample}_R2.fastq.gz" : "-s ${sample}_R1.fastq.gz";
	"""
	mkdir -p ${sample}
	motus profile -t $task.cpus -k genus -c -v 7 -l ${params.motus2_min_length} -g ${params.motus2_n_marker_genes} ${motus_database} ${motus_input} > ${sample}/${sample}.motus.txt
	"""
}


process pathseq {
	publishDir "$output_dir", mode: params.publish_mode

	input:
	tuple val(sample), path(bam)

	output:
	tuple val(sample), path("${sample}/${sample}.pathseq.score*"), emit: scores
	tuple val(sample), path("${sample}/${sample}.pathseq.bam*"), emit: bam

	script:
	"""
	mkdir -p ${sample}
	maxmem=\$(echo \"$task.memory\"| sed 's/ GB/g/g')
	gatk --java-options \"-Xmx\$maxmem\" PathSeqPipelineSpark \\
		--input $bam \\
		--filter-bwa-image ${params.pathseq_database}/reference.fasta.img \\
		--kmer-file ${params.pathseq_database}/host.hss \\
		--min-clipped-read-length ${params.pathseq_min_clipped_read_length} \\
		--microbe-fasta ${params.pathseq_database}/microbe.fasta \\
		--microbe-bwa-image ${params.pathseq_database}/microbe.fasta.img \\
		--taxonomy-file ${params.pathseq_database}/microbe.db \\
		--output ${sample}/${sample}.pathseq.bam \\
		--scores-output ${sample}/${sample}.pathseq.scores \\
		--score-metrics ${sample}/${sample}.pathseq.score_metrics
	"""
}


process collate_data {
	publishDir "$output_dir", mode: params.publish_mode

	input:
	path(outputs)

	output:
	path("combined_profiles_rds/*.rds"), emit: collated_data

	script:
	def kraken2_output = run_kraken2 ? "--kraken2_res_path kraken2/" : ""
	def motus_output = run_motus2 ? "--mOTUs_res_path motus/" : ""
	def mtags_output = run_mtags ? "--mTAGs_res_path mtags_tables/" : ""
	def counts_output = run_count_reads ? "--libsize_res_path libsize/ --lib_layout_res_path lib_layout/" : ""
	def pathseq_output = run_pathseq ? "--PathSeq_res_path pathseq/" : ""
	"""
	mkdir -p combined_profiles_rds/
	mkdir -p {kraken2,pathseq,motus,libsize,lib_layout,otu_tables,mtags_tables}

	find \$(pwd) -maxdepth 1 -name '*kraken2_report.txt' -exec ln -sf {} kraken2/ \\;
	find \$(pwd) -maxdepth 1 -name '*pathseq.scores' -exec ln -sf {} pathseq/ \\;
	find \$(pwd) -maxdepth 1 -name '*motus.txt' -exec ln -sf {} motus/ \\;
	find \$(pwd) -maxdepth 1 -name '*libsize.txt' -exec ln -sf {} libsize/ \\;
	find \$(pwd) -maxdepth 1 -name '*is_paired.txt' -exec ln -sf {} lib_layout/ \\;
	find \$(pwd) -maxdepth 1 -name '*bac_ssu.tsv' -exec ln -sf {} otu_tables/ \\;
	find \$(pwd) -maxdepth 1 -name 'merged_profile.genus.tsv' -exec ln -sf {} mtags_tables/ \\;

	Rscript ${params.collate_script} \\
	${kraken2_output} \\
	${motus_output} \\
	${mtags_output} \\
	${pathseq_output} \\
	${counts_output} \\
	--out_folder combined_profiles_rds/
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

	if (convert_fastq2bam) {

		/*
			Convert all fastqs to bam files
		*/

		fq2bam(combined_fastq_ch)

		combined_bam_ch = fq2bam.out.reads

		/* perform bam-based analyses */

		if (run_count_reads) {
			count_reads(combined_bam_ch)
		}

		if (run_pathseq) {
			pathseq(combined_bam_ch)
		}

	}

	/* perform fastq-based analyses */

	if (run_kraken2) {
		kraken2(combined_fastq_ch)
	}

	if (run_motus2) {
		motus2(combined_fastq_ch)
	}

	if (run_mtags) {
		mtags_extract(combined_fastq_ch)

		mtags_annotate(mtags_extract.out.mtags_out)

		mtags_merge(mtags_annotate.out.mtags_bins.collect())

		if (run_mapseq) {
			mapseq(mtags_extract.out.mtags_out)

			collate_mapseq_tables(mapseq.out.bac_ssu.collect())
		}
	}

	/* collate data */

	if (params.collate_script != null && params.collate_script != "") {
		data_to_collate_ch = Channel.empty()

		if (run_kraken2) {
			data_to_collate_ch = data_to_collate_ch.concat(kraken2.out.kraken2_out)
		}

		if (run_count_reads) {
			data_to_collate_ch = data_to_collate_ch.concat(count_reads.out.counts)
				.concat(count_reads.out.is_paired)
		}

		if (run_motus2) {
			data_to_collate_ch = data_to_collate_ch.concat(motus2.out.motus_out)
		}

		if (run_pathseq) {
			data_to_collate_ch = data_to_collate_ch.concat(pathseq.out.scores)
		}

		data_to_collate_ch = data_to_collate_ch
			.map { sample, files -> return files }

		if (run_mtags) {
			data_to_collate_ch = data_to_collate_ch.concat(mtags_merge.out.mtags_tables)
		}

		collate_data(data_to_collate_ch.collect())
	}
}
