process starmap {
	cpus 2
	memory { 8.GB * task.attempt }
	time { 3.d * task.attempt }
	// container "dceoy/star:latest"
	container "quay.io/biocontainers/star:2.7.0b--0"
	publishDir "${params.output_dir}/host_gene_counts", mode: "copy", pattern: "${sample}/*.tsv"

    input:
    tuple val(sample), path(reads)
    path(star_index)

    output:
    // tuple val(sample), path("${sample}/*.fastq.gz"), emit: reads
	tuple val(sample), path("${sample}/${sample}.bam"), emit: bam
	tuple val(sample), path("${sample}/*.tsv"), emit: gene_counts

    script:
    """
	set -e -o pipefail

    mkdir -p ${sample}

    STAR --runThreadN ${task.cpus} --genomeDir ${star_index} --readFilesIn ${reads} --readFilesCommand zcat --outFileNamePrefix ./${sample}. --outReadsUnmapped Fastx --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts

	mv -v ${sample}.Aligned.sortedByCoord.out.bam ${sample}/${sample}.bam
	mv -v ${sample}.ReadsPerGene.out.tab ${sample}/${sample}.ReadsPerGene.tsv
    """

}