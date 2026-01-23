process starmap {
	cpus 2
	memory { 8.GB * task.attempt }
	time { 3.d * task.attempt }
	container "quay.io/biocontainers/star:2.7.0b--0"
	publishDir "${params.output_dir}/host_gene_counts", mode: "copy", pattern: "${sample.id}/*.tsv"

    input:
    tuple val(sample), path(reads)
    path(star_index)

    output:
    // tuple val(sample), path("${sample}/*.fastq.gz"), emit: reads
	tuple val(sample), path("${sample.id}/${sample.id}.bam"), emit: bam
	tuple val(sample), path("${sample.id}/*.tsv"), emit: gene_counts

    script:
    """
	set -e -o pipefail

    mkdir -p ${sample.id}

    STAR --runThreadN ${task.cpus} --genomeDir ${star_index} --readFilesIn ${reads} --readFilesCommand zcat --outFileNamePrefix ./${sample.id}. --outSAMunmapped Within KeepPairs --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts

	mv -v ${sample.id}.Aligned.sortedByCoord.out.bam ${sample.id}/${sample.id}.bam
	mv -v ${sample.id}.ReadsPerGene.out.tab ${sample.id}/${sample.id}.ReadsPerGene.tsv
    """

}