process mtags_extract {
    container "ghcr.io/cschu/vknight_profilers:main"
    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("*_bac_ssu.fasta"), optional: true, emit: mtags_out

    script:
    def mtags_input = (sample.is_paired) ? "-f ${sample.id}_R1.fastq.gz -r ${sample.id}_R2.fastq.gz" : "-s ${sample.id}_R1.fastq.gz";

    """
    mtags extract -t $task.cpus -o . ${mtags_input}
    """
}


process mtags_annotate {
    container "ghcr.io/cschu/vknight_profilers:main"
    input:
    tuple val(sample), path(seqs)

    output:
    path("${sample.id}.bins"), emit: mtags_bins

    script:
    def mtags_input = (seqs.size() == 2) ? "-f ${sample.id}_R1.fastq.gz_bac_ssu.fasta -r ${sample.id}_R2.fastq.gz_bac_ssu.fasta" : "-s ${sample.id}_R1.fastq.gz_bac_ssu.fasta";
    """
    mtags annotate -t $task.cpus -o . -n ${sample.id} ${mtags_input}
    """
}


process mtags_merge {
    container "ghcr.io/cschu/vknight_profilers:main"
    publishDir params.output_dir, mode: "copy"

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
