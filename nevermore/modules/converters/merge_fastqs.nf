process merge_single_fastqs {
    input:
    tuple val(sample), path(fastqs)

    output:
    tuple val(sample), path("merged/${sample.id}_R1.fastq.gz"), emit: fastq

    script:
    """
    mkdir -p merged/

    cat *.fastq.gz | sortbyname.sh in=stdin.gz out=merged/${sample.id}_R1.fastq.gz
    """

}
