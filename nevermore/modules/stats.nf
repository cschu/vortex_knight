process flagstats {
    container "registry.git.embl.org/schudoma/align-docker:latest"
    label "default"

    input:
    tuple val(sample), path(bam)
    val(stage)

    output:
    tuple val(sample), path("${stage}/${sample.id}/${sample.id}.flagstats.txt"), emit: flagstats
    tuple val(sample), path("${stage}/${sample.id}/${sample.id}.libsize.txt"), emit: counts
    tuple val(sample), path("${stage}/${sample.id}/${sample.id}.is_paired.txt"), emit: is_paired
    
    script:
    """
    mkdir -p ${stage}/${sample.id}
    samtools flagstat $bam > "${stage}/${sample.id}/${sample.id}.flagstats.txt"
    head -n 1 "${stage}/${sample.id}/${sample.id}.flagstats.txt" | awk '{print \$1 + \$3}' > "${stage}/${sample.id}/${sample.id}.libsize.txt"
    grep -m 1 "paired in sequencing" "${stage}/${sample.id}/${sample.id}.flagstats.txt" | awk '{npaired = \$1 + \$3; if (npaired==0) {print "unpaired"} else {print "paired"};}' > "${stage}/${sample.id}/${sample.id}.is_paired.txt"
    """
}
