process merge_and_sort {
    label 'samtools'

    input:
    tuple val(sample), path(bamfiles)
    val(do_name_sort)

    output:
    tuple val(sample), path("bam/${sample}.bam"), emit: bam
    tuple val(sample), path("stats/bam/${sample}.flagstats.txt"), emit: flagstats

    script:
    def sort_order = (do_name_sort) ? "-n" : ""
    def merge_cmd = ""

    // need a better detection for this
    if (bamfiles instanceof Collection && bamfiles.size() >= 2) {
        merge_cmd = "samtools merge -@ $task.cpus ${sort_order} bam/${sample}.bam ${bamfiles}"
    } else {
        merge_cmd = "ln -s ../${bamfiles[0]} bam/${sample}.bam"
    }

    """
    mkdir -p bam/ stats/bam/
    ${merge_cmd}
    samtools flagstats bam/${sample}.bam > stats/bam/${sample}.flagstats.txt
    """
}


process db_filter {
    label 'samtools'

    input:
    tuple val(sample), path(bam)
    path(db_bedfile)

    output:
    tuple val(sample), path("filtered_bam/${sample}.bam"), emit: bam
    tuple val(sample), path("stats/filtered_bam/${sample}.flagstats.txt"), emit: flagstats

    script:
    """
    mkdir -p filtered_bam/ stats/filtered_bam/
    bedtools intersect -u -ubam -a ${bam} -b ${db_bedfile} > filtered_bam/${sample}.bam
    samtools flagstats filtered_bam/${sample}.bam > stats/filtered_bam/${sample}.flagstats.txt
    """
}

process db2bed3 {
    input:
    path(db)

    output:
    path("db.bed3"), emit: db

    script:
    """
    sqlite3 ${db} 'select seqid,start,end from annotatedsequence;' '.exit' | awk -F \\| -v OFS='\\t' '{print \$1,\$2-1,\$3}' | sort -k1,1 -k2,2g -k3,3g > db.bed3
    """
}
