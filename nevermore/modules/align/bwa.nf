process bwa_mem_align {
    label 'align'

    input:
    tuple val(sample), path(reads)
    path(reference)
    val(do_name_sort)

    output:
    tuple val(sample), path("${sample.id}.bam"), emit: bam

    script:
    def maxmem = task.memory.toGiga()
    def align_cpus = 4 // figure out the groovy division garbage later (task.cpus >= 8) ?
    def sort_cpus = 4
    def reads2 = (sample.is_paired) ? "${sample.id}_R2.sorted.fastq.gz" : ""
    def sort_reads2 = (sample.is_paired) ? "sortbyname.sh -Xmx${maxmem}g in=${sample.id}_R2.fastq.gz out=${sample.id}_R2.sorted.fastq.gz" : ""
    def blocksize = "-K 10000000"  // shamelessly taken from NGLess

    def sort_cmd = (do_name_sort) ? "samtools collate -@ ${sort_cpus} -o ${sample.id}.bam -" : "samtools sort -@ ${sort_cpus} -o ${sample.id}.bam -"

    """
    set -e -o pipefail
    sortbyname.sh -Xmx${maxmem}g in=${sample.id}_R1.fastq.gz out=${sample.id}_R1.sorted.fastq.gz
    ${sort_reads2}
    bwa mem -a -t ${align_cpus} ${blocksize} \$(readlink ${reference}) ${sample.id}_R1.sorted.fastq.gz ${reads2} | samtools view -F 4 -buSh - | ${sort_cmd}
    """
    // samtools sort -@ ${sort_cpus} -o ${sample.id}.bam

}
