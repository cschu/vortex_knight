process bwa_mem_align {
    cpus 4
    memory { 20.GB * task.attempt }
    time { 3.d * task.attempt }
    container "registry.git.embl.de/schudoma/align-docker"
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
    def reads2 = (sample.is_paired) ? "${sample.id}_R2.fastq.gz" : ""
    def sort_reads2 = (sample.is_paired) ? "sortbyname.sh -Xmx${maxmem}g in=${sample.id}_R2.fastq.gz out=${sample.id}_R2.sorted.fastq.gz" : ""
    def blocksize = "-K 10000000"  // shamelessly taken from NGLess

    def sort_cmd = "samtools collate -@ ${sort_cpus} -o ${sample.id}.bam - tmp/collated_bam"

    // sortbyname.sh -Xmx${maxmem}g in=${sample.id}_R1.fastq.gz out=${sample.id}_R1.sorted.fastq.gz
    // bwa mem -a -t ${align_cpus} ${blocksize} \$(readlink ${reference}) ${sample.id}_R1.sorted.fastq.gz ${reads2} | samtools view -F 4 -buSh - | ${sort_cmd}
    // ${sort_reads2}
    """
    set -e -o pipefail
    mkdir -p tmp/
    bwa mem -a -t ${align_cpus} ${blocksize} \$(readlink ${reference}) ${reads} | samtools view -buSh - | ${sort_cmd}
    rm -rvf tmp/ *.sorted.fastq.gz
    """
}
