process merge_and_sort {
    label 'samtools'
    publishDir params.output_dir, mode: params.publish_mode

    input:
    tuple val(sample), path(bamfiles)
    val(do_name_sort)

    output:
    tuple val(sample), path("bam/${sample}.bam"), emit: bam
    tuple val(sample), path("stats/bam/${sample}.flagstats.txt"), emit: flagstats

    script:
    def sort_order = (do_name_sort) ? "-n" : ""
    // need a better detection for this
    if (bamfiles instanceof Collection && bamfiles.size() >= 2) {
        """
        mkdir -p bam/ stats/bam
        samtools merge -@ $task.cpus ${sort_order} bam/${sample}.bam ${bamfiles}
        samtools flagstats bam/${sample}.bam > stats/bam/${sample}.flagstats.txt
        """
    } else {
        // i don't like this solution
        """
        mkdir -p bam/ stats/bam
        ln -s ../${bamfiles[0]} bam/${sample}.bam
        samtools flagstats bam/${sample}.bam > stats/bam/${sample}.flagstats.txt
        """
    }
}
