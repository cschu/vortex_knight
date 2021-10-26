process mapseq {
    input:
    tuple val(sample), path(seqs)

    output:
    path("${sample.id}/${sample.id}_R*bac_ssu.mseq"), emit: bac_ssu

    script:
    if (sample.is_paired) {
        """
        mkdir -p ${sample.id}
        ${params.mapseq_bin} ${sample.id}_R1.fastq.gz_bac_ssu.fasta > ${sample.id}/${sample.id}_R1_bac_ssu.mseq
        ${params.mapseq_bin} ${sample.id}_R2.fastq.gz_bac_ssu.fasta > ${sample.id}/${sample.id}_R2_bac_ssu.mseq
        """
    } else {
        """
        mkdir -p ${sample.id}
        ${params.mapseq_bin} ${sample.id}_R1.fastq.gz_bac_ssu.fasta > ${sample.id}/${sample.id}_R1_bac_ssu.mseq
        """
    }
}


process collate_mapseq_tables {
    publishDir params.output_dir, mode: params.publish_mode

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
