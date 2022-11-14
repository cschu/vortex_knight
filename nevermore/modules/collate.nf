process collate_stats {
    // publishDir params.output_dir, mode: params.publish_mode

    input:
    path(stats_files)

    output:
    path("reports/read_count_table.txt")

    script:
    """
    mkdir -p reports/
    collate_stats.py . > reports/read_count_table.txt
    """
}

