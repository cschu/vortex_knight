process collate_stats {
    container "ghcr.io/astral-sh/uv:python3.14-trixie-slim"
    label "default"

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

