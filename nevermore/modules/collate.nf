process collate_stats {
    // container "ghcr.io/astral-sh/uv:python3.14-trixie-slim"
    container "registry.git.embl.org/schudoma/portraits_metatraits:latest"
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

