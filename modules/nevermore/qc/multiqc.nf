process multiqc {
    publishDir params.output_dir, mode: params.publish_mode

    input:
    path(reports)

    output:
    path("multiqc_report.html")

    script:
    def send_report = (params.email) ? "echo . | mailx -s 'multiqc_report' -a multiqc_report.html ${params.email}" : ""
    """
    multiqc -c ${projectDir}/config/multiqc.config .
    ${send_report}
    """
}
