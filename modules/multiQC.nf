process MULTIQC {
    
    publishDir "${params.outdir}/quality_reports", mode:'copy'

    input:
    path '*'

    output:
    path 'multiqc_report.html'

    script:
    """
    multiqc .
    """

}