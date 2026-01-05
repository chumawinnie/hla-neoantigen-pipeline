/*
========================================================================================
    MULTIQC
========================================================================================
    Aggregate QC reports across all samples
*/

process MULTIQC {
    label 'process_low'
    
    publishDir "${params.outdir}/multiqc", mode: params.publish_dir_mode

    input:
    path logs

    output:
    path "multiqc_report.html",     emit: report
    path "multiqc_data/",           emit: data

    script:
    """
    multiqc . \\
        --force \\
        --title "Neoantigen Pipeline QC Report" \\
        --comment "Quality control summary for neoantigen prediction pipeline"
    """
}
