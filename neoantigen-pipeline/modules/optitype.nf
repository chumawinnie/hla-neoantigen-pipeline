/*
========================================================================================
    OPTITYPE - HLA CLASS I TYPING
========================================================================================
    4-digit HLA Class I typing (HLA-A, HLA-B, HLA-C) from WES/WGS data
*/

process OPTITYPE {
    tag "$meta.id"
    label 'process_high'
    
    publishDir "${params.outdir}/${meta.id}/hla_typing/class_i", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}_optitype_result.tsv"),    emit: hla_types
    tuple val(meta), path("${meta.id}_optitype_coverage.pdf"),  emit: coverage_plot, optional: true
    path "*.log",                                                emit: logs

    script:
    def r1 = reads[0]
    def r2 = reads[1]
    """
    # Create output directory
    mkdir -p optitype_out
    
    # Run OptiType
    OptiTypePipeline.py \\
        -i ${r1} ${r2} \\
        --dna \\
        --outdir optitype_out \\
        --prefix ${meta.id} \\
        --verbose
    
    # Rename outputs
    mv optitype_out/${meta.id}_result.tsv ${meta.id}_optitype_result.tsv
    if [ -f optitype_out/${meta.id}_coverage_plot.pdf ]; then
        mv optitype_out/${meta.id}_coverage_plot.pdf ${meta.id}_optitype_coverage.pdf
    fi
    
    # Create log
    echo "OptiType HLA Class I typing for ${meta.id}" > ${meta.id}_optitype.log
    cat ${meta.id}_optitype_result.tsv >> ${meta.id}_optitype.log
    
    # Cleanup
    rm -rf optitype_out
    """
}
