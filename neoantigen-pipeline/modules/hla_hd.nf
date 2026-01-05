/*
========================================================================================
    HLA-HD - HLA CLASS II TYPING
========================================================================================
    High-resolution HLA typing including Class II (DRB1, DQB1, DPB1) from WES/WGS
*/

process HLA_HD {
    tag "$meta.id"
    label 'process_high'
    
    publishDir "${params.outdir}/${meta.id}/hla_typing/class_ii", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}_hlahd_result.txt"),   emit: hla_types
    path "*.log",                                            emit: logs

    script:
    def r1 = reads[0]
    def r2 = reads[1]
    def hla_ref = params.hla_reference ?: '/usr/local/share/hlahd/data'
    """
    # Run HLA-HD
    hlahd.sh \\
        -t ${task.cpus} \\
        -m 50 \\
        -f ${hla_ref}/freq_data \\
        ${r1} \\
        ${r2} \\
        ${hla_ref}/gene_split.txt \\
        ${hla_ref}/dictionary \\
        ${meta.id} \\
        .
    
    # HLA-HD outputs to a subdirectory
    mv ${meta.id}/result/${meta.id}_final.result.txt ${meta.id}_hlahd_result.txt
    
    # Create log
    echo "HLA-HD typing results for ${meta.id}" > ${meta.id}_hlahd.log
    cat ${meta.id}_hlahd_result.txt >> ${meta.id}_hlahd.log
    
    # Cleanup
    rm -rf ${meta.id}
    """
}
