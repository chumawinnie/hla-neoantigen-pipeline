/*
========================================================================================
    PVACSEQ - NEOANTIGEN PREDICTION
========================================================================================
    Run pVACseq for neoantigen prediction using NetMHCpan and NetMHCIIpan
*/

process PVACSEQ {
    tag "$meta.id"
    label 'process_high'
    
    publishDir "${params.outdir}/${meta.id}/neoantigens", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(hla_alleles), path(vcf), path(vcf_tbi)

    output:
    tuple val(meta), path("${meta.id}_pvacseq/"),                       emit: predictions
    tuple val(meta), path("${meta.id}_pvacseq/MHC_Class_I/*.tsv"),      emit: class_i_results, optional: true
    tuple val(meta), path("${meta.id}_pvacseq/MHC_Class_II/*.tsv"),     emit: class_ii_results, optional: true
    path "*.log",                                                        emit: logs

    script:
    def hla_list = hla_alleles.text.trim()
    def epitope_class1 = params.epitope_lengths_class1
    def epitope_class2 = params.epitope_lengths_class2
    def binding_thresh = params.binding_threshold
    def percentile_thresh = params.percentile_threshold
    def top_score = params.top_score_metric
    def net_chop = params.net_chop_method
    def netmhc_stab = params.netmhc_stab ? '--netmhc-stab' : ''
    
    // Separate Class I and Class II alleles
    """
    #!/bin/bash
    set -euo pipefail
    
    # Read HLA alleles
    HLA_ALLELES="${hla_list}"
    
    # Separate Class I (A, B, C) and Class II (D*) alleles
    CLASS_I_ALLELES=\$(echo "\$HLA_ALLELES" | tr ',' '\\n' | grep -E 'HLA-[ABC]' | tr '\\n' ',' | sed 's/,\$//')
    CLASS_II_ALLELES=\$(echo "\$HLA_ALLELES" | tr ',' '\\n' | grep -E 'HLA-D' | tr '\\n' ',' | sed 's/,\$//')
    
    echo "Class I alleles: \$CLASS_I_ALLELES" > ${meta.id}_pvacseq.log
    echo "Class II alleles: \$CLASS_II_ALLELES" >> ${meta.id}_pvacseq.log
    
    # Create output directory
    mkdir -p ${meta.id}_pvacseq
    
    # Run pVACseq for Class I predictions
    if [ -n "\$CLASS_I_ALLELES" ]; then
        echo "Running Class I prediction..." >> ${meta.id}_pvacseq.log
        
        pvacseq run \\
            ${vcf} \\
            ${meta.id} \\
            "\$CLASS_I_ALLELES" \\
            NetMHCpan \\
            ${meta.id}_pvacseq \\
            --iedb-install-directory /opt/iedb \\
            -e1 ${epitope_class1} \\
            --binding-threshold ${binding_thresh} \\
            --percentile-threshold ${percentile_thresh} \\
            --top-score-metric ${top_score} \\
            --net-chop-method ${net_chop} \\
            ${netmhc_stab} \\
            -t ${task.cpus} \\
            --keep-tmp-files
            
        echo "Class I prediction completed" >> ${meta.id}_pvacseq.log
    fi
    
    # Run pVACseq for Class II predictions
    if [ -n "\$CLASS_II_ALLELES" ]; then
        echo "Running Class II prediction..." >> ${meta.id}_pvacseq.log
        
        pvacseq run \\
            ${vcf} \\
            ${meta.id} \\
            "\$CLASS_II_ALLELES" \\
            NetMHCIIpan \\
            ${meta.id}_pvacseq \\
            --iedb-install-directory /opt/iedb \\
            -e2 ${epitope_class2} \\
            --binding-threshold ${binding_thresh} \\
            --percentile-threshold ${percentile_thresh} \\
            --top-score-metric ${top_score} \\
            -t ${task.cpus} \\
            --keep-tmp-files
            
        echo "Class II prediction completed" >> ${meta.id}_pvacseq.log
    fi
    
    # Count results
    echo "\\n=== Results Summary ===" >> ${meta.id}_pvacseq.log
    
    if [ -d "${meta.id}_pvacseq/MHC_Class_I" ]; then
        CLASS_I_COUNT=\$(find ${meta.id}_pvacseq/MHC_Class_I -name "*.filtered.tsv" -exec wc -l {} + 2>/dev/null | tail -1 | awk '{print \$1-1}' || echo "0")
        echo "Class I neoantigens (filtered): \$CLASS_I_COUNT" >> ${meta.id}_pvacseq.log
    fi
    
    if [ -d "${meta.id}_pvacseq/MHC_Class_II" ]; then
        CLASS_II_COUNT=\$(find ${meta.id}_pvacseq/MHC_Class_II -name "*.filtered.tsv" -exec wc -l {} + 2>/dev/null | tail -1 | awk '{print \$1-1}' || echo "0")
        echo "Class II neoantigens (filtered): \$CLASS_II_COUNT" >> ${meta.id}_pvacseq.log
    fi
    """
}
