/*
========================================================================================
    VEP ANNOTATE
========================================================================================
    Annotate somatic VCF with Ensembl VEP for protein consequence prediction
*/

process VEP_ANNOTATE {
    tag "$meta.id"
    label 'process_medium'
    
    publishDir "${params.outdir}/${meta.id}/annotation", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("${meta.id}.vep.vcf.gz"),     emit: annotated_vcf
    tuple val(meta), path("${meta.id}.vep.vcf.gz.tbi"), emit: annotated_vcf_tbi
    path "${meta.id}_vep_summary.html",                  emit: summary
    path "*.log",                                        emit: logs

    script:
    def vep_cache = params.vep_cache
    def species   = params.vep_species
    def assembly  = params.vep_assembly
    def cache_ver = params.vep_cache_version
    """
    # Run VEP with annotations required for pVACseq
    vep \\
        --input_file ${vcf} \\
        --output_file ${meta.id}.vep.vcf \\
        --format vcf \\
        --vcf \\
        --symbol \\
        --terms SO \\
        --tsl \\
        --hgvs \\
        --fasta ${vep_cache}/homo_sapiens/${cache_ver}_${assembly}/Homo_sapiens.${assembly}.dna.toplevel.fa.gz \\
        --offline \\
        --cache \\
        --dir_cache ${vep_cache} \\
        --species ${species} \\
        --assembly ${assembly} \\
        --cache_version ${cache_ver} \\
        --pick \\
        --plugin Frameshift \\
        --plugin Wildtype \\
        --fork ${task.cpus} \\
        --stats_file ${meta.id}_vep_summary.html \\
        --warning_file ${meta.id}_vep_warnings.log
    
    # Compress and index
    bgzip -c ${meta.id}.vep.vcf > ${meta.id}.vep.vcf.gz
    tabix -p vcf ${meta.id}.vep.vcf.gz
    
    # Create log
    echo "VEP annotation completed for ${meta.id}" > ${meta.id}_vep.log
    grep -c "^#" ${meta.id}.vep.vcf || true >> ${meta.id}_vep.log
    echo "Total variants annotated:" >> ${meta.id}_vep.log
    zgrep -v "^#" ${meta.id}.vep.vcf.gz | wc -l >> ${meta.id}_vep.log
    
    # Cleanup
    rm -f ${meta.id}.vep.vcf
    """
}
