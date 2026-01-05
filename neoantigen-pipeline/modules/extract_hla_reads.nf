/*
========================================================================================
    EXTRACT HLA READS
========================================================================================
    Extract reads mapping to HLA region from BAM for HLA typing
*/

process EXTRACT_HLA_READS {
    tag "$meta.id"
    label 'process_medium'
    
    publishDir "${params.outdir}/${meta.id}/hla_reads", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}_hla_R{1,2}.fastq.gz"), emit: fastq
    path "*.log", emit: logs

    script:
    // HLA region on chr6 (GRCh38 coordinates)
    def hla_region = "chr6:28510120-33480577"
    """
    # Extract reads from HLA region
    samtools view -b -h ${bam} ${hla_region} > ${meta.id}_hla.bam
    
    # Also extract unmapped reads (may contain HLA reads not mapped to reference)
    samtools view -b -f 4 ${bam} > ${meta.id}_unmapped.bam
    
    # Merge
    samtools merge -f ${meta.id}_hla_combined.bam ${meta.id}_hla.bam ${meta.id}_unmapped.bam
    
    # Sort by name for paired-end extraction
    samtools sort -n -@ ${task.cpus} ${meta.id}_hla_combined.bam -o ${meta.id}_hla_sorted.bam
    
    # Extract to FASTQ
    samtools fastq \\
        -@ ${task.cpus} \\
        -1 ${meta.id}_hla_R1.fastq.gz \\
        -2 ${meta.id}_hla_R2.fastq.gz \\
        -0 /dev/null \\
        -s /dev/null \\
        -n \\
        ${meta.id}_hla_sorted.bam
    
    # Log stats
    echo "HLA read extraction for ${meta.id}" > ${meta.id}_extract_hla.log
    echo "HLA region reads:" >> ${meta.id}_extract_hla.log
    samtools view -c ${meta.id}_hla.bam >> ${meta.id}_extract_hla.log
    echo "Unmapped reads:" >> ${meta.id}_extract_hla.log
    samtools view -c ${meta.id}_unmapped.bam >> ${meta.id}_extract_hla.log
    
    # Cleanup
    rm -f ${meta.id}_hla.bam ${meta.id}_unmapped.bam ${meta.id}_hla_combined.bam ${meta.id}_hla_sorted.bam
    """
}
