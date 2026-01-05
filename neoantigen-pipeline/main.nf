#!/usr/bin/env nextflow

/*
========================================================================================
    NEOANTIGEN PREDICTION PIPELINE
========================================================================================
    A Nextflow pipeline for HLA typing and neoantigen prediction from WES tumor-normal pairs
    
    Components:
    - HLA Class I typing (OptiType)
    - HLA Class II typing (HLA-HD)
    - Neoantigen prediction (pVACtools + NetMHCpan/NetMHCIIpan)
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    PARAMETER DEFINITIONS
========================================================================================
*/

params.input            = null                      // Samplesheet CSV
params.outdir           = './results'               // Output directory
params.genome           = 'GRCh38'                  // Reference genome
params.vep_cache        = null                      // Path to VEP cache
params.vep_species      = 'homo_sapiens'            // VEP species
params.vep_assembly     = 'GRCh38'                  // VEP assembly

// HLA typing parameters
params.hla_reference    = null                      // HLA-HD reference path
params.optitype_data    = null                      // OptiType data directory

// pVACtools parameters
params.epitope_lengths_class1  = '8,9,10,11'        // Class I peptide lengths
params.epitope_lengths_class2  = '15'               // Class II peptide lengths
params.binding_threshold       = 500                // IC50 binding threshold (nM)
params.percentile_threshold    = 2                  // Percentile rank threshold
params.top_score_metric        = 'lowest'           // 'lowest' or 'median'
params.net_chop_method         = 'cterm'            // NetChop cleavage method
params.netmhc_stab             = false              // Include stability predictions

// Filtering parameters
params.min_vaf          = 0.05                      // Minimum variant allele frequency
params.min_depth        = 10                        // Minimum read depth
params.min_alt_reads    = 4                         // Minimum alt allele reads

// Resource defaults
params.max_cpus         = 16
params.max_memory       = '64.GB'
params.max_time         = '48.h'

/*
========================================================================================
    INPUT VALIDATION
========================================================================================
*/

if (!params.input) {
    error "Please provide a samplesheet with --input"
}

if (!params.vep_cache) {
    error "Please provide VEP cache path with --vep_cache"
}

/*
========================================================================================
    IMPORT MODULES
========================================================================================
*/

include { EXTRACT_HLA_READS       } from './modules/extract_hla_reads'
include { OPTITYPE                } from './modules/optitype'
include { HLA_HD                  } from './modules/hla_hd'
include { MERGE_HLA_TYPES         } from './modules/merge_hla_types'
include { VEP_ANNOTATE            } from './modules/vep_annotate'
include { PVACSEQ                 } from './modules/pvacseq'
include { FILTER_NEOANTIGENS      } from './modules/filter_neoantigens'
include { GENERATE_REPORT         } from './modules/generate_report'
include { MULTIQC                 } from './modules/multiqc'

/*
========================================================================================
    HELPER FUNCTIONS
========================================================================================
*/

def validateInputSamplesheet(LinkedHashMap row) {
    // Validate required columns exist
    def required = ['sample_id', 'tumor_bam', 'normal_bam', 'vcf']
    def missing = required.findAll { !row.containsKey(it) || !row[it] }
    if (missing) {
        error "Missing required columns in samplesheet: ${missing.join(', ')}"
    }
    return row
}

/*
========================================================================================
    MAIN WORKFLOW
========================================================================================
*/

workflow {

    // Parse samplesheet
    Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        .map { row -> validateInputSamplesheet(row) }
        .map { row ->
            def meta = [
                id: row.sample_id,
                tumor_id: row.tumor_id ?: "${row.sample_id}_tumor",
                normal_id: row.normal_id ?: "${row.sample_id}_normal"
            ]
            return [
                meta,
                file(row.tumor_bam, checkIfExists: true),
                file(row.tumor_bam + '.bai', checkIfExists: true),
                file(row.normal_bam, checkIfExists: true),
                file(row.normal_bam + '.bai', checkIfExists: true),
                file(row.vcf, checkIfExists: true)
            ]
        }
        .set { ch_input }

    // Split channels for different processes
    ch_input
        .map { meta, tumor_bam, tumor_bai, normal_bam, normal_bai, vcf ->
            [meta, normal_bam, normal_bai]
        }
        .set { ch_normal_bam }

    ch_input
        .map { meta, tumor_bam, tumor_bai, normal_bam, normal_bai, vcf ->
            [meta, vcf]
        }
        .set { ch_vcf }

    /*
    ================================================================================
        HLA TYPING
    ================================================================================
    */

    // Extract HLA region reads from normal BAM (germline HLA typing)
    EXTRACT_HLA_READS(ch_normal_bam)

    // Class I HLA typing with OptiType
    OPTITYPE(EXTRACT_HLA_READS.out.fastq)

    // Class II HLA typing with HLA-HD
    HLA_HD(EXTRACT_HLA_READS.out.fastq)

    // Merge HLA types from both tools
    OPTITYPE.out.hla_types
        .join(HLA_HD.out.hla_types)
        .set { ch_hla_both }

    MERGE_HLA_TYPES(ch_hla_both)

    /*
    ================================================================================
        VARIANT ANNOTATION
    ================================================================================
    */

    // Annotate VCF with VEP (transcript, protein consequences)
    VEP_ANNOTATE(ch_vcf)

    /*
    ================================================================================
        NEOANTIGEN PREDICTION
    ================================================================================
    */

    // Combine HLA types with annotated VCF
    MERGE_HLA_TYPES.out.merged_hla
        .join(VEP_ANNOTATE.out.annotated_vcf)
        .set { ch_pvacseq_input }

    // Run pVACseq for neoantigen prediction
    PVACSEQ(ch_pvacseq_input)

    // Filter and rank neoantigens
    FILTER_NEOANTIGENS(PVACSEQ.out.predictions)

    /*
    ================================================================================
        REPORTING
    ================================================================================
    */

    // Collect all results for report generation
    FILTER_NEOANTIGENS.out.filtered
        .join(MERGE_HLA_TYPES.out.merged_hla)
        .set { ch_report_input }

    GENERATE_REPORT(ch_report_input)

    // MultiQC summary
    Channel.empty()
        .mix(OPTITYPE.out.logs)
        .mix(HLA_HD.out.logs)
        .mix(VEP_ANNOTATE.out.logs)
        .mix(PVACSEQ.out.logs)
        .collect()
        .set { ch_multiqc }

    MULTIQC(ch_multiqc)
}

/*
========================================================================================
    WORKFLOW COMPLETION
========================================================================================
*/

workflow.onComplete {
    log.info ""
    log.info "Pipeline completed at: ${workflow.complete}"
    log.info "Execution status: ${workflow.success ? 'SUCCESS' : 'FAILED'}"
    log.info "Duration: ${workflow.duration}"
    log.info "Results: ${params.outdir}"
    log.info ""
}

workflow.onError {
    log.error "Pipeline execution failed"
    log.error "Error message: ${workflow.errorMessage}"
}
