# Usage Guide

This guide covers how to run the neoantigen prediction pipeline with your data.

## Table of Contents

- [Preparing Input Data](#preparing-input-data)
- [Creating a Samplesheet](#creating-a-samplesheet)
- [Running the Pipeline](#running-the-pipeline)
- [Configuration Options](#configuration-options)
- [Running on HPC](#running-on-hpc)
- [Monitoring Execution](#monitoring-execution)

## Preparing Input Data

### Required Input Files

For each sample, you need:

1. **Tumor BAM** - Aligned tumor WES reads
2. **Normal BAM** - Aligned normal/germline WES reads  
3. **Somatic VCF** - Called somatic variants

### BAM File Requirements

```bash
# BAM files must be:
# 1. Sorted by coordinate
samtools sort -o tumor_sorted.bam tumor.bam

# 2. Indexed
samtools index tumor_sorted.bam  # Creates .bai file

# 3. Aligned to GRCh38 (hg38)
```

### VCF Requirements

The somatic VCF should:

- Contain only somatic variants (tumor-specific)
- Be generated from tumor-normal variant calling
- Include standard VCF fields (CHROM, POS, REF, ALT)
- Optionally include VAF/AF in FORMAT or INFO

Supported callers:
- Mutect2 (recommended)
- Strelka2
- VarDict
- Any standard VCF format

```bash
# If VCF is not compressed
bgzip somatic_variants.vcf
tabix -p vcf somatic_variants.vcf.gz
```

## Creating a Samplesheet

Create a CSV file listing your samples:

### Basic Samplesheet

```csv
sample_id,tumor_id,normal_id,tumor_bam,normal_bam,vcf
PATIENT001,PATIENT001_T,PATIENT001_N,/data/bams/PATIENT001_tumor.bam,/data/bams/PATIENT001_normal.bam,/data/vcfs/PATIENT001_somatic.vcf.gz
PATIENT002,PATIENT002_T,PATIENT002_N,/data/bams/PATIENT002_tumor.bam,/data/bams/PATIENT002_normal.bam,/data/vcfs/PATIENT002_somatic.vcf.gz
```

### Column Descriptions

| Column | Required | Description |
|--------|----------|-------------|
| `sample_id` | Yes | Unique identifier for the sample |
| `tumor_id` | No | Tumor sample name in BAM (defaults to `{sample_id}_tumor`) |
| `normal_id` | No | Normal sample name in BAM (defaults to `{sample_id}_normal`) |
| `tumor_bam` | Yes | Full path to tumor BAM file |
| `normal_bam` | Yes | Full path to normal BAM file |
| `vcf` | Yes | Full path to somatic VCF file |

### Integration with Existing Pipeline

If you're using output from your clinical genomics pipeline:

```csv
sample_id,tumor_id,normal_id,tumor_bam,normal_bam,vcf
SAMPLE001,SAMPLE001_tumor,SAMPLE001_normal,/results/preprocessing/recal/SAMPLE001_tumor.recal.bam,/results/preprocessing/recal/SAMPLE001_normal.recal.bam,/results/variant_calling/mutect2/SAMPLE001/SAMPLE001.filtered.vcf.gz
```

## Running the Pipeline

### Basic Execution

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --outdir results \
    --vep_cache /path/to/vep_cache \
    -profile singularity
```

### With All Options

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --outdir results \
    --genome GRCh38 \
    --vep_cache /path/to/vep_cache \
    --vep_cache_version 110 \
    --hla_reference /path/to/hlahd/data \
    --epitope_lengths_class1 '8,9,10,11' \
    --epitope_lengths_class2 '15' \
    --binding_threshold 500 \
    --percentile_threshold 2 \
    --min_vaf 0.05 \
    --max_cpus 16 \
    --max_memory 64.GB \
    -profile singularity
```

### Resume After Failure

```bash
# Use -resume to continue from where it stopped
nextflow run main.nf \
    --input samplesheet.csv \
    --outdir results \
    --vep_cache /path/to/vep_cache \
    -profile singularity \
    -resume
```

### Running Specific Samples

Create a filtered samplesheet with only the samples you want to process:

```bash
# Extract specific samples
head -1 samplesheet.csv > subset.csv
grep "PATIENT001\|PATIENT002" samplesheet.csv >> subset.csv

nextflow run main.nf --input subset.csv ...
```

## Configuration Options

### Using a Config File

Create `my_params.config`:

```groovy
params {
    input       = 'samplesheet.csv'
    outdir      = 'results'
    vep_cache   = '/data/references/vep_cache'
    
    // Prediction settings
    binding_threshold     = 500
    percentile_threshold  = 2
    min_vaf              = 0.05
}
```

Run with:

```bash
nextflow run main.nf -c my_params.config -profile singularity
```

### Environment Variables

```bash
# Set default paths
export NXF_VEP_CACHE=/path/to/vep_cache
export NXF_WORK=/scratch/nextflow_work

nextflow run main.nf \
    --input samplesheet.csv \
    --vep_cache $NXF_VEP_CACHE \
    -w $NXF_WORK \
    -profile singularity
```

## Running on HPC

### SLURM Configuration

Use the `augsburg` profile or create a custom profile:

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --outdir results \
    --vep_cache /path/to/vep_cache \
    -profile singularity,augsburg
```

### Custom SLURM Profile

Add to `nextflow.config`:

```groovy
profiles {
    my_cluster {
        process {
            executor = 'slurm'
            queue = 'batch'
            clusterOptions = '--account=my_account'
        }
        
        executor {
            queueSize = 50
            submitRateLimit = '10 sec'
        }
        
        params {
            max_cpus   = 32
            max_memory = '256.GB'
            max_time   = '168.h'
        }
    }
}
```

### Singularity on HPC

```bash
# Set Singularity bind paths
export SINGULARITY_BIND="/data:/data,/scratch:/scratch,/home:/home"

# Set cache directory
export SINGULARITY_CACHEDIR=/scratch/$USER/singularity_cache
mkdir -p $SINGULARITY_CACHEDIR

nextflow run main.nf \
    --input samplesheet.csv \
    -profile singularity,my_cluster
```

## Monitoring Execution

### Real-time Monitoring

```bash
# Nextflow log
tail -f .nextflow.log

# Process status
nextflow log last
```

### Nextflow Tower

For advanced monitoring, use [Nextflow Tower](https://tower.nf/):

```bash
# Set Tower access token
export TOWER_ACCESS_TOKEN=your_token

nextflow run main.nf \
    --input samplesheet.csv \
    -profile singularity \
    -with-tower
```

### Execution Reports

After completion, check:

- `results/pipeline_info/execution_report.html` - Resource usage
- `results/pipeline_info/execution_timeline.html` - Timeline
- `results/pipeline_info/pipeline_dag.svg` - Workflow DAG
- `results/multiqc/multiqc_report.html` - QC summary

## Common Use Cases

### Running a Single Sample

```bash
# Create single-sample samplesheet
echo "sample_id,tumor_id,normal_id,tumor_bam,normal_bam,vcf" > single.csv
echo "PATIENT001,PATIENT001_T,PATIENT001_N,/path/tumor.bam,/path/normal.bam,/path/somatic.vcf.gz" >> single.csv

nextflow run main.nf --input single.csv -profile singularity
```

### Batch Processing

```bash
# Process all samples in a directory
for vcf in /data/vcfs/*.vcf.gz; do
    sample=$(basename $vcf .somatic.vcf.gz)
    echo "$sample,${sample}_T,${sample}_N,/data/bams/${sample}_tumor.bam,/data/bams/${sample}_normal.bam,$vcf"
done > batch_samplesheet.csv

# Add header
sed -i '1i sample_id,tumor_id,normal_id,tumor_bam,normal_bam,vcf' batch_samplesheet.csv

nextflow run main.nf --input batch_samplesheet.csv -profile singularity
```

### Re-running with Different Parameters

```bash
# Run with stricter filtering
nextflow run main.nf \
    --input samplesheet.csv \
    --outdir results_strict \
    --binding_threshold 150 \
    --percentile_threshold 0.5 \
    --min_vaf 0.1 \
    -profile singularity
```

## Next Steps

- Review [Output Documentation](output.md) to understand results
- Check [Troubleshooting Guide](troubleshooting.md) for common issues
- See [Parameters Reference](parameters.md) for all options
