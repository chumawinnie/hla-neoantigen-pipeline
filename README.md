# 🧬 Neoantigen Prediction Pipeline

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg)](https://sylabs.io/docs/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A Nextflow pipeline for **HLA typing** and **neoantigen prediction** from whole exome sequencing (WES) tumor-normal pairs. Designed for clinical genomics workflows and molecular tumor board (MTB) presentations.

## 📋 Table of Contents

- [Overview](#overview)
- [Pipeline Summary](#pipeline-summary)
- [Quick Start](#quick-start)
- [Installation](#installation)
- [Usage](#usage)
- [Input Requirements](#input-requirements)
- [Output](#output)
- [Parameters](#parameters)
- [Neoantigen Prioritization](#neoantigen-prioritization)
- [Troubleshooting](#troubleshooting)
- [Citations](#citations)
- [License](#license)

## Overview

This pipeline integrates HLA typing with neoantigen prediction to identify potential immunotherapy targets from somatic mutations. It combines best-in-class tools for both MHC Class I and Class II epitope prediction.

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                                                                             │
│   Tumor/Normal BAMs ──┬──► HLA Typing ──┬──► Merge HLA ──► pVACseq ──►     │
│                       │    ├─ OptiType   │                    │             │
│   Somatic VCF ────────┤    └─ HLA-HD     │                    ▼             │
│                       │                  │              Filtering &         │
│                       └──► VEP Annotation┘              Ranking             │
│                                                              │              │
│                                                              ▼              │
│                                                      Clinical Report        │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
```

## Pipeline Summary

| Step | Tool | Description |
|------|------|-------------|
| 1. HLA Read Extraction | samtools | Extract reads from HLA region (chr6) |
| 2. HLA Class I Typing | OptiType | 4-digit typing for HLA-A, -B, -C |
| 3. HLA Class II Typing | HLA-HD | Typing for DRB1, DQB1, DPB1 |
| 4. Variant Annotation | VEP | Transcript and protein consequences |
| 5. Neoantigen Prediction | pVACseq | Binding prediction with NetMHCpan/IIpan |
| 6. Filtering & Ranking | Custom | Priority scoring algorithm |
| 7. Report Generation | Custom | MTB-ready HTML reports |

## Quick Start

```bash
# Clone the repository
git clone https://github.com/chumawinnie/neoantigen-pipeline.git
cd neoantigen-pipeline

# Run with test data
nextflow run main.nf \
    --input samplesheet.csv \
    --outdir results \
    --vep_cache /path/to/vep_cache \
    -profile singularity
```

## Installation

### Prerequisites

- [Nextflow](https://www.nextflow.io/) (≥23.04.0)
- [Docker](https://www.docker.com/) or [Singularity](https://sylabs.io/singularity/)
- NetMHCpan 4.1 and NetMHCIIpan 4.0 (academic license from DTU)

### 1. Install Nextflow

```bash
curl -s https://get.nextflow.io | bash
mv nextflow ~/bin/  # or anywhere in your PATH
```

### 2. Obtain NetMHCpan License

NetMHCpan and NetMHCIIpan require an academic license:

1. Visit [DTU Health Tech](https://services.healthtech.dtu.dk/software.php)
2. Register and request licenses for:
   - NetMHCpan 4.1
   - NetMHCIIpan 4.0
3. Download and install following DTU instructions

### 3. Download Reference Data

```bash
# VEP cache (GRCh38)
cd /path/to/vep_cache
curl -O https://ftp.ensembl.org/pub/release-110/variation/indexed_vep_cache/homo_sapiens_vep_110_GRCh38.tar.gz
tar -xzf homo_sapiens_vep_110_GRCh38.tar.gz

# VEP FASTA
curl -O https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
```

### 4. Clone Pipeline

```bash
git clone https://github.com/chumawinnie/neoantigen-pipeline.git
cd neoantigen-pipeline
```

## Usage

### Basic Usage

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --outdir results \
    --vep_cache /path/to/vep_cache \
    -profile singularity
```

### With SLURM (HPC)

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --outdir results \
    --vep_cache /path/to/vep_cache \
    -profile singularity,augsburg
```

### Resume Failed Run

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --outdir results \
    --vep_cache /path/to/vep_cache \
    -profile singularity \
    -resume
```

## Input Requirements

### Samplesheet Format

Create a CSV file with the following columns:

```csv
sample_id,tumor_id,normal_id,tumor_bam,normal_bam,vcf
PATIENT001,PATIENT001_T,PATIENT001_N,/data/PATIENT001_tumor.bam,/data/PATIENT001_normal.bam,/data/PATIENT001_somatic.vcf.gz
PATIENT002,PATIENT002_T,PATIENT002_N,/data/PATIENT002_tumor.bam,/data/PATIENT002_normal.bam,/data/PATIENT002_somatic.vcf.gz
```

| Column | Required | Description |
|--------|----------|-------------|
| `sample_id` | ✓ | Unique sample identifier |
| `tumor_bam` | ✓ | Path to tumor BAM file |
| `normal_bam` | ✓ | Path to normal/germline BAM file |
| `vcf` | ✓ | Path to somatic variant VCF |
| `tumor_id` | | Tumor sample name (default: `{sample_id}_tumor`) |
| `normal_id` | | Normal sample name (default: `{sample_id}_normal`) |

### Input File Requirements

- **BAM files**: Must be sorted and indexed (`.bai` file present)
- **VCF files**: Should contain somatic variants; can be gzipped with tabix index
- **Reference genome**: GRCh38 (hg38) is the default

## Output

```
results/
├── PATIENT001/
│   ├── hla_reads/
│   │   └── PATIENT001_hla_R{1,2}.fastq.gz      # Extracted HLA reads
│   ├── hla_typing/
│   │   ├── class_i/
│   │   │   ├── PATIENT001_optitype_result.tsv  # OptiType results
│   │   │   └── PATIENT001_optitype_coverage.pdf
│   │   ├── class_ii/
│   │   │   └── PATIENT001_hlahd_result.txt     # HLA-HD results
│   │   ├── PATIENT001_hla_alleles.txt          # Combined alleles (pVACseq format)
│   │   └── PATIENT001_hla_summary.json         # Structured HLA summary
│   ├── annotation/
│   │   ├── PATIENT001.vep.vcf.gz               # VEP-annotated VCF
│   │   └── PATIENT001_vep_summary.html         # VEP statistics
│   ├── neoantigens/
│   │   ├── PATIENT001_pvacseq/                 # Full pVACseq output
│   │   │   ├── MHC_Class_I/
│   │   │   └── MHC_Class_II/
│   │   ├── PATIENT001_neoantigens_filtered.tsv # Filtered candidates
│   │   ├── PATIENT001_neoantigens_ranked.tsv   # Priority-ranked candidates
│   │   └── PATIENT001_neoantigen_summary.json  # Summary statistics
│   └── report/
│       └── PATIENT001_neoantigen_report.html   # Clinical report
├── multiqc/
│   └── multiqc_report.html                     # Aggregated QC report
└── pipeline_info/
    ├── execution_timeline.html
    ├── execution_report.html
    └── pipeline_dag.svg
```

### Key Output Files

| File | Description |
|------|-------------|
| `*_hla_alleles.txt` | HLA alleles in pVACseq format (comma-separated) |
| `*_neoantigens_ranked.tsv` | Priority-ranked neoantigen candidates |
| `*_neoantigen_report.html` | Clinical summary for MTB presentation |

## Parameters

### Required Parameters

| Parameter | Description |
|-----------|-------------|
| `--input` | Path to samplesheet CSV |
| `--vep_cache` | Path to VEP cache directory |

### Optional Parameters

#### Output Options
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--outdir` | `./results` | Output directory |
| `--publish_dir_mode` | `copy` | Method to publish results (`copy`, `symlink`, `link`) |

#### Reference Genome
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--genome` | `GRCh38` | Reference genome build |
| `--vep_species` | `homo_sapiens` | VEP species |
| `--vep_assembly` | `GRCh38` | VEP assembly version |
| `--vep_cache_version` | `110` | VEP cache version |

#### HLA Typing
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--hla_reference` | `null` | Path to HLA-HD reference directory |
| `--optitype_data` | `null` | Path to OptiType data directory |

#### Neoantigen Prediction
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--epitope_lengths_class1` | `8,9,10,11` | Peptide lengths for Class I |
| `--epitope_lengths_class2` | `15` | Peptide lengths for Class II |
| `--binding_threshold` | `500` | IC50 binding threshold (nM) |
| `--percentile_threshold` | `2` | Percentile rank threshold |
| `--top_score_metric` | `lowest` | Score selection method |
| `--net_chop_method` | `cterm` | NetChop cleavage method |
| `--netmhc_stab` | `false` | Include stability predictions |

#### Filtering
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--min_vaf` | `0.05` | Minimum variant allele frequency |
| `--min_depth` | `10` | Minimum read depth |
| `--min_alt_reads` | `4` | Minimum alternative allele reads |

#### Resources
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--max_cpus` | `16` | Maximum CPUs per process |
| `--max_memory` | `128.GB` | Maximum memory per process |
| `--max_time` | `72.h` | Maximum time per process |

## Neoantigen Prioritization

Candidates are scored using a multi-factor priority algorithm:

### Scoring Criteria

| Factor | Score | Criteria |
|--------|-------|----------|
| **Binding Affinity** | | |
| | +40 | IC50 < 50 nM (strong binder) |
| | +30 | IC50 < 150 nM |
| | +20 | IC50 < 500 nM |
| | +10 | IC50 < 1000 nM |
| **Percentile Rank** | | |
| | +30 | %Rank < 0.5 |
| | +20 | %Rank < 1.0 |
| | +10 | %Rank < 2.0 |
| **Differential Agretopicity** | | |
| | +20 | WT/MT IC50 ratio > 10 |
| | +15 | WT/MT IC50 ratio > 5 |
| | +10 | WT/MT IC50 ratio > 2 |
| **Clonality (VAF)** | | |
| | +10 | VAF > 30% |
| | +5 | VAF > 10% |

### Interpretation

- **Score ≥ 70**: High-priority candidate
- **Score 40-69**: Medium-priority candidate  
- **Score < 40**: Lower-priority candidate

## Troubleshooting

### Common Issues

#### OptiType Memory Error

```bash
# Increase memory allocation
nextflow run main.nf ... --max_memory 64.GB
```

#### NetMHCpan Not Found

Ensure NetMHCpan is properly installed and in PATH:

```bash
# Check installation
which netMHCpan

# Or specify path in config
params.netmhcpan_path = '/path/to/netMHCpan-4.1'
```

#### VEP Cache Issues

```bash
# Verify cache exists and version matches
ls /path/to/vep_cache/homo_sapiens/110_GRCh38/

# Specify version explicitly
nextflow run main.nf ... --vep_cache_version 110
```

#### Singularity Bind Mounts

```bash
# Set bind paths for your data
export SINGULARITY_BIND="/data:/data,/scratch:/scratch"
```

### Getting Help

```bash
# Show all parameters
nextflow run main.nf --help

# Show pipeline version
nextflow run main.nf --version
```

### Logs

Check these locations for debugging:

- `.nextflow.log` - Nextflow execution log
- `results/pipeline_info/` - Execution reports
- `work/` - Process working directories

## Citations

If you use this pipeline, please cite:

### Tools

- **OptiType**: Szolek A, et al. OptiType: precision HLA typing from next-generation sequencing data. *Bioinformatics*. 2014;30(23):3310-3316. doi:[10.1093/bioinformatics/btu548](https://doi.org/10.1093/bioinformatics/btu548)

- **HLA-HD**: Kawaguchi S, et al. HLA-HD: An accurate HLA typing algorithm for next-generation sequencing data. *Human Mutation*. 2017;38(7):788-797. doi:[10.1002/humu.23230](https://doi.org/10.1002/humu.23230)

- **pVACtools**: Hundal J, et al. pVACtools: A computational toolkit to identify and visualize cancer neoantigens. *Cancer Immunology Research*. 2020;8(3):409-420. doi:[10.1158/2326-6066.CIR-19-0401](https://doi.org/10.1158/2326-6066.CIR-19-0401)

- **NetMHCpan**: Reynisson B, et al. NetMHCpan-4.1 and NetMHCIIpan-4.0: improved predictions of MHC antigen presentation by concurrent motif deconvolution and integration of MS MHC eluted ligand data. *Nucleic Acids Research*. 2020;48(W1):W449-W454. doi:[10.1093/nar/gkaa379](https://doi.org/10.1093/nar/gkaa379)

- **Ensembl VEP**: McLaren W, et al. The Ensembl Variant Effect Predictor. *Genome Biology*. 2016;17:122. doi:[10.1186/s13059-016-0974-4](https://doi.org/10.1186/s13059-016-0974-4)

### Workflow

- **Nextflow**: Di Tommaso P, et al. Nextflow enables reproducible computational workflows. *Nature Biotechnology*. 2017;35:316-319. doi:[10.1038/nbt.3820](https://doi.org/10.1038/nbt.3820)

## Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit changes (`git commit -m 'Add amazing feature'`)
4. Push to branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Authors

**Bioinformatics Core Facility**  
University of Augsburg

---

<p align="center">
  <i>Built for clinical genomics and molecular tumor board workflows</i>
</p>
