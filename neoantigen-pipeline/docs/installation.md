# Installation Guide

This guide provides detailed instructions for setting up the neoantigen prediction pipeline.

## Table of Contents

- [System Requirements](#system-requirements)
- [Installing Nextflow](#installing-nextflow)
- [Container Setup](#container-setup)
- [Reference Data](#reference-data)
- [NetMHCpan Setup](#netmhcpan-setup)
- [HLA-HD Setup](#hla-hd-setup)
- [Testing Installation](#testing-installation)

## System Requirements

### Minimum Requirements

| Resource | Minimum | Recommended |
|----------|---------|-------------|
| CPU | 8 cores | 16+ cores |
| RAM | 32 GB | 64+ GB |
| Storage | 100 GB | 500+ GB |
| OS | Linux (64-bit) | CentOS 7+, Ubuntu 18.04+ |

### Software Dependencies

- Nextflow ≥23.04.0
- Java 11 or later
- Docker 20.10+ OR Singularity 3.8+

## Installing Nextflow

### Option 1: Quick Install

```bash
# Download and install
curl -s https://get.nextflow.io | bash

# Move to PATH
sudo mv nextflow /usr/local/bin/

# Verify installation
nextflow -version
```

### Option 2: Conda Install

```bash
conda create -n nextflow -c bioconda nextflow
conda activate nextflow
nextflow -version
```

### Option 3: Manual Install

```bash
# Download specific version
wget https://github.com/nextflow-io/nextflow/releases/download/v23.10.0/nextflow
chmod +x nextflow
mv nextflow ~/bin/

# Add to PATH
echo 'export PATH=$PATH:~/bin' >> ~/.bashrc
source ~/.bashrc
```

## Container Setup

### Docker Setup

```bash
# Install Docker (Ubuntu)
sudo apt-get update
sudo apt-get install docker-ce docker-ce-cli containerd.io

# Add user to docker group
sudo usermod -aG docker $USER

# Test Docker
docker run hello-world
```

### Singularity Setup (HPC)

```bash
# On CentOS/RHEL
sudo yum install singularity

# On Ubuntu
sudo apt-get install singularity-container

# Test Singularity
singularity --version
```

### Pre-pulling Containers

```bash
# Docker
docker pull fred2/optitype:1.3.5
docker pull ensemblorg/ensembl-vep:release_110.1
docker pull griffithlab/pvactools:4.1.1

# Singularity
singularity pull docker://fred2/optitype:1.3.5
singularity pull docker://ensemblorg/ensembl-vep:release_110.1
singularity pull docker://griffithlab/pvactools:4.1.1
```

## Reference Data

### VEP Cache

The pipeline requires the Ensembl VEP cache for variant annotation.

```bash
# Create cache directory
mkdir -p /path/to/vep_cache
cd /path/to/vep_cache

# Download GRCh38 cache (approximately 15 GB)
curl -O https://ftp.ensembl.org/pub/release-110/variation/indexed_vep_cache/homo_sapiens_vep_110_GRCh38.tar.gz

# Extract
tar -xzf homo_sapiens_vep_110_GRCh38.tar.gz

# Download reference FASTA
curl -O https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
gunzip Homo_sapiens.GRCh38.dna.toplevel.fa.gz

# Index FASTA
samtools faidx Homo_sapiens.GRCh38.dna.toplevel.fa
```

### VEP Plugins

pVACseq requires specific VEP plugins:

```bash
# Clone VEP plugins
cd /path/to/vep_cache
git clone https://github.com/Ensembl/VEP_plugins.git

# Required plugins:
# - Frameshift
# - Wildtype
```

## NetMHCpan Setup

NetMHCpan requires an academic license from DTU.

### 1. Obtain License

1. Visit https://services.healthtech.dtu.dk/software.php
2. Register for an academic account
3. Request licenses for:
   - NetMHCpan 4.1
   - NetMHCIIpan 4.0
4. Download the archives

### 2. Install NetMHCpan

```bash
# Create installation directory
mkdir -p /opt/netmhc
cd /opt/netmhc

# Extract NetMHCpan
tar -xzf netMHCpan-4.1.Linux.tar.gz
cd netMHCpan-4.1

# Edit configuration
# Open 'netMHCpan' script and set:
#   NMHOME=/opt/netmhc/netMHCpan-4.1
#   TMPDIR=/tmp

# Download and extract data
cd /opt/netmhc/netMHCpan-4.1
wget https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/data.tar.gz
tar -xzf data.tar.gz

# Test installation
./netMHCpan -h
```

### 3. Install NetMHCIIpan

```bash
cd /opt/netmhc

# Extract NetMHCIIpan
tar -xzf netMHCIIpan-4.0.Linux.tar.gz
cd netMHCIIpan-4.0

# Edit configuration similar to NetMHCpan

# Download data
wget https://services.healthtech.dtu.dk/services/NetMHCIIpan-4.0/data.tar.gz
tar -xzf data.tar.gz

# Test
./netMHCIIpan -h
```

### 4. Add to PATH

```bash
echo 'export PATH=$PATH:/opt/netmhc/netMHCpan-4.1:/opt/netmhc/netMHCIIpan-4.0' >> ~/.bashrc
source ~/.bashrc
```

## HLA-HD Setup

HLA-HD requires manual installation.

### 1. Download HLA-HD

Request HLA-HD from: https://www.genome.med.kyoto-u.ac.jp/HLA-HD/

### 2. Install

```bash
# Extract
tar -xzf hlahd.version.tar.gz
cd hlahd.version

# Install dependencies
# HLA-HD requires: bowtie2, python3

# Run install script
sh install.sh

# Add to PATH
echo 'export PATH=$PATH:/path/to/hlahd/bin' >> ~/.bashrc
```

### 3. Download Reference Data

```bash
# Download IPD-IMGT/HLA database
cd /path/to/hlahd
./fetch_dictionary.sh
```

## Testing Installation

### Quick Test

```bash
# Clone pipeline
git clone https://github.com/YOUR_USERNAME/neoantigen-pipeline.git
cd neoantigen-pipeline

# Run with test profile
nextflow run main.nf -profile test,docker

# Or with Singularity
nextflow run main.nf -profile test,singularity
```

### Verify Components

```bash
# Check Nextflow
nextflow -version

# Check containers
docker images | grep -E "optitype|vep|pvactools"

# Check NetMHCpan
netMHCpan -h
netMHCIIpan -h

# Check VEP cache
ls /path/to/vep_cache/homo_sapiens/110_GRCh38/
```

## Troubleshooting

### Java Version Issues

```bash
# Check Java version
java -version

# Install Java 11 if needed
sudo apt install openjdk-11-jdk
```

### Container Pull Failures

```bash
# For Docker rate limits
docker login

# For Singularity, use cache
export SINGULARITY_CACHEDIR=/path/to/cache
```

### Memory Issues

```bash
# Increase Nextflow memory
export NXF_OPTS='-Xms1g -Xmx4g'
```

## Next Steps

After installation:

1. Create a samplesheet with your data
2. Configure paths in `nextflow.config`
3. Run a test sample
4. Scale to full cohort

See [Usage Guide](usage.md) for detailed instructions.
