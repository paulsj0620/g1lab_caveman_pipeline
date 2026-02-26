# 🧬 Caveman Variant Calling Pipeline
**Somatic SNV Calling Workflow (CaVEMan-based)**

---

## 📌 Overview

This repository provides a reproducible and modular Snakemake-based somatic variant calling pipeline built around the CaVEMan (Cancer Variants through Expectation Maximisation) algorithm.

The workflow is designed for:

- Tumor–Normal paired analysis
- Whole-genome sequencing (WGS)
- High-depth somatic SNV detection
- HPC cluster execution

The pipeline emphasizes reproducibility, scalability, and clean rule modularization.

---

## 🏗 Workflow Architecture

FASTQ / BAM  
↓  
Alignment (if required)  
↓  
Pre-processing (sorting, indexing, QC)  
↓  
CaVEMan Somatic SNV Calling  
↓  
Post-processing & Filtering  
↓  
Final VCF Output

---

## 📂 Repository Structure

.
├── Snakefile \
├── config/ \
│   └── config.yaml \
├── modules/ \
│   ├── alignment.snakefile \
│   ├── caveman.snakefile \
│   ├── filtering.snakefile \
│   └── utils.snakefile \
├── scripts/ \
├── logs/ \
└── benchmarks/ \

---

## ⚙️ Requirements

Core Software:

- Python ≥ 3.7
- Snakemake ≥ 6
- CaVEMan
- Samtools
- BWA (if alignment included)
- GATK (optional)
- Conda (recommended)

Install Snakemake:

conda install -c bioconda snakemake

---

## 🧪 Input Requirements

Typical inputs:

- Tumor BAM (or FASTQ)
- Matched Normal BAM (or FASTQ)
- Reference genome (indexed)
- Associated reference files required by CaVEMan

Example configuration snippet:

samples:
  sample1:
    tumor: path/to/tumor.bam
    normal: path/to/normal.bam

reference:
  fasta: path/to/reference.fa

---

## 🚀 Quick Start

1️⃣ Configure Samples

Edit config/config.yaml and define tumor-normal pairs and reference paths.

2️⃣ Dry Run

snakemake -np

3️⃣ Execute Pipeline

Local execution:

snakemake --cores all

Cluster execution:

snakemake --cluster-config cluster.json --jobs 200 --max-jobs-per-second 5

---

## 📊 Outputs

- Somatic SNV VCF files
- Intermediate processed BAM files
- Log files (logs/)
- Benchmark files (benchmarks/)

---

## 🖥 HPC Support

Designed for SLURM/SGE-like environments.

Recommended options:

--rerun-incomplete  
--latency-wait 60  
--use-conda

---

## 🧠 Design Philosophy

- Modular Snakemake rules
- Clear separation of alignment and calling
- Tumor–Normal explicit pairing
- Scalable for large WGS cohorts
- Production-ready HPC execution

---

## 🛠 Customization

You may customize:

- Filtering thresholds in filtering.snakefile
- Resource allocation in cluster.json
- Reference configuration in config.yaml

---

## 📜 License

Add your preferred license here.

---

## 👨‍🔬 Maintainer

Developed for research use.  
Maintained by: SungJoon
