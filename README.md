# RiceLncRNA: An Optimized Pipeline for Rice Long Non-coding RNA Identification

## Project Overview

This project provides a comprehensive pipeline for lncRNA identification and analysis, from CodingRNA database construction to lncRNA identification and analysis. The pipeline integrates multiple bioinformatics tools to achieve full-process analysis from raw data processing to functional annotation.


This is the English version of the README. You can view the Chinese version of the README by clicking the link below:
[Chinese README](https://github.com/njausxl/RiceLncRNA/edit/main/README.md)
 

The project consists of four main modules:
- Module 001: Basic data processing environment for raw data QC, alignment, and assembly
- Module 002: LncRNA identification environment focusing on non-coding RNA prediction and screening
- Module 003: Functional analysis environment for expression analysis and functional annotation
- Module 004: Basic data processing environment (snakemake one-click run) for raw data QC, alignment, and assembly

⚠️ Important Notes:
- This pipeline involves multiple programming languages including Python, R, and Perl
- Integrates over 30 professional bioinformatics software tools
- Due to complex dependencies, one-click installation and usage is not possible
- Recommended to install and run each module step by step
- Verify output results before proceeding to the next step

## Features
- CodingRNA database construction
- Raw data quality control and preprocessing
- Reference genome-based transcript assembly
- Strict screening process for long non-coding RNA
- Transcript functional annotation and analysis
- Support for multi-sample parallel processing
- Integration of multiple professional bioinformatics tools

## System Requirements

- Linux/Unix operating system
- Python ≥ 3.7
- R ≥ 4.0
- Conda package manager

## Dependencies

Major bioinformatics tools:
- HISAT2 (transcriptome alignment)
- StringTie (transcript assembly)
- TACO (transcript integration)
- FastQC (sequencing data QC)
- Samtools (SAM/BAM file processing)
- MultiQC (QC report integration)
- Bowtie2 (transcriptome alignment)
- GFFread (GFF/GTF file conversion)
- Gffcompare (transcript assembly evaluation)
- Fastp (raw data QC)
- TrimGalore (raw data QC)
- seqkit (sequence processing)
- transeq (sequence translation)
- pfam_scan.pl (Pfam database alignment)
- cmscan (Rfam database alignment)
- Diamond (NR database alignment)
- Snakemake (workflow management)
- PLEK (lncRNA identification)
- cnci (lncRNA identification)
- cpc2 (lncRNA identification)
- FeatureCounts (expression quantification)
- DESeq2 (differential expression analysis)
- wgcna (WGCNA analysis)
- clusterProfiler (functional enrichment analysis)
- gseGO (GO enrichment analysis)
- gseKEGG (KEGG enrichment analysis)
- enrichplot (enrichment result visualization)
- ggplot2 (result visualization)

## Installation Guide

1. Create three independent conda environments:

### Environment 001: Basic Data Processing Environment
```bash
conda create -n lncrna_001 python=3.7
conda activate lncrna_001
conda install -c bioconda hisat2 stringtie taco fastqc samtools multiqc bowtie2 gffread gffcompare fastp trim-galore seqkit emboss
```

### Environment 002: LncRNA Identification Environment
```bash
conda create -n lncrna_002 python=3.7
conda activate lncrna_002
conda install -c bioconda plek cnci cpc2
conda install -c bioconda diamond pfam_scan hmmer infernal
```

### Environment 003: Functional Analysis Environment
```bash
conda create -n lncrna_003 r=4.0
conda activate lncrna_003
conda install -c bioconda bioconductor-deseq2 bioconductor-wgcna bioconductor-clusterprofiler
conda install -c conda-forge r-ggplot2
conda install -c bioconda subread  # for FeatureCounts
```

2. Create working directories:
```bash
# Create main directory
mkdir -p ~/rice_lncrna_pipeline
cd ~/rice_lncrna_pipeline

# Create raw data directory
mkdir -p raw_data/{sra,fastq}

# Create QC-related directories
mkdir -p qc/{fastqc,multiqc}

# Create alignment and assembly directories
mkdir -p alignment/{hisat2_index,hisat2_output,stringtie_output}
mkdir -p assembly/{taco_output,merged_transcripts}

# Create lncRNA identification directory
mkdir -p lncrna/{plek,cnci,cpc2,merged_results}

# Create functional analysis directory
mkdir -p analysis/{expression,deseq2,wgcna,enrichment}

# Create results and log directories
mkdir -p results
mkdir -p logs
```

3. Download reference genome:
```bash
cd ~/rice_lncrna_pipeline
wget http://rice.plantbiology.msu.edu/pub/data/Eukaryotic_Projects/o_sativa/annotation_dbs/pseudomolecules/version_7.0/all.dir/Osativa_204_v7.0.fa.gz
gunzip Osativa_204_v7.0.fa.gz
mv Osativa_204_v7.0.fa reference/
cd alignment/hisat2_index
hisat2-build ../../reference/Osativa_204_v7.0.fa rice7_index
```

## Analysis Pipeline

### 1. CodingRNA Database Construction (Environment 001)
- Data collection and download
- Data quality control
- rRNA sequence removal
- Reference genome alignment
- Transcript assembly and merging
- Comparison with reference genome annotation
- Final fusion with original genome annotation, database generation

### 2. RNA-seq Data Analysis (Environment 002)

#### 2.1 Data Collection and Download
- Collect published rice RNA-seq data
- Download raw data using SRA toolkit
- Convert to FASTQ format



#### 2.2 Data Quality Control
- Quality assessment using FastQC
- QC report integration using MultiQC
- Low-quality sequence and adapter removal using Fastp

#### 2.3 rRNA Sequence Removal
- Bowtie2 alignment to rRNA database
- Non-rRNA sequence extraction
- Quality reassessment

#### 2.4 Reference Genome Alignment
- Alignment using HISAT2
- SAM/BAM file conversion and sorting using SAMtools
- Generate alignment statistics report

#### 2.5 Transcript Assembly and Merging
- Transcript assembly using StringTie
- Merge GTF files from multiple samples
- Generate non-redundant transcript set

#### 2.6 Reference Genome Annotation Comparison
- Compare annotation files using GFFcompare
- Extract known coding genes
- Construct coding RNA database

### 3. Non-coding RNA Database Construction (Environment 002)

#### 3.1 Data Preprocessing
- Extract transcripts with class codes i,u,x,o,p
- Length filtering (>200nt)
- Exon number filtering

#### 3.2 Coding Potential Assessment
- CPC2 analysis
- CNCI analysis
- PLEK analysis
- Rfam/Pfam/NR database alignment
- Integration of all software identification results

#### 3.3 LncRNA Classification
- Long non-coding RNA identification
- Classification (antisense, intronic, intergenic, etc.)
- Integration with known lncRNA database information

### 4. Functional Analysis and Annotation (Environment 003)

#### 4.1 Expression Analysis
- Calculate expression levels in each sample
- Expression normalization
- Sample correlation analysis

#### 4.2 Differential Expression Analysis
- DESeq2 differential expression analysis
- Differential gene screening

#### 4.3 Co-expression Network Analysis
- WGCNA co-expression network construction
- Module identification and analysis
- Key gene screening

#### 4.4 Functional Enrichment Analysis
- GO functional enrichment analysis
- KEGG pathway enrichment analysis
- Functional annotation visualization

#### 4.5 Individual lncRNA GSEA Analysis
- Functional enrichment analysis
- Result visualization

## Precautions

### Environment Configuration
1. Ensure all dependencies are correctly installed
2. Check software version compatibility
3. Properly allocate computational resources

### Data Processing
1. Regularly backup important data
2. Record key parameter settings
3. Save analysis log files

### Result Interpretation
1. Pay attention to biological replication
2. Focus on data quality control results
3. Validate reliability of key findings


## RiceLncRNA Repository

### 1. Overview

This repository contains essential datasets and tools for the analysis of long non-coding RNAs (lncRNAs) in rice, focusing on their roles in plant immunity.

### 2. Dataset Contents

#### 2.1 OSA rRNA
- **File**: `OSA_rRNA.fa`
- **Description**: Contains the ribosomal RNA (rRNA) sequences of rice (*Oryza sativa*). This file is used for filtering rRNA contamination during RNA-seq preprocessing.

#### 2.2 OSA Known LncRNA
- **File**: `OSA_Known_LncRNA.zip`
- **Description**: A zip file containing known lncRNA sequences from multiple lncRNA databases specifically for rice. It is useful for the identification and annotation of lncRNAs in rice.

#### 2.3 Rice7 Genomic Data
- **File**: `rice7_all_fa.zip`
- **Description**: This archive contains the rice genome sequence (MSU v7) in FASTA format. It serves as a reference for genome-wide RNA-seq alignment and other genomic analyses.

- **File**: `rice7_all_gtf.zip`
- **Description**: Contains the MSU v7 gene annotation in GTF format. It is essential for the functional annotation of RNA-seq data.

#### 2.4 Org RiceDB (GO/KEGG Enrichment)
- **File**: `org_riceDB_AH94060.zip`
- **Description**: This archive contains the rice-specific annotation database (OrgDb) used for performing GO and KEGG enrichment analyses. It provides necessary gene annotations for functional analysis of rice-related data.

### 3. Dataset Usage

These datasets provide fundamental resources for transcriptomic and functional analyses, such as the identification of lncRNAs, alignment of RNA-seq data, and gene annotation. You can use these files for a variety of bioinformatics analyses, such as:

- **rRNA Filtering**: Use `OSA_rRNA.fa` to filter out rRNA sequences during RNA-seq preprocessing.
- **LncRNA Identification**: Use `OSA_Known_LncRNA.zip` to compare with newly identified lncRNAs in your RNA-seq dataset.
- **Genomic Alignment**: Align your RNA-seq reads to the reference genome found in `rice7_all_fa.zip`.
- **Gene Annotation**: Annotate your aligned reads using the `rice7_all_gtf.zip` file.
- **GO/KEGG Enrichment**: For functional enrichment analysis, `org_riceDB_AH94060.zip` is used to provide necessary gene annotations for rice.

### 4. Acknowledgements

The datasets and annotations used in this repository have been generously provided by various contributors. These datasets are crucial for the research and development of the RiceLncRNA pipeline and facilitate the study of lncRNAs in rice.

## License
This repository is licensed under the MIT License. See the LICENSE file for more details.

## Citation

If you use this pipeline, please cite the following paper:

Shan, X., Xia, S., Peng, L., Tang, C., Tao, S., Baig, A., & Zhao, H. (2024). Identification of Rice LncRNAs and Their Roles in the Rice Blast Resistance Network Using Transcriptome and Translatome. *Plants* (Under Review).

## Contact Information
- Issue Feedback: [GitHub Issues]
- Email Contact: [sdausxl@126.com]

## Change Log
- 2024-02-17: Initial Version V1.0 Release
- Upcoming: Continuous pipeline optimization
- Upcoming: Addition of new analysis tools
```
