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


## RiceLncRNA Repository

This repository contains various datasets and resources related to the identification and functional analysis of long noncoding RNAs (lncRNAs) in rice (Oryza sativa). All data here are based on **Oryza sativa Japonica** (commonly known as Japonica rice). If you plan to use other rice varieties (e.g., Indica rice), please carefully cross-check the data for compatibility.

## Data Files

Below is a list of files in this repository and their descriptions:

### 1. **s001-OSA_Known_LncRNA.zip**
   - Contains known rice lncRNAs from multiple databases. This file is crucial for the identification and analysis of rice lncRNAs.

### 2. **s002-OSA_rRNA.fa**
   - This file contains the ribosomal RNA (rRNA) sequences for **Oryza sativa**. These sequences are used as references in RNA-seq analyses to filter out rRNA from the data.

### 3. **s003-rice7_all_gtf.zip**
   - Contains the GTF (General Feature Format) file of the **MSU V7 rice genome annotation**. This is a reference genome annotation file for use in RNA-seq analysis and transcript annotation.

### 4. **s004-rice7_all_fa.zip**
   - Contains the **MSU V7 rice genome sequence** in FASTA format. This file includes all genomic sequences that are used for alignment during transcriptome analysis.

### 5. **s005-RiceSymbolGoAnnotation.xlsx**
   - This file contains **symbol IDs and Gene Ontology (GO) annotations** for rice genes. It helps in the functional annotation of genes and the analysis of biological functions using GO terms.

### 6. **s006-org_riceDB_AH94060.zip**
   - A zip file containing the **org.db package for rice**. This package is used for gene enrichment analysis (GO, KEGG, etc.) and provides key rice gene information in the appropriate format for enrichment tools.

### 7. **s007-OryzabaseGeneListGOKEGGSymbolTraitsRice.xlsx**
   - This file contains detailed information about **rice gene GO and KEGG annotations**, as well as **symbol IDs and trait information**. It is essential for understanding gene function and performing further analysis such as pathway enrichment.

### 8. **s008-ALL_oryza_sativa_hormone.xlsx**
   - Contains information about **hormones in rice**. It includes data relevant to plant hormone signaling pathways and their effects on rice physiology.

### 9. **s009-ANNA_2024-03-23-rice-NLR.fa**
   - This file contains **NLR (Nucleotide Binding Leucine-Rich Repeat) genes** sequences in rice. NLR genes play critical roles in immunity and are important for studying rice disease resistance.

### 10. **s010-Osj_pep.fas.gz**
   - Contains **protein sequences** for rice. This file is important for downstream analysis, including protein structure prediction and functional studies.

### 11. **s011-Osj_cds.fas.gz**
   - Contains **CDS (coding sequences)** for rice. This file includes the translated gene sequences for functional annotation and gene expression analysis.

### 12. **s012-20230225_qPCR_R.md**
   - A **R script** for calculating qPCR results. This script is useful for analyzing quantitative PCR data and conducting gene expression studies.

### 13. **s013-20230221_qPCR_Python.py**
   - A **Python script** for calculating qPCR results. This Python-based tool provides an alternative to the R script for analyzing qPCR data.

### 14. **s014-Osj_TF_list.xlsx**
   - Contains a list of **transcription factor (TF) IDs** for rice. This file is important for understanding the regulation of gene expression in rice and is used in transcription factor-based analyses.

## Usage Instructions

1. **Data Processing:**
   - Use the rice genome annotations and lncRNA databases (found in `s001`, `s003`, `s004`) for transcript assembly and RNA-seq analysis.
   - The ribosomal RNA file (`s002`) is used for filtering out rRNA from RNA-seq data before lncRNA identification.

2. **Gene Enrichment Analysis:**
   - For functional annotation and enrichment analysis, refer to the **GO**, **KEGG**, and **symbol ID** annotations provided in `s005` and `s007`.
   - For enrichment analysis, use the `s006` package, which supports GO and KEGG enrichment tools.

3. **Hormone Pathway Analysis:**
   - Use `s008` for the analysis of plant hormones and their pathways, which are critical for studying rice immunity and other physiological traits.

4. **NLR and Transcription Factor Studies:**
   - For disease resistance and regulatory studies, refer to the **NLR gene sequences** in `s009` and the **transcription factors** in `s014`.

5. **qPCR Analysis:**
   - Analyze qPCR data using the provided R (`s012`) or Python (`s013`) scripts.

```
