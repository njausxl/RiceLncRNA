# RiceLncRNA
## 1. conda
conda create --name CodingRNA_database_pipeline python=3.7

conda activate CodingRNA_database_pipeline

conda install -c bioconda taco
## 2. Create folders
mkdir -p ~/CodingRNA_database_pipeline/
cd ~/CodingRNA_database_pipeline
mkdir -p sra
mkdir -p fastq
mkdir -p hisat2_index
mkdir -p fastqc
mkdir -p taco
mkdir -p hisat2_index

## 3. Download genome files    
wget http://rice.plantbiology.msu.edu/download/Genome_MSU_r7.tar.gz
tar -xvf Genome_MSU_r7.tar.gz
hisat2-build Genome_MSU_r7.tar.gz hisat2_index/rice7_index

## 4. snakemake
vim snakemake_coding.py

```python
workdir:"~/CodingRNA_database_pipeline/"

### 4.1 Get all sample names
import glob

#### Add a variable to specify file format
file_format = "sra"  # can be "sra" or "fastq.gz"

if file_format == "sra":
    samples = [x.split('/')[-1].split('.')[0] for x in glob.glob('sra/*.sra')]
elif file_format == "fastq.gz":
    samples = [x.split('/')[-1].split('_')[0] for x in glob.glob('fastq/*.fastq.gz')]



#### Specify specific samples for analysis
# samples=["Sample1","Sample2","Sample3"]

#### 4.2 Activate environment
onstart:
    shell("""
        source ~/miniconda3/bin/activate CodingRNA_database_pipeline
    """)

#### 4.3 Check if folder exists
rule check_folder:
    output:
        "path/to/your/folder/"
    shell:
        "[ -d {output} ] || mkdir -p {output}"

#### 4.4 Total output files
rule all:
    input:
        expand("hisat2/{sample}_sorted.bam", sample=samples),
        "multiqc_report/multiqc_report.html",
        "multiqc_report2/multiqc_report.html",
        expand("stringtie/{sample}.gtf",sample=samples)



def check_samples(samples):
    """Check if sample names are correct"""
    for sample in samples:
        if sample not in samples:
            raise ValueError(f"Sample name {sample} does not exist")

### 4.5 Start running
if file_format == "sra":
#### 4.5.1 fasterq-dump
    rule fasterq_dump:
        input:
            "sra/{sample}.sra"
        output:
            r1 = "read1/{sample}_1.fastq",
            r2 = "read1/{sample}_2.fastq"
        threads: 21
        run:
            import os.path as path
            output_dir = path.dirname(output.r1)
            shell(
                f"fasterq-dump.3.0.2 -e {threads} -f --split-files {input} -O {output_dir}"
            )

#### 4.5.2 pigz compress
    rule pigz_compress:
        input:
            r1 = "read1/{sample}_1.fastq",
            r2 = "read1/{sample}_2.fastq"
        output:
            r1 = "read1/{sample}_1.fastq.gz",
            r2 = "read1/{sample}_2.fastq.gz"
        threads: 21
        shell:
            """
            pigz -7 -p {threads} {input.r1}
            pigz -7 -p {threads} {input.r2}
            """

#### 4.5.3 fastqc
rule fastqc:
    input:
        r1="read1/{sample}_1.fastq.gz",
        r2="read1/{sample}_2.fastq.gz"
    output:
        html1="qc/fastqc/{sample}_1_fastqc.html",
        zip1 ="qc/fastqc/{sample}_1_fastqc.zip",
        html2="qc/fastqc/{sample}_2_fastqc.html",
        zip2 ="qc/fastqc/{sample}_2_fastqc.zip"
    params:
        extra="--quiet"
    log:
        "log/fastqc/{sample}.log"
    run:
            shell("""
                fastqc --threads {threads} {params.extra} --outdir qc/fastqc {input.r1} {input.r2}
            """)


#### 4.5.4 hisat2

## --rna-strandness FR  # Forward data
## --rna-strandness RF  # Reverse data

rule hisat2:
    input:
        r1="bowtie2/{sample}_clean.1.fastq.gz",
        r2="bowtie2/{sample}_clean.2.fastq.gz",
        idx=expand("hisat2_index/rice7_all.{i}.ht2", i=range(1,9))
    output:
        bam="hisat2/{sample}_sorted.bam"
    log:
        "log/hisat2/{sample}.log"
    threads: 26
    run:
        import os
        os.environ['TMPDIR'] = '~/temp/'
        idx_prefix = os.path.commonprefix(input.idx).rstrip(".1.ht2")

        shell(
            f"""
            hisat2 \
            --threads {threads} \
            -x {idx_prefix} \
            -1 {input.r1} \
            -2 {input.r2} \
            | samtools sort -O bam -@ {threads} -o {output.bam} -
            2> {log}
            """
        )
#### 4.5.5 stringtie

rule stringtie:
    input:
        bam="hisat2/{sample}_sorted.bam"
    output:
        gtf="stringtie/{sample}.gtf"
    params:
        threads=28
    shell:
        """
        stringtie \
        -p {params.threads} -v \
        --rf \
        -o {output.gtf} \
        {input.bam}
        """


#### 4.5.6. taco_run

rule taco_run:
    input:
        gtf=expand("stringtie/pe/{sample}.gtf", sample=SAMPLES)
    output:
        directory("taco_run/pe/")
    params:
        threads=8
    conda:
        "/gss1/home/hzhao/20230411/2023-09-17_translatome/py2_environment.yaml"
    shell:
        """
        taco_run -p {params.threads} \
        -o {output} \
        {input.gtf}
```
## 开始运行snakemake
snakemake --latency-wait 60 -ps snakemake.py --cores 22

