# Part 1: CodingRNA Database Construction

## 1) Data Collection
- Collect  translatome data of rice varieties, including both conventional ,resistant and susceptible. (the same applies to other species.)

## 2) Data Download and Conversion
- Use `fastq-dump` and `fasterq-dump` to download SRA data and convert it into **FASTQ** format for further analysis.
```bash
## downloda data
aria2 url
wget -c url

## conversion data
fasterq-dump --split-3 SRR_ID -O Output_folder/ -e 16 -m 1000

## compress
pigz -p 16 *fastq
 
```

## 3) Quality Control and Removal of Low-Quality Sequences and Adapters

- Perform quality control using `fastp` to remove low-quality reads and adapter sequences from the data.
```bash
## Quality control
fastqc -t 40 SRR_ID.fastq.gz -o output_folder/ 2> SRR_ID.log

## Integration
multiqc output_folder/

## Remove adapters and low-quality sequences.
ls *fastq.gz|cut -d"_" -f 1|sort -u | while read id;

do
fastp \
--detect_adapter_for_pe \
--in1 ./${id}_1.fastq.gz \
--in2 ./${id}_2.fastq.gz \
--out1 ./1_pair_fastp/${id}_R1.fastq.gz \
--out2 ./1_pair_fastp/${id}_R2.fastq.gz \
--html ./1_pair_fastp/${id}_fastp.html \
-w 16 \
-z 7 \
-q 30 \
-u 20 \
-c \
-n 4 \
2> ./1_pair_fastp/${id}_fastp.log
done

## Second Remove adapters and low-quality sequences.
for id in {95..98}
do
trim_galore --gzip -q 20 -j 20 --length 15 --stringency 3 -e 0.1 --paired \
-o /gss1/home/hzhao/20230411/coding_data/SRR/unzip_SRR/paired_reads/1_pair_fastp \
/gss1/home/hzhao/20230411/coding_data/SRR/unzip_SRR/paired_reads/1_pair_fastp/SRR92180${id}_R1.fastq.gz \
/gss1/home/hzhao/20230411/coding_data/SRR/unzip_SRR/paired_reads/1_pair_fastp/SRR92180${id}_R2.fastq.gz \
&> SRR92180${id}_trimGal.log
done
```


## 4) Data Removal of rRNA

- Use `Bowtie2` to align and remove rRNA sequences from the data to ensure only non-rRNA data is used for transcript assembly.
```bash
for id in $(ls *.fastq.gz | cut -d'_' -f 1); do
    bowtie2 \
    --very-sensitive-local \
    --no-unal \
    --no-head \
    -I 1 -X 1000 -p 28 \
    -x /2_reference_rRNA \
    -1 ./${id}_R1.fastq.gz \
    -2 ./${id}_R2.fastq.gz \
    --un-conc-gz ../3_rmRNA/${id}_clean.fastq.gz \
    2> ${id}_MaprRNAStat.xls

    mv *.xls ../3_rmRNA/
done
```


## 5) Alignment to Reference Genome

- Align the RNA-seq data to the **rice reference genome (MSU v7)** using `Hisat2`.

```bash

    ls *gz|cut -d"_" -f 1|sort -u | while read id;
    do
    indexPrefix=/2_reference_genome/rice7
    hisat2 -p 28 \
    -x $indexPrefix \
    -1 ./${id}_clean.fastq.1.gz \
    -2 ./${id}_clean.fastq.1.gz \
    | samtools view -bS - \
    | samtools sort -O bam -@ 28 -o ./4_alignment/${id}_sorted.bam -
    2> ${id}.log
    mv *.log ./4_alignment/
    done
```
    

## 6) Transcript Assembly using StringTie and TACO Fusion

- Use `StringTie` for transcript assembly and `TACO` for transcript fusion to improve the completeness of the assemblies.
    
```bash
data=/input_folder/;
    for file in $data/*.bam;
    do
    # 获取文件名的前缀，去掉后缀和路径
    prefix=$(basename $file .bam)
    # 运行stringtie，指定输出目录和文件名
    stringtie \
    -p 28 -v \
    -o /output_folder/$prefix.gtf $file
    done
    
taco_run -p 8 \
-o ./TACO_results \
./coding_gtf.list \
--ref-genome-fasta /reference_genome.fa \
--filter-splice-juncs    

```
    

## 7) Comparison with MSU v7 Reference

- Use `GFFcompare` to compare assembled transcripts against the MSU v7 annotation to assess the quality and accuracy of the assembly.
    
```bash

gffcompare \
-r /codingRNA.gtf \
-o MSU_codingRNA_Compare \
codingRNA_merged.gtf
    
```
    

## 8) Merge with MSU v7 and Generate Final CodingRNA Database

- Merge the new transcript annotations with the MSU v7 database and generate the **final CodingRNA** database.
```bash

stringtie --merge -p 12 \
-G codingRNA_merged.gtf \
-o codingRNA_MSU.gtf \
-l Coding \
/MSU_codingRNA_gtf.list

```


### Final Output:

The final **CodingRNA** database includes high-quality protein-coding gene annotations derived from integrated translatome data, significantly improving the accuracy and completeness of the rice genome annotation. This database serves as the foundation for identifying and filtering potential **lncRNAs** in the subsequent pipeline.




---
