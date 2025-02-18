# Part 2: lncRNA Identification

In this section, we focus on identifying long non-coding RNAs (lncRNAs). We assume that all datasets and reference files have been prepared or generated in the previous step (Part 1: CodingRNA Database Construction).

---

## 1) Data Collection
- Gather **strand-specific RNA-seq** datasets for various rice samples (including resistant and susceptible varieties).
- Ensure these data cover the time points or treatments that will be analyzed for lncRNA expression.


## 2) Data Download and Conversion
- Obtain raw SRA files and convert them into FASTQ format using `fastq-dump` or `fasterq-dump`.
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

- Use `fastp` to perform adapter trimming and remove low-quality reads.
    
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
    --out1 ./${id}_R1.fastq.gz \
    --out2 ./${id}_R2.fastq.gz \
    --html ./${id}_fastp.html \
    -w 16 \
    -z 7 \
    -q 30 \
    -u 20 \
    -c \
    -n 4 \
    2> ./${id}_fastp.log
    done


  
## Second Remove adapters and low-quality sequences.
for id in {95..98}
do
trim_galore --gzip -q 20 -j 20 --length 15 --stringency 3 -e 0.1 --paired \
-o /fastp \
/SRR${id}_R1.fastq.gz \
/SRR${id}_R2.fastq.gz \
&> SRR${id}_trimGal.log
done
   
```
    

## 4) Removal of rRNA

- Align to an rRNA reference and remove matching reads to retain only non-rRNA sequences.
    
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
        --un-conc-gz ../rmRNA/${id}_clean.fastq.gz \
        2> ${id}_MaprRNAStat.xls
    done
    ```
    

## 5) Alignment to Reference Genome

- Align clean reads to the **rice reference genome** (e.g., MSU v7 or IRGSP) using `Hisat2`.
    
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
    
    > **Note**: Use the appropriate `--rna-strandness` option (e.g., `RF` or `FR`) based on the library preparation.
    

## 6) Transcript Assembly and Merging

- Assemble transcripts for each sample using `StringTie`, then optionally merge them (e.g., via `TACO`) if needed for lncRNA discovery.
    
    ```bash
data=/input_folder/;
        for file in $data/*.bam;
        do
        prefix=$(basename $file .bam)
        # 运行stringtie，
        stringtie \
        -p 28 -v \
        -o /output_folder/$prefix.gtf $file
        done
        
    taco_run -p 8 \
    -o ./TACO_results \
    ./lncRNA_gtf.list \
    --ref-genome-fasta /reference_genome.fa \
    --filter-splice-juncs   
    ```
    

## 7) Using the CodingRNA Database to Extract Class Codes

- Compare the assembled transcripts with the **CodingRNA** reference (constructed in Part 1) to label transcripts by class code.
    
```bash

## 1. compare gtf
gffcompare \
-r /codingRNA.gtf\
-o LncRNA_codingRNA_Compare \
../LncRNA_merged.gtf    


## 2. extract i u p x o 
awk '$3=="i" || $3=="u" || $3=="x" || $3=="o" || $3=="p" {print "\"" $5 "\""}' LncRNA_codingRNA_Compare.gtf.tmap > LncRNA_iuxop.id

## 3. Merge gtf annotated files 
LC_ALL=C fgrep -f LncRNA_iuxop_MSU.id ../01_gtf/2024-03-18_LncRNA_merged.gtf > LncRNA_ixoup.gtf

## 4. Extract the sequence of the lncRNA from the gtf file
gffread LncRNA_ixoup.gtf -g /reference_genome.fa -w LncRNA_ixoup.fa

```
    
    > **Note**: Focus on transcripts labeled **i, x, o, u, p** for potential lncRNAs.
    

## 8) Filtering with Pfam, Rfam, and NR Databases

- Remove transcripts that align to known protein domains or non-coding RNA families:
    1. **Pfam** to filter out transcripts with known protein domains.
    2. **Rfam** to remove tRNAs, rRNAs, snRNAs, snoRNAs, etc.
    3. **NR** (via `Diamond`) to exclude transcripts with significant similarity to known proteins.
```bash
# Pfam

## 1. Translated into amino acid sequence
source activate gffcompare

transeq ../LncRNA_ixoup.fa LncRNA_ixoup.pep3.fa -frame=3

## 2. Remove sequences larger than 100k
seqkit seq -m 0 -M 100000 -i LncRNA_ixoup.pep3.fa > LncRNA_ixoup.pep3_100k.fa

## 3. Sequence to remove pfam database
pfam_scan.pl \
-cpu 24 \
-fasta LncRNA_ixoup.pep3.fa \
-dir ~/lncRNA_analysis/2_reference/pfamdata \
-outfile LncRNA_ixoup_pfam.out

# Rfam

cmscan --cpu 24 \
-Z 720 \
--cut_ga --rfam --nohmmonly --fmt 2 \
--tblout LncRNA_ixoup_rfam.tblout \
-o LncRNA_ixoup_Rfam.result \
--clanin ~/Rfam.cm \
LncRNA_ixoup.fa

# NR

## 1. database building
diamond makedb --in nr.fa -p 24 --db nr_diamond

## 2. align
diamond blastx \
--db /nr_diamond \
--query LncRNA_ixoup.fa \
-e 1e-5 --outfmt 6 \
--more-sensitive \
--max-target-seqs 1 \
--threads 16 \
--quiet \
--id 80 \
--subject-cover 50 \
--query-cover 50 \
--out ./LncRNA_ixoup_NR_matches.txt

## 3. Filter align results
awk 'BEGIN { FS = "\t" } $11 < 1e-10 && $3 > 90 && $4 > 100 && $12 > 100 && ($5/($4-$5-$6)) <= 0.05 && ($6/($4-$5-$6)) <= 0.02 { print $1 }' LncRNA_ixoup_NR_matches.txt | sort | uniq > LncRNA_ixoup_NR_matches_noncoding_ids.txt


```


## 9) Coding Potential Assessment (CNCI, PLEK, CPC2)

- Evaluate the remaining transcripts using three coding-potential prediction tools:
    
    1. **CNCI**
    2. **PLEK**
    3. **CPC2**
    
    Only transcripts predicted as **non-coding** by **all three** tools are retained. For example, to run **CPC2**:
    
```bash

# CNCI
python2 CNCI.py \
-f $INPUT_FILE \
-o $OUTPUT_DIR \
-m pl \
-p 10

# CPC2
~/CPC2_standalone-1.0.1/bin/CPC2.py -i $INPUT_FILE -o $OUTPUT_FILE 2>cpc2.log

# PLEK
PLEK.py -fasta $INPUT_FILE -out $OUTPUT_FILE -thread 10

```
    

## 10) Final Intersection of Non-coding Transcripts

- Take the intersection of transcripts that pass both **step 8** (Pfam/Rfam/NR) and **step 9** (non-coding by all tools).  
    These represent **high-confidence lncRNAs**.

```bash
## 1. cpc2
grep 'noncoding' cpc2.out | awk '{print $1}' > cpc2_msu_id.txt
wc -l cpc2_msu_id.txt
15175 cpc2_msu_id.txt

  
## 2. plek   
grep -w 'Non-coding' plek.out | awk '{print $3}' | sed 's/>//g' > plek_msu_id.txt
 
wc -l plek_msu_id.txt   
13513 plek_msu_id.txt

## 3. cnci
awk '$2=="noncoding"{print $1}' .cnci.out > cnci_msu_id.txt  

wc -l cnci_msu_id.txt  
12261 cnci_msu_id.txt


## 4. The intersection of three software.
cat *txt | sort | uniq -c | awk '{if($1==3){print}}' | wc -l
10875

cat *txt | sort | uniq -c | awk '{if($1==3){print $2}}' > 3_noncoding_msu.id

## 5. The identification results

### Three software identifications ultimately resulted in x pending lncRNAs.


### Quotation marks need to be added.
sed 's/^/"/;s/$/"/' 3_noncoding_msu.id > 3_noncoding_msuv7.id


LC_ALL=C fgrep -w -f 3_noncoding_msuv7.id LncRNA.gtf > LncRNA_noncoding3.gtf

## Verify whether the number of transcripts is equal.
cut -f 3 LncRNA_noncoding3.gtf | grep "transcript" | wc -l

## 6. Organize the identification results from the Rfam database.

perl infernal-tblout2gff.pl --cmscan --fmt2 LncRNA_Rfam.tblout > LncRNA_Rfam.gff3

### 6.1 Filter out sequences with 1e-5 and transcript lengths greater than 200.

#### This saves all the results.
awk -F'\t' '($5-$4+1 > 200 && $9+0 < 1e-5) {print}' LncRNA_Rfam.gff3 > LncRNA_Rfam.id

#### This is only the ID that needs to be removed from rfam.
awk -F'\t' '($5-$4+1 > 200 && $9+0 < 1e-5) {print $1}' LncRNA_Rfam.gff3 | sort | uniq > LncRNA_Rfam.id
wc -l LncRNA_Rfam.id




## 7. Organize the results of the Pfam database.

grep -v '^#' LncRNA_pfam.out | grep -v '^\s*$' | awk '($13 < 1e-5){print $1}' | awk -F "_" '{print $1}' | sort | uniq > LncRNA_pfam.ID

wc -l LncRNA_pfam.ID


## 8. We need to further remove the transcripts from the nr database.


### Statistic the number of NR database comparisons.
#### Number before filtering:
cut -f 1 LncRNA_NR_matches.txt | uniq | wc -l


#### Filtering and comparing the results.
awk 'BEGIN { FS = "\t" } $11 < 1e-10 && $3 > 90 && $4 > 100 && $12 > 100 && ($5/($4-$5-$6)) <= 0.05 && ($6/($4-$5-$6)) <= 0.02 { print $1 }' LncRNA_NR_matches.txt | sort | uniq > LncRNA_NR_matches.id
#### The number after filtering:
1075 LncRNA_NR_matches.id


#### 9. Integrate Rfam-Pfam-NR cat together and reorder them.

cat ../002_pfam/LncRNA_pfam.ID ../003_rfam/LncRNA_Rfam.id ../005_nr/LncRNA_NR_matches.id | sort | uniq > LncRNA_NR-Rfam-Pfam.id

wc -l LncRNA_NR-Rfam-Pfam.id
2596 LncRNA_NR-Rfam-Pfam.id


#### 10. Remove these IDs from the LncRNA_noncoding3.gtf file.

LC_ALL=C fgrep -v -w -f LncRNA_NR-Rfam-Pfam.id LncRNA_noncoding3.gtf > LncRNA_ixoup_msuv7_noncoding3-rfam-pfam-nr.gtf

cut -f 3 LncRNA_ixoup_msuv7_noncoding3-rfam-pfam-nr.gtf | grep "transcript" | wc -l


#### 11. Extract transcripts.
gffread LncRNA_ixoup_msuv7_noncoding3-rfam-pfam-nr.gtf -g ~/lncRNA_analysis/2_reference/all.con.fa -w LncRNA_ixoup_msuv7_noncoding3-rfam-pfam-nr.fa

grep "^>" LncRNA_ixoup_msuv7_noncoding3-rfam-pfam-nr.fa | wc -l



```

## 11) lncRNA Classification

- Classify the final lncRNAs based on their genomic positions relative to protein-coding genes:
    
    - **Intergenic lncRNAs (lincRNAs)**
    - **Antisense lncRNAs**
    - **Intronic lncRNAs**
    - **Bidirectional** or **divergent** lncRNAs
    - **Sense overlapping** lncRNAs
    
    > This classification often leverages class code labels (i, x, o, u, p) or other annotation overlap strategies.
```bash

## 1. Identify based on the classcode results of gffcompare.
x reverse lncRNA

i Intronic lncRNAs

p,u  (including Intergenic and bidirectional lncRNA) - bidirectional = Intergenic numbers

Bidirectional lncRNA 	
### 2. Continue to identify bidirectional lncRNAs.
    
# 2.1  First, go to tmap to extract the ids of p and u, and then use the ids to extract gtf.
    
awk '$3=="u" || $3=="p" {print "\"" $5 "\""}' LncRNA_codingbase_type.LncRNA_ixoup_rfampfam_noncoding3_nr.gtf.tmap > LncRNA_up_codingbase.id
    
wc -l 3218 LncRNA_up_codingbase.id
    
# 2.2 Extract GTF file based on ID
LC_ALL=C fgrep -f LncRNA_up_codingbase.id ../01_gtf/2024-03-18_LncRNA_merged.gtf > LncRNA_up_codingbase.gtf

# Count the number of transcripts in the extracted GTF file
cut -f 3 LncRNA_up_codingbase.gtf | grep "transcript" | wc -l  # Expected result: 

# 2.3 Sort the lncRNA GTF file
sort -k1,1 -k4,4n LncRNA_up_codingbase.gtf > sorted_LncRNA_up_codingbase.gtf

# Sort the genome annotation GTF file
sort -k1,1 -k4,4n ~/rice7_genome/rice7_all.gtf > ~/rice7_genome/sorted_rice7_all.gtf

# Use bedtools to find the closest features
bedtools closest \
  -a sorted_LncRNA_up_codingbase.gtf \
  -b ~/rice7_genome/sorted_rice7_all.gtf \
  -S \
  -D b \
  -id > bid_filter_type.txt

# Filter the table
# 1. First, keep only the rows with "transcript" in the third column. Both lncRNA and mRNA need to be filtered.
# 2. Then, filter for values between -1000 and -1 in the 19th column.
# 3. The final result will have 460 entries.

awk -F'\t' '$3 == "transcript" && $12 == "transcript" && $19 >= -1000 && $19 <= -1' bid_filter_type.txt > bid_after_filter_type.txt

# Count the number of rows in the filtered file
wc -l  # Expected result: 

    ```


## 12 ) Identify known lncRNA

- The obtained results were compared with multiple public rice lncRNA databases (PlantNATdb, PNRD, RNAcentral, NONCODE, CANTATAdb, GreeNC) using an E-value of <1e-10, sequence consistency of >80%, and coverage of >50%.
```bash

blastn \
-db OSA_Known_LncRNA \
-evalue 1e-10 \
-num_threads 18 \
-max_target_seqs 5 \
-query LncRNA.fa \
-outfmt '6 qseqid sseqid pident qcovs length mismatch gapopen qstart qend sstart send qseq sseq evalue bitscore' \
-out blastn_LncRNA.txt

```



---

### Final Output

After this pipeline, you will have a reliable **lncRNA dataset** for rice:

- Unambiguously **non-coding** transcripts validated by multiple databases (Pfam, Rfam, NR) and coding-potential tools (CNCI, PLEK, CPC2).
- Systematically **classified** based on genomic context.
- Ready for **functional analyses** in the next pipeline phase.
