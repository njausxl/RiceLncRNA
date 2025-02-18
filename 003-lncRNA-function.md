# Part 3: Functional Analysis and Screening of lncRNAs

In this section, we focus on functional analysis of the identified lncRNAs. This includes quantification, differential expression analysis, target gene prediction, and enrichment analysis.

---

## 1) lncRNA and mRNA Quantification
- Use **FeatureCounts** to quantify both lncRNAs and mRNAs. This step ensures we have the expression levels of both types of RNAs.
  ```bash
# Set the input folder path
  input_folder=./input_folder
  
  # Set the output folder path
  output_folder=./output_folder
  
  # Loop through all .bam files in the input folder
  for file in $input_folder/*.bam
  do
    # Get the prefix of the filename, remove the path and extension
    prefix=$(basename $file .bam)
  
    # Run featureCounts, specify the output directory and filename
    featureCounts -T 18 \
      -p \
      -s 2 \
      -t exon \
      -g transcript_id \
      -a coding_gtf \
      -o $output_folder/${prefix}_lncRNA.txt \
      $file
  done
  
```

## 2) Differential Expression Analysis with DESeq2

- Perform differential expression analysis using **DESeq2**. This will identify differentially expressed lncRNAs (DELs) across different conditions or time points.
### DESeq 2
```R
# ################################################################### 
# ############### 1. Set paths, load R packages and datasets #########################
# ###################################################################

# 1.1 Set working directory

rm(list=ls())

# sink()
########### ！！！！！Uncomment to run second time ！！！！#######
########### ！！！！！Uncomment to run second time ！！！！#######
########### ！！！！！Uncomment to run second time ！！！！#######

# sink("*.txt", type="output")

# Set working directory
setwd("./working_directory")

# Load necessary libraries
library(DESeq2)
library(ggplot2)

# Set file paths
file_path_example <- "example_group/example_file.txt" # 1

file_path <- file_path_example

# Read data
data <- read.table(file_path, header = TRUE, row.names = 1)

# Parse file name to get sample names
file_name <- basename(file_path)
matches <- regmatches(file_name, regexec("(.+)_vs_(.+).txt", file_name))
treated_sample <- matches[[1]][2] # Treatment group sample name
untreated_sample <- matches[[1]][3] # Control group sample name

# Get column names
col_names <- colnames(data)

# Initialize condition vector
condition_labels <- rep("", length(col_names))

# Fill condition vector
for (i in seq_along(col_names)) {
  if (grepl(treated_sample, col_names[i])) {
    condition_labels[i] <- "treated"
  } else if (grepl(untreated_sample, col_names[i])) {
    condition_labels[i] <- "untreated"
  }
}

# Convert condition to factor and ensure "untreated" comes first, "treated" second
condition <- factor(condition_labels, levels = c("untreated", "treated"))

# Print condition vector to confirm
print(condition)

## Batch effect
# batch <- factor(c(rep("batch1", 3), rep("batch2", 3)))

head(condition)
length(condition)

# ## 1.5 Create a dataframe called coldata

coldata <- data.frame(row.names = colnames(data), condition)

colnames(data)

# ## 1.6 Experimental design matrix

# design <- model.matrix( ~condition+batch, data=coldata)

## Batch effect ## 

design <- model.matrix( ~condition, data=coldata)

head(design)

# ## 1.7 Check data for rank deficiency, if data is linearly correlated, it can't be used
#library(caret)

# linearCombos <- findLinearCombos(design)
# View(linearCombos)

# design_clean <- design[, -3]

# ## 1.8 Build expression matrix, store RNA-seq count data and sample information.
## Remove rank deficiency ## 

dataExpr_deseq <- DESeqDataSetFromMatrix(countData = data,
                                         colData = coldata,
                                         design = design)

### 1.9 Filter low-expressed data
### Keep genes that have counts >= 10 in at least 2 samples and total count >= 50
keep <- rowSums(counts(dataExpr_deseq) >= 10) >= 2 &
  rowSums(counts(dataExpr_deseq)) >= 50

ddsHTSeq <- dataExpr_deseq[keep,]
dim(ddsHTSeq)
#[1] 1949   6
dim(dataExpr_deseq)
#[1] 9003   6


###################################################################
############### 2. Perform differential expression analysis ###########################
###################################################################

## 2.1 Differential expression analysis

dds <- DESeq(ddsHTSeq)

## 2.2 Check sizeFactor
## SizeFactor near 1 is ideal, if it's too large, it should be removed.
# Check size factors, keep samples with size factors between 0.5 and 1.5

dds$sizeFactor

dds_filtered <- dds[, dds$sizeFactor >= 0.2 & dds$sizeFactor <= 2.5]

# Check filtered samples
colnames(dds_filtered)

# Run DESeq2 differential expression analysis

# Get results
res <- results(dds_filtered)
# res <- results(dds)

resultsNames(dds)

sum(res$padj < 0.05, na.rm=TRUE)
#[1] 201

colData(dds_filtered)
# colData(dds)
## The variables and test information used can be found by calling mcols on the result object

mcols(res)$description
summary(res)

design(dds_filtered)


## 2.6 Save results
# Save results to CSV file
## Create a general variable name

# Extract file name from file path
file_name <- basename(file_path)

# Use regular expression to extract specific name format (example_name_vs_example.txt)
matches <- regmatches(file_name, regexec("(.+)_vs_(.+).txt", file_name))
base_name <- paste(matches[[1]][2], "vs", matches[[1]][3], sep="_")  # Create variable

# File name for saving padj results
padj_file_name <- paste0("00-", base_name, "_lncRNA_padj.csv")

# When you need to save different results, just add different suffixes to the base name
# For example, saving padj results

write.csv(as.data.frame(res), file = padj_file_name)


## 2.7 Convert results to dataframe, and save normalized count files  [Not needed for repeated use]
# Extract normalized read counts

normalized_counts <- counts(dds_filtered, normalized=TRUE)

# File name for saving normalized data
normalized_file_name <- paste0("01-", base_name, "_lncRNA_normalized_counts.csv")

write.csv(normalized_counts, file=normalized_file_name,
          row.names=TRUE)

resdata <- as.data.frame(res)

## padj < 0.05 
resdata$padj.signif <- ifelse(resdata$padj < 0.05, 
                              "padj < 0.05", 
                              "padj >= 0.05")

sig_genes <- subset(resdata, padj < 0.05 & abs(log2FoldChange) >= 1.5)
### Check the number of genes with p < 0.05 and significant expression

nrow(sig_genes)

## 2.10 Count upregulated and downregulated genes | log2FC | >= 2 

num_upregulated <- sum(sig_genes$log2FoldChange >=1.5)
num_downregulated <- sum(sig_genes$log2FoldChange < -1.5)

### 2.11 Output the number of upregulated and downregulated genes
print(paste("Number of significantly upregulated genes: ", num_upregulated))
# 80
print(paste("Number of significantly downregulated genes: ", num_downregulated))
# 98 

### 2.12 Add new column indicating if the gene is upregulated or downregulated

sig_genes$direction <- ifelse(sig_genes$log2FoldChange >= 1.5, "up", "down")


## 2.13 Write significant differential genes to CSV file
# File name for saving significant genes

sig_file_name <- paste0("02-", base_name, "_lncRNA_sig_padj.csv")

write.csv(sig_genes, file= sig_file_name, row.names = TRUE)


###################################################################
############### 3. Generate plots for differential expression ###########################
###################################################################

## 3.1 Volcano plot
### Create a volcano plot, x-axis is log2 fold change, y-axis is -log10 adjusted p-adj,
### Point colors indicate whether adjusted p-value is less than 0.05
# File name for saving volcano plot
volcano_file_name <- paste0("03-", base_name, "_lncRNA_volcano_padj.pdf")

pdf(volcano_file_name, width = 8, height = 6.5)

resdata <- na.omit(res)  ## Needed when plotting with padj
resdata <- as.data.frame(resdata)
label = subset(resdata, padj < 0.05 & abs(log2FoldChange) > 1.5)
label1 = rownames(label)

Significant=ifelse((resdata$padj < 0.05 & abs(resdata$log2FoldChange) > 1.5), 
                   ifelse(resdata$log2FoldChange > 1,"Up","Down"), "No Change")

# Create volcano plot, ensuring symmetry
# Calculate maximum log2FoldChange absolute value
max_log2FC <- max(abs(resdata$log2FoldChange))

# Create symmetric volcano plot
ggplot(resdata, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(col = Significant)) +
  scale_color_manual(values = c("#0072B5", "grey", "#BC3C28")) +
  labs(title = base_name) +
  geom_vline(xintercept = c(-1.5, 1.5), colour = "black", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), colour = "black", linetype = "dashed") +
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold")) +
  labs(x = "log2(Fold Change)", y = "-log10(p-adj)") +
  theme(axis.text = element_text(size = 13), axis.title = element_text(size = 13)) +
  theme_bw() +
  coord_cartesian(xlim = c(-max_log2FC, max_log2FC)) # Use maximum log2FoldChange absolute value to set symmetric x-axis range
dev.off()

################################################################################
### 3.2 Results Interpretation
## log2 fold change (MLE): condition treated【treated】 vs untreated【control】 
## Wald test p-value: condition treated vs untreated 
# baseMean: Mean expression value
# log2FoldChange: log2 fold change of treatment vs control
# lfcSE: Estimated standard deviation of log2FoldChange
# stat: Test statistic used
# pvalue: p-value
# padj: Adjusted p-value after multiple testing correction
###############################################################################

## 3.3 MA plot for overall results
## y-axis is log2 fold change, x-axis is mean of normalized counts. This plot visualizes the distribution of differential expressed genes.
## Function plotMA shows log2 fold change against the normalized count average from DESeqDataSet across all samples.
## Points are shown as red if adjusted p-value < 0.1. Points outside the window are plotted as up or down triangles.
library(ggplotify)

# Generate the plot using plotMA() function
# Save the plot as a PNG file
# File name for saving MA plot
MA_file_name <- paste0("04-", base_name, "_lncRNA_mA_padj.pdf")

pdf(MA_file_name)
plotMA(res, alpha = 0.05,
       xlab = "Mean of normalized Counts",
       ylab = "Log2 Fold Change", colSig = "red")
dev.off()


## 3.4 Sample distance heatmap

## dist function is applied to the transposed normalized count matrix to get the sample-to-sample distance
# dist function calculates Euclidean distance between these vectors to get the sample distance matrix.
# Calculate sample distance
vsd <- vst(dds_filtered, blind=FALSE)

# Calculate sample distance for the raw count matrix
sampleDists <- dist(t(assay(vsd)))

# Perform average linkage hierarchical clustering on the distance matrix
tree <- hclust(sampleDists, method = 'average')

# Open PDF device to save the plot as a PDF file
# File name for saving hclust plot
hclust_file_name <- paste0("05-", base_name, "_lncRNA_hclust_padj.pdf")
pdf(file= hclust_file_name, w=10, h=7)

# Plot the dendrogram with no x-axis labels and subtitle, adjust text size
plot(tree, xlab="", sub="", cex=0.7, main="Sample Clustering")

# Add a red line on the dendrogram as a reference for clustering
abline(h=120, col="red")

# Close PDF device to save the plot
dev.off()


## 3.5 PCA plot

# returnData=TRUE to return the data for plotting
# percentVar is the percentage of variance explained by each principal component

pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
# Create PCA plot, adjust text positions and size
pcaPlot <- ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  geom_text(aes(label=name), nudge_x = 2, nudge_y = 3, size=3) +  # Adjust position and font size
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed()

# Display the plot
print(pcaPlot)

# Save the plot as a PDF file
# File name for saving PCA plot
PCA_file_name <- paste0("06-", base_name, "_lncRNA_PCA_padj.pdf")
ggsave(PCA_file_name , plot = pcaPlot, 
       width = 7, height = 7, device = "pdf")
dev.off()


#### 3.10 Gene expression heatmap
### Normalized gene expression levels

# signii <- subset(resdata, pvalue< 0.05 & abs(log2FoldChange) > 1.5) # Earlier it was res and padj
# sig_genes2
allsig <- merge(normalized_counts, sig_genes, by=0)

x <- length(dds$sizeFactor)

sigCount <- allsig[, 2:(x + 1)]

library(pheatmap)

## File name for saving heatmap
heatmap_file_name <- paste0("07-", base_name, "_lncRNA_heatmap_padj.pdf")
pdf(heatmap_file_name)

# Generate heatmap
pheatmap(log2(sigCount + 1), scale = 'row', # row
         show_rownames = F,
         treeheight_col = 1,
         treeheight_row = 0)

# Close plot device
dev.off()

###########################################################


# sink()
########### ！！！！！Uncomment to run second time ！！！！#######
########### ！！！！！Uncomment to run second time ！！！！#######
########### ！！！！！Uncomment to run second time ！！！！#######
################################### End ##################################
################################## End #################################

```

### 2. Statistical analysis of the number of differentially expressed genes.
```PYTHON
import os
import pandas as pd
import glob

# Change directory to the working folder
os.chdir("./working_directory")

# Define the file path pattern to match
file_pattern = 'example_folder/*_sig_padj.csv'

# Get all files matching the pattern
file_list = glob.glob(file_pattern)

# Initialize a list to store the results
results = []

# Iterate through each file, read it, and count 'up' and 'down' values in the last column
for file_path in file_list:
    df = pd.read_csv(file_path)
    last_col_name = df.columns[-1]  # Get the last column's name
    up_count = (df[last_col_name] == 'up').sum()  # Count 'up' values
    down_count = (df[last_col_name] == 'down').sum()  # Count 'down' values
    file_name = file_path.split('/')[-1]  # Extract the file name from the path
    results.append([file_name, up_count, down_count])  # Add the result to the results list

# Convert the results list to a DataFrame
results_df = pd.DataFrame(results, columns=['FileName', 'Upregulated', 'Downregulated'])

# Define the output file path
output_file = 'expression_summary.csv'

# Save the results DataFrame to a CSV file
results_df.to_csv(output_file, index=False)

# Output the file path for download
output_file

```

## 3) Correlation analysis

```bash
# ###################################################################
# ############### Step 1: Copy DEG and DEL Gene ID Files #############
# ###################################################################

# Define the base path where your files are located
base_path="input_folder"

# Define the directories for DEG and DEL gene ID files
DEG_directory="${base_path}DEG-Re-ID"
DEL_directory="${base_path}DEL-Re-ID"

# Print message to indicate the start of the copying process
echo "Copying DEG and DEL gene ID files to their respective directories."

# Copy DEG and DEL files to their respective directories
# Modify the source paths based on where your files are located
cp /path/to/DEG_files/*.csv $DEG_directory/
cp /path/to/DEL_files/*.csv $DEL_directory/

# ###################################################################
# ############### Step 2: Remove Double Quotes from Gene IDs #########
# ###################################################################

# Go to the DEL folder and remove quotes from the CSV files
cd $DEL_directory
echo "Removing double quotes from DEL gene IDs files in ${DEL_directory}."
for file in *.csv; do
    sed -i 's/\"//g' "$file"
done

# Go to the DEG folder and remove quotes from the CSV files
cd $DEG_directory
echo "Removing double quotes from DEG gene IDs files in ${DEG_directory}."
for file in *.csv; do
    sed -i 's/\"//g' "$file"
done

# ###################################################################
# ############### Step 3: Run R Script to Extract High Correlation Genes
# ###################################################################

# Print message for running R script (You need to execute this separately in R)
echo "Please run the R script to extract genes with correlation > 0.5."
echo "R script: Cor.R"

# Example R command (ensure R script is set up correctly):
# Rscript Cor.R

# ###################################################################
# ############### Step 4: Count Unique High-Correlation Gene Pairs ####
# ###################################################################

# Go to the folder where the correlation results are saved
echo "Counting the unique high-correlation gene pairs in all_high_Re_correlations.tsv."
cd $base_path

# Count unique gene pairs (removes duplicates and counts)
cat all_high_Re_correlations.tsv | sort | uniq | wc -l

# ###################################################################
# ############### Step 5: Extract Specific lncRNA Correlation Results
# ###################################################################

# Extract correlation results for a specific lncRNA (e.g., "LncRNA.13491.1")
echo "Extracting correlations for LncRNA.13491.1 from all_high_Re_correlations.tsv."
grep "LncRNA.13491.1" all_high_Re_correlations.tsv | sort -k3,3nr

# ###################################################################
# ############### Final Message #####################################
# ###################################################################

echo "Process completed successfully."
echo "You can now find the extracted and processed files in the respective folders."

```

### Cor. R
```R
# Clear environment variables
rm(list = ls())

# Set working directory
setwd("./working_directory")  # Example directory

# Load required packages
library(data.table)
library(dplyr)
library(pheatmap)
library(ggplot2)

# Define function to read and process data
process_correlation <- function(mrna_lncrna_file, deg_file, del_file, output_file) {
  # Read mRNA and lncRNA count matrix
  mRNA_lncRNA_counts <- fread(mrna_lncrna_file, key = "gene")
  
  # Read DEG and DEL gene lists
  mRNA_list <- readLines(deg_file)
  lncRNA_list <- readLines(del_file)
  
  # Extract matching gene expression matrices
  matching_lncRNA_matrix <- mRNA_lncRNA_counts[gene %in% lncRNA_list, ]
  matching_mRNA_matrix <- mRNA_lncRNA_counts[gene %in% mRNA_list, ]
  
  # Combine the matrices
  mRNA_lncRNA_matrix <- rbind(matching_lncRNA_matrix, matching_mRNA_matrix)
  
  # Convert data.table to data.frame
  mRNA_lncRNA_matrix_df <- as.data.frame(mRNA_lncRNA_matrix)
  
  # Group by gene ID and calculate the mean for each gene
  merged_df <- mRNA_lncRNA_matrix_df %>%
    group_by(gene) %>%
    summarise(across(everything(), mean, na.rm = TRUE))
  
  # Convert tibble to data.frame
  merged_df <- as.data.frame(merged_df)
  
  # Save gene IDs for row and column name settings
  gene_ids <- merged_df$gene
  
  # Remove the gene ID column
  merged_df <- merged_df[, -1]  # Removed the gene ID column after saving it
  
  # Transpose and log-transform
  merged_df_t <- t(merged_df)
  merged_df_t <- log(merged_df_t + 1)
  
  # Calculate correlation
  merged_df_tt <- cor(merged_df_t)
  
  # Use saved gene IDs for row and column names
  rownames(merged_df_tt) <- gene_ids
  colnames(merged_df_tt) <- gene_ids
  
  # Add gene name as part of the data to ensure correct alignment of column names during output
  merged_df_tt <- cbind(Gene = rownames(merged_df_tt), merged_df_tt)
  
  # Save the result with row and column names, using the newly added column as a marker
  write.table(merged_df_tt, output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

# Get all DEG and DEL file names
deg_files <- list.files("DEG-Re-ID/", pattern = "*.csv", full.names = TRUE)  # Example folder
del_files <- list.files("DEL-Re-ID/", pattern = "*.csv", full.names = TRUE)  # Example folder

# Ensure that DEG and DEL files match up
for (i in seq_along(deg_files)) {
  deg_file <- deg_files[i]
  del_file <- del_files[i]  # Ensure DEL and DEG files are in the correct order or sorted by name
  
  # Get file name and construct output file name
  base_name <- sub("mRNA_sig_padjDEG.id.csv", "", basename(deg_file))
  output_file <- sprintf("correlation_results_%s.tsv", base_name)
  
  # Run processing function
  process_correlation("./WGCNA_example/2024-3-27_normal_count_lncrna_mrna.tsv",  # Example file
                      deg_file, del_file, output_file)
}

### Extract high correlation gene pairs greater than 0.5

# Set working directory to the folder containing the correlation result files
setwd("./working_directory/")  # Example folder

# Load required package
library(data.table)

# Get a list of all correlation result files
cor_files <- list.files(pattern = "correlation_results.*\\.tsv")

# Initialize an empty list to store high correlation results for each file
high_cor_list <- list()

# Set high correlation threshold
threshold <- 0.5

# Process each file
for (file in cor_files) {
  # Read correlation matrix
  cor_data <- fread(file, data.table = FALSE)
  
  # Set the first column as row names
  rownames(cor_data) <- cor_data$Gene
  cor_data <- cor_data[ , -1]
  
  # Convert correlation data to matrix
  cor_matrix <- as.matrix(cor_data)
  
  # Extract correlation values above the threshold
  high_cor_indices <- which(abs(cor_matrix) > threshold, arr.ind = TRUE)
  high_cor_data <- cor_matrix[high_cor_indices]
  
  # Create a data frame to store high correlation gene pairs
  high_cor_df <- data.frame(
    gene1 = rownames(cor_matrix)[high_cor_indices[, 'row']],
    gene2 = colnames(cor_matrix)[high_cor_indices[, 'col']],
    correlation = high_cor_data,
    stringsAsFactors = FALSE
  )
  
  # Save high correlation gene pairs for each file
  output_filename <- sub("correlation_results", "high_cor", file)
  write.table(high_cor_df, output_filename, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  
  # Add data frame to the list
  high_cor_list[[file]] <- high_cor_df
}

# Combine all high correlation results
all_high_cor <- do.call(rbind, high_cor_list)

# Save the combined results
write.table(all_high_cor, "all_high_Re_correlations.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


```
## 4) Cis Target Gene Prediction for DELs

- Use **BedTools** or custom scripts to predict the cis-regulatory target genes of differentially expressed lncRNAs (DELs). This typically involves identifying DEGs located within 20 kb upstream or downstream of DELs, and their correlation is greater than 0.5. 
- *(According to the size of rice TAD: The median size of TAD is 35 kb, covering 69.7% of the rice genome, while the median size of the regions located outside the identified TAD (inter-TAD) is 25 kb, covering 30.3% of the rice genome. (The specific data depends on the TAD of your species.)).*
### Cis -Target genes
```bash

# Step 1: Extract DEG and DEL gene IDs from the result files
# Modify this command with actual DEG and DEL result files

# Example: extracting DEG gene IDs from a CSV or other file
awk -F',' '{print $1}' DEG_results.csv > DEG_gene_ids.txt   # Extract gene IDs from DEG results

# Example: extracting DEL gene IDs from a CSV or other file
awk -F',' '{print $1}' DEL_results.csv > DEL_gene_ids.txt   # Extract gene IDs from DEL results

# Step 2: Extract mRNA GTF using DEG gene IDs
# The fgrep command is used to extract lines that match the DEG gene IDs from the mRNA GTF file
fgrep -f DEG_gene_ids.txt /path/to/mRNA_annotation.gtf > mRNA_DEG.gtf   # Modify the file path accordingly

# Step 3: Extract lncRNA GTF using DEL gene IDs
# Use fgrep to extract lines that match the DEL gene IDs from the lncRNA GTF file
fgrep -f DEL_gene_ids.txt /path/to/lncRNA_annotation.gtf > lncRNA_DEL.gtf   # Modify the file path accordingly

# Step 4: Use bedtools for annotation (cis-regulatory analysis)
# Use the window function to find overlapping regions of lncRNA and mRNA annotations with a specified window size
bedtools window \
  -a lncRNA_DEL.gtf \
  -b mRNA_DEG.gtf \
  -l 20000 -r 20000 > cis_DE-lncmRNA.txt   # Annotation with window size of 20,000 bp on both sides

# Output results
echo "Annotation completed. Results saved in cis_DE-lncmRNA.txt."


```
## 5) Trans Target Gene Prediction for DELs

- Predict trans-target genes by using **RIblast** and calculating correlation between DELs and potential target DEGs.
```bash


# Set base path
base_path="/share/home/hnu_panyix/data/sxl-trans-2024-05-14/"

# Print all DEG files to confirm their existence
echo "Listing all DEG files in the directory:"
ls ${base_path}/*_sig_padjDEG.fa

# Loop through each DEL file
for del_file in ${base_path}/*_sig_padjDEL.fa; do
    # Use basename to get the file name
    base_name=$(basename "${del_file}" "_sig_padjDEL.fa")

    # Find the corresponding DEG file
    deg_file="${base_path}/${base_name}_mRNA_sig_padjDEG.fa"

    # Check if the DEG file exists
    if [[ -f "${deg_file}" ]]; then
        echo "Found DEG file: ${deg_file}"
        
        # Create the database for RIblast
        RIblast db -i "${deg_file}" -o "${deg_file%.fa}.db" -r 1

        # Perform prediction using RIblast
        RIblast ris -i "${del_file}" -o "${base_path}/${base_name}_RIblast.out" -d "${deg_file%.fa}.db" -s 0
    else
        echo "DEG file not found for ${base_name}"
    fi
done


```

## 6) venn Diagram of Stat lncRNA Expression Across Different Groups

- Use **VennDiagram** in python to visualize the overlap of differentially expressed lncRNAs across multiple conditions or time points.
- We recommend directly using this online website ( https://www.bioinformatics.com.cn/static/others/jvenn/example.html ).
```python

import os
import pandas as pd
import numpy as np

# Set the working directory
os.chdir("./working_directory")  # Example directory

# Define function to read and process data
def process_correlation(mrna_lncrna_file, deg_file, del_file, output_file):
    # Read mRNA and lncRNA count matrix
    mRNA_lncRNA_counts = pd.read_csv(mrna_lncrna_file, sep="\t", index_col="gene")
    
    # Read DEG and DEL gene lists
    mRNA_list = pd.read_csv(deg_file, header=None).squeeze().tolist()
    lncRNA_list = pd.read_csv(del_file, header=None).squeeze().tolist()
    
    # Extract matching gene expression matrices
    matching_lncRNA_matrix = mRNA_lncRNA_counts.loc[mRNA_lncRNA_counts.index.isin(lncRNA_list), :]
    matching_mRNA_matrix = mRNA_lncRNA_counts.loc[mRNA_lncRNA_counts.index.isin(mRNA_list), :]
    
    # Combine the matrices
    mRNA_lncRNA_matrix = pd.concat([matching_lncRNA_matrix, matching_mRNA_matrix])
    
    # Group by gene and calculate the mean for each gene
    merged_df = mRNA_lncRNA_matrix.groupby('gene').mean()
    
    # Log-transform and transpose the matrix
    merged_df_t = np.log1p(merged_df.T)
    
    # Calculate correlation
    correlation_matrix = merged_df_t.corr()
    
    # Save the result
    correlation_matrix.to_csv(output_file, sep="\t", header=True, index=True)

# Get all DEG and DEL file names
deg_files = [os.path.join("DEG-Re-ID/", f) for f in os.listdir("DEG-Re-ID/") if f.endswith(".csv")]  # Example folder
del_files = [os.path.join("DEL-Re-ID/", f) for f in os.listdir("DEL-Re-ID/") if f.endswith(".csv")]  # Example folder

# Ensure that DEG and DEL files match up
for deg_file, del_file in zip(deg_files, del_files):
    # Get file name and construct output file name
    base_name = deg_file.split("/")[-1].replace("mRNA_sig_padjDEG.id.csv", "")
    output_file = f"correlation_results_{base_name}.tsv"
    
    # Run processing function
    process_correlation("./WGCNA_example/2024-3-27_normal_count_lncrna_mrna.tsv",  # Example file
                        deg_file, del_file, output_file)

### Extract high correlation gene pairs greater than 0.5

# Set working directory to the folder containing the correlation result files
os.chdir("./working_directory/")  # Example folder

# Get a list of all correlation result files
cor_files = [f for f in os.listdir() if f.startswith("correlation_results") and f.endswith(".tsv")]

# Initialize an empty list to store high correlation results for each file
high_cor_list = []

# Set high correlation threshold
threshold = 0.5

# Process each file
for file in cor_files:
    # Read correlation matrix
    cor_data = pd.read_csv(file, sep="\t", index_col="Gene")
    
    # Convert correlation data to matrix
    cor_matrix = cor_data.to_numpy()
    
    # Extract correlation values above the threshold
    high_cor_indices = np.where(np.abs(cor_matrix) > threshold)
    high_cor_data = cor_matrix[high_cor_indices]
    
    # Create a DataFrame to store high correlation gene pairs
    high_cor_df = pd.DataFrame({
        "gene1": np.array(cor_data.index)[high_cor_indices[0]],
        "gene2": np.array(cor_data.columns)[high_cor_indices[1]],
        "correlation": high_cor_data
    })
    
    # Save high correlation gene pairs for each file
    output_filename = file.replace("correlation_results", "high_cor")
    high_cor_df.to_csv(output_filename, sep="\t", index=False, header=True)

    # Add data frame to the list
    high_cor_list.append(high_cor_df)

# Combine all high correlation results
all_high_cor = pd.concat(high_cor_list)

# Save the combined results
all_high_cor.to_csv("all_high_Re_correlations.tsv", sep="\t", index=False, header=True)


```
    

## 7) GO and KEGG Enrichment Analysis for Cis and Trans Target Genes

- Perform **GO** and **KEGG** enrichment analysis on the cis and trans target genes of the DELs. This helps to understand the biological roles of lncRNAs.

```R

################## 1. Load Data and R Packages ###########################

# Clear environment variables to avoid previous data affecting current analysis
rm(list=ls())

# Set stringsAsFactors to FALSE to prevent automatic conversion of strings to factors
options(stringsAsFactors = FALSE)

# Load required libraries
library(org.Osativa.eg.db)
library(clusterProfiler)
library(AnnotationHub)
library(AnnotationDbi)
library(riceidconverter)
library(ggplot2)
library(pathview)
library(DOSE)
library(topGO)
library(enrichplot)
library(GOplot)
library(tidyverse)
library(ggsci)
library(ggthemes)
library(viridis)
library(dplyr)
library(createKEGGdb)
library(GseaVis)
library(ReactomePA)

### 1.1 Read Original DEG Files

# Set working directory to the specific folder on your desktop for easy file management
setwd("./working_directory/")  # Example directory

# Read gene ID text file with headers and tab-separated format
filename <- "./example_folder/two_inter_21_DEL_target_gene_id.txt"  # Example file path

# Read the data
gene_info <- read.csv(filename, header=TRUE)

# Set column name for gene IDs
names(gene_info) <- c("MSU")

# Extract file name without extension
filename_no_ext <- sub("\\.txt$", "", basename(filename))

# Define regular expression pattern to remove prefix "02-" and suffix "_sig_padj"
pattern <- "^02-(.*)_sig_padj$"
matches <- sub(pattern, "\\1", filename_no_ext)

matches <- "two_inter_21_DEL_target_gene_id"

############################################################ 

### 1.4 Convert MSU ID to SYMBOL and RAP ID
# Use RiceIDConvert function to convert MSU ID to SYMBOL ID and RAP ID

msu_symbol_id <- RiceIDConvert(myID = gene_info, 
                               fromType = 'MSU', 
                               toType = 'SYMBOL')

# Ensure column names are correct
colnames(msu_symbol_id) <- c("MSU", "SYMBOL")

### 1.5 Load orgDB Database

rice <- loadDb("./example_folder/ricedb")  # Example database path

### 1.6 Extract SYMBOL ID
# Assign the converted SYMBOL ID to rapgene variable for further enrichment analysis
symbol_id <- as.character(msu_symbol_id[,2])

symbol_id <- unique(symbol_id)
length(symbol_id)

# Check available keytypes in rice database
keytypes(rice)

### 1.7 Convert SYMBOL to ENTREZID
# Use AnnotationDbi::select function to select required columns from rice database
gene_annotations <- AnnotationDbi::select(rice, 
                                          keys = symbol_id, 
                                          columns = c("ENTREZID", 
                                                      "SYMBOL", 
                                                      "GENENAME", 
                                                      "GOALL", 
                                                      "REFSEQ",
                                                      "ONTOLOGYALL",
                                                      "ALIAS"), 
                                          keytype = "SYMBOL")
dim(gene_annotations)

#### 1.8 Merge DataFrames
# Merge data frames based on SYMBOL
merged_data <- merge(gene_annotations,
                     msu_symbol_id, by = "SYMBOL")
dim(merged_data)

# Remove duplicate entries
merged_data <- unique(merged_data)
dim(merged_data)

# Remove rows with "#N/A" or NA values
merged_data <- merged_data[!apply(merged_data, 1, function(row) any(row == "#N/A")), ]
merged_data <- na.omit(merged_data)

dim(merged_data)

# Continue merging with gene_info
merged_data2 <- merge(merged_data,
                      gene_info, by = "MSU")
dim(merged_data2)

## 1.8 Save the merged DataFrame
# Generate file name for merged data
merged_data_filename <- paste0("001-", matches, "_annotation.tsv")

write.table(merged_data2, 
            file =  merged_data_filename,
            sep = "\t", 
            row.names = FALSE, 
            col.names = TRUE, 
            quote = FALSE)

### 1.9 Extract Gene Set and Genes
### Gene set used for GSEA and with log2FC, genes for KEGG/GO enrichment analysis input

# Extract ENTREZID and log2FoldChange columns, and remove duplicates
gene <- unique(merged_data[, c("ENTREZID")])

str(gene)

#############################################################

###################################################################
#
#######                2. GO Enrichment Analysis              ###########
#
###################################################################

### 2.1 Perform GO Enrichment Analysis

Go_enrichment <- enrichGO(gene,
                          OrgDb = rice, # Replace with the appropriate organism database
                          keyType = "ENTREZID", # Gene ID type
                          ont = "ALL", # Or "BP", "MF", "CC"
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.12, # Default 0.05
                          qvalueCutoff = 0.2)  # Default 0.2

# Generate dot plot for GO enrichment
Go_enrichment_dotplot <- dotplot(Go_enrichment, 
                                 showCategory = 10,
                                 split="ONTOLOGY") +
  facet_grid(ONTOLOGY~.,scale='free')

print(Go_enrichment_dotplot)

### 2.2 Save Enrichment Analysis Results
# Convert enrichment results to data frame and save to file
Go_enrichment_dataframe <- as.data.frame(Go_enrichment)

# Generate output file name
Go_enrichment_name <- paste0("002-", "_Cis_GO_enrichment-", matches, ".tsv")

write.table(Go_enrichment_dataframe, 
            file = Go_enrichment_name, 
            sep = "\t", 
            row.names = FALSE, 
            col.names = TRUE, 
            quote = FALSE)

#### 2.3 Plot GO Enrichment Dot Plot

# Calculate enrichment factor
Go_enrichment_tsv = read.table(Go_enrichment_name,
                               header=TRUE,
                               sep="\t",
                               quote = "",
                               fill = TRUE) 

Go_enrichment_tsv <- separate(data=Go_enrichment_tsv, 
                              col=GeneRatio,
                              into = c("GR1", "GR2"), 
                              sep = "/")  # Split GeneRatio into 2 columns (GR1, GR2)

Go_enrichment_tsv <- separate(data=Go_enrichment_tsv, 
                              col=BgRatio, 
                              into = c("BR1", "BR2"), 
                              sep = "/")  # Split BgRatio into 2 columns (BR1, BR2)

Go_enrichment_tsv <- mutate(Go_enrichment_tsv,
                            enrichment_factor = (as.numeric(GR1)/as.numeric(GR2))/(as.numeric(BR1)/as.numeric(BR2))) # Calculate Enrichment Factor 

### Save Enrichment Factor
go_enrichment_factor_filename <- paste0("003-", "_Cis_GO_enrichment-", matches, "_enrichmentFactor.tsv") 

write.table(Go_enrichment_tsv, 
            file = go_enrichment_factor_filename,
            sep = "\t", 
            col.names = TRUE,
            quote = FALSE,
            row.names = FALSE)

### 2.3.2 Plot GO Enrichment (Top 10 terms)

eGoBP <- Go_enrichment_tsv %>% 
  filter(ONTOLOGY=="BP") %>%
  filter(row_number() >= 1, row_number() <= 10)

eGoCC <- Go_enrichment_tsv %>% 
  filter(ONTOLOGY=="CC") %>%
  filter(row_number() >= 1, row_number() <= 10)

eGoMF <- Go_enrichment_tsv %>% 
  filter(ONTOLOGY=="MF") %>%
  filter(row_number() >= 1, row_number() <= 10)

eGo10 <- rbind(eGoBP, eGoMF, eGoCC)

# Plot GO enrichment results
pal = c("#5d8cc2", "#f1bc3a", "#e7673a", "#50a07d", "#0d294e")

go_dotpot <- ggplot(eGo10, aes(enrichment_factor, fct_reorder(factor(Description), enrichment_factor))) + 
  geom_point(aes(size=Count, color=-1*log10(pvalue), shape=ONTOLOGY)) +
  scale_color_gradient(low="green", high ="red") + 
  labs(color=expression(-log[10](p_value)), size="Count", shape="Ontology",
       x="Enrichment Factor", y="GO term", title="GO enrichment") + 
  theme_bw()

print(go_dotpot)

### 2.3.3 Save GO Enrichment Dot Plot

go_dotplot_name <- paste0("004", "_Cis_GO_enrichment-", matches, "_godotplot.pdf")

ggsave(go_dotplot_name, go_dotpot, width = 10, height = 6, units = "in")

#### 2.4 GO Enrichment Bubble Plot

# Create the gene set for GO bubble plot
GOterms = data.frame(category = eGo10$ONTOLOGY,
                     ID = eGo10$ID,
                     term = eGo10$Description, 
                     genes = gsub("/", ",", eGo10$geneID), 
                     adj_pval = eGo10$p.adjust)

genelist <- data.frame(ID = merged_data2$ENTREZID,
                       logFC = merged_data2$log2FoldChange.x)

# Remove rows with "#N/A" in the ID column
genelist <- genelist[!grepl("#N/A", genelist$ID), ]

circ <- circle_dat(GOterms, genelist)  # Generate bubble plot data

go_bubble_plot <- GOBubble(circ, 
                           labels = 3, 
                           table.legend = FALSE, 
                           table.col = TRUE,
                           ID = TRUE, 
                           display = 'single')

print(go_bubble_plot)

### 2.4.2 Save GO Enrichment Bubble Plot

Go_enrichment_gobabble_name <- paste0("005", "_Cis_GO_enrichment-", matches, "_gobabble.pdf")

ggsave(Go_enrichment_gobabble_name,
       go_bubble_plot,
       width = 14, height = 10, 
       units = "in")

#########################################################################

###################################################################
### 3.7 Additional Plots for KEGG Enrichment                      ###
###################################################################

### 3.7.1 KEGG Enrichment Cnet Plot

# Generate KEGG enrichment Cnet plot
kegg_cnet_plot <- cnetplot(kegg_enrichment, 
                           circular = TRUE, 
                           colorEdge = TRUE, 
                           color.params = list(foldChange = GSEA_genelist,
                                               edge = TRUE))

# Print the KEGG cnet plot
print(kegg_cnet_plot)

### 3.7.2 Save KEGG Enrichment Cnet Plot

# Generate file name for KEGG enrichment Cnet plot
kegg_enrich_cnet_name <- paste0("019", "_Cis_KEGG_enrichment-", matches, "_gseaKeggCnet.pdf")

# Save the plot to a file
ggsave(kegg_enrich_cnet_name, 
       kegg_cnet_plot, 
       width = 18, height = 10, 
       units = "in")

###################################################################
### 3.8 KEGG Enrichment Pathview Plot                             ###
###################################################################

### 3.8.1 Pathway Analysis using Pathview

# Generate pathway data for KEGG
pathway_ids <- kegg_enrichment$ID  # Get KEGG pathway IDs

# Loop through each pathway and generate the pathview plot
for (pathway_id in pathway_ids) {
  pathview(gene.data = GSEA_genelist,
           pathway.id = pathway_id,
           species = "osa",
           limit = list(gene = max(abs(GSEA_genelist)), cpd = 1))
}

### 3.8.2 Save KEGG Enrichment Pathview Plots

# Generate file name for saving KEGG pathview results
pathview_file_name <- paste0("020-", matches, "_KEGG_pathview.pdf")

# Save the pathview plot
ggsave(pathview_file_name, 
       width = 14, height = 10, 
       units = "in")

###################################################################
### 4. GSEA Enrichment Analysis                                     ###
###################################################################

# Perform GSEA analysis for genes
GSEA_genelist <- GSEA_input$log2FoldChange
names(GSEA_genelist) <- GSEA_input$ENTREZID

# Run GSEA analysis for GO
gsea_go <- gseGO(geneList = GSEA_genelist,
                 OrgDb = rice,
                 keyType = "ENTREZID",
                 ont = "ALL",
                 pvalueCutoff = 0.05)

### 4.1 GSEA GO Ridgeplot

# Plot GSEA GO Ridgeplot
gsea_go_plot <- ridgeplot(gsea_go)

# Print GSEA GO Ridgeplot
print(gsea_go_plot)

### 4.1.1 Save GSEA GO Ridgeplot

# Generate file name for GSEA GO Ridgeplot
gsea_go_file_name <- paste0("025-", "_GSEA_GO_", matches, "_gseaGoRidge.pdf")

# Save the plot to a file
ggsave(gsea_go_file_name, 
       gsea_go_plot, 
       width = 10, height = 10)

### 4.2 GSEA GO Cnetplot

# Create and plot the GSEA GO Cnetplot
gsea_go_cnet_plot <- cnetplot(gsea_go, 
                              showCategory = 10,
                              foldChange = GSEA_genelist, 
                              circular = TRUE, 
                              colorEdge = TRUE)

# Print the GSEA GO Cnetplot
print(gsea_go_cnet_plot)

### 4.2.1 Save GSEA GO Cnet Plot

# Save GSEA GO Cnet plot to a file
gsea_go_cnet_file_name <- paste0("026-", "_GSEA_GO_", matches, "_gseaGOcnet.pdf")

# Save the plot to a file
ggsave(gsea_go_cnet_file_name, 
       gsea_go_cnet_plot, 
       width = 14, height = 10)

### 4.3 GSEA GO Circle Plot

# Create GSEA GO Circle plot data
gsea_go_circ_data <- circle_dat(gsea_goterms, GSEA_babble_data)

# Generate and print the GSEA GO Circle plot
gsea_go_circle_plot <- GOCircle(gsea_go_circ_data, 
                                rad1 = 2, 
                                rad2 = 3, 
                                label.size = 3, 
                                label.fontface = 'bold', 
                                nsub = 10, 
                                zsc.col = c("red", "#BD378680", "#FCFDBFFF"), 
                                lfc.col = c("#7301A880", "#FDC92680"))

# Print the GSEA GO Circle plot
print(gsea_go_circle_plot)

### 4.3.1 Save GSEA GO Circle Plot

# Save the GSEA GO Circle plot to a file
gsea_go_circle_file_name <- paste0("027-", "_GSEA_GO_", matches, "_gseaGoCircle.pdf")

# Save the plot
ggsave(gsea_go_circle_file_name,
       gsea_go_circle_plot,
       width = 14, height = 10, 
       units = "in")

###################################################################
### 5. KEGG GSEA Enrichment Analysis                                ###
###################################################################

# Perform KEGG GSEA analysis
kegg_gsea <- gseKEGG(GSEA_genelist,
                     organism = "osa",
                     pvalueCutoff = 0.05)

### 5.1 KEGG GSEA Ridgeplot

# Plot KEGG GSEA Ridgeplot
kegg_gsea_plot <- ridgeplot(kegg_gsea)

# Print the KEGG GSEA Ridgeplot
print(kegg_gsea_plot)

### 5.1.1 Save KEGG GSEA Ridgeplot

# Generate file name for KEGG GSEA Ridgeplot
gsea_kegg_file_name <- paste0("030-", "_GSEA_KEGG_", matches, "_gseaKeggRidge.pdf")

# Save the plot to a file
ggsave(gsea_kegg_file_name, 
       gsea_kegg_plot, 
       width = 10, height = 10)

### 5.2 Save KEGG GSEA Results

# Generate and save the KEGG GSEA results to a file
gsea_kegg_file <- paste0("031-", "_GSEA_KEGG_", matches, "_gseaKeggSingle.pdf")

# Save the plot to a file
ggsave(gsea_kegg_file, 
       gsea_kegg_plot, 
       width = 15, height = 7)

###################################################################
### 6. Additional KEGG Enrichment Plots                            ###
###################################################################

# For additional KEGG enrichment plots like bubble plots, heatmaps, and chord plots, you can continue following the pattern above. 

# Example of KEGG enrichment chord plot:
kegg_chord_plot <- GOChord(kegg_cnet_data, 
                           space = 0.02,  # Adjust space between chords
                           gene.order = 'logFC',  # Sort genes based on logFC
                           gene.size = 3,  # Gene size
                           nlfc = 1,  # Use logFC values
                           border.size = NULL, 
                           lfc.col = npg_colors)

# Print and save the KEGG enrichment chord plot
print(kegg_chord_plot)

# Save the plot
kegg_chord_name <- paste0("032", "_GSEA_KEGG_", matches, "_gseaKeggCnet.pdf")

ggsave(kegg_chord_name, 
       kegg_chord_plot, 
       width = 12, height = 10)

###################################################################

```
## 8) WGCNA Analysis for Key lncRNAs

- Use **WGCNA** to perform co-expression network analysis. This identifies modules of lncRNAs that are co-expressed with mRNAs and may help discover key regulatory lncRNAs.
```R
################## 1. Load R Packages and Data ###########################

# Clear environment variables to avoid previous data affecting current analysis
rm(list=ls())

# Set stringsAsFactors to FALSE to prevent automatic conversion of strings to factors
options(stringsAsFactors = FALSE)

# Load required libraries
library(WGCNA)
library(ggplot2)
library(pheatmap)
library(tidyverse)
library(DESeq2)

## 1.2 Set working directory and read gene expression data matrix
setwd("./working_directory/")  # Example directory

# Read the gene expression data file (tab-separated)
filename <- "./example_folder/merged_gene_expression_data_matrix.txt"  # Example file path

# Read the data
gene_info <- read.csv(filename, header=TRUE)

# Set column name for gene IDs
names(gene_info) <- c("MSU")

# Extract file name without extension
filename_no_ext <- sub("\\.txt$", "", basename(filename))

# Define regular expression pattern to remove prefix "02-" and suffix "_sig_padj"
pattern <- "^02-(.*)_sig_padj$"
matches <- sub(pattern, "\\1", filename_no_ext)

matches <- "two_inter_21_DEL_target_gene_id"

############################################################ 

### 1.4 Convert MSU ID to SYMBOL and RAP ID
# Use RiceIDConvert function to convert MSU ID to SYMBOL ID and RAP ID

msu_symbol_id <- RiceIDConvert(myID = gene_info, 
                               fromType = 'MSU', 
                               toType = 'SYMBOL')

# Ensure column names are correct
colnames(msu_symbol_id) <- c("MSU", "SYMBOL")

### 1.5 Load orgDB Database

rice <- loadDb("./example_folder/ricedb")  # Example database path

### 1.6 Extract SYMBOL ID
# Assign the converted SYMBOL ID to rapgene variable for further enrichment analysis
symbol_id <- as.character(msu_symbol_id[,2])

symbol_id <- unique(symbol_id)
length(symbol_id)

# Check available keytypes in rice database
keytypes(rice)

### 1.7 Convert SYMBOL to ENTREZID
# Use AnnotationDbi::select function to select required columns from rice database
gene_annotations <- AnnotationDbi::select(rice, 
                                          keys = symbol_id, 
                                          columns = c("ENTREZID", 
                                                      "SYMBOL", 
                                                      "GENENAME", 
                                                      "GOALL", 
                                                      "REFSEQ",
                                                      "ONTOLOGYALL",
                                                      "ALIAS"), 
                                          keytype = "SYMBOL")
dim(gene_annotations)

#### 1.8 Merge DataFrames
# Merge data frames based on SYMBOL
merged_data <- merge(gene_annotations,
                     msu_symbol_id, by = "SYMBOL")
dim(merged_data)

# Remove duplicate entries
merged_data <- unique(merged_data)
dim(merged_data)

# Remove rows with "#N/A" or NA values
merged_data <- merged_data[!apply(merged_data, 1, function(row) any(row == "#N/A")), ]
merged_data <- na.omit(merged_data)

dim(merged_data)

# Continue merging with gene_info
merged_data2 <- merge(merged_data,
                      gene_info, by = "MSU")
dim(merged_data2)

## 1.8 Save the merged DataFrame
# Generate file name for merged data
merged_data_filename <- paste0("001-", matches, "_annotation.tsv")

write.table(merged_data2, 
            file =  merged_data_filename,
            sep = "\t", 
            row.names = FALSE, 
            col.names = TRUE, 
            quote = FALSE)

### 1.9 Extract Gene Set and Genes
### Gene set used for GSEA and with log2FC, genes for KEGG/GO enrichment analysis input

# Extract ENTREZID and log2FoldChange columns, and remove duplicates
gene <- unique(merged_data[, c("ENTREZID")])

str(gene)

#############################################################

###################################################################
#
#######                2. GO Enrichment Analysis              ###########
#
###################################################################

### 2.1 Perform GO Enrichment Analysis

Go_enrichment <- enrichGO(gene,
                          OrgDb = rice, # Replace with the appropriate organism database
                          keyType = "ENTREZID", # Gene ID type
                          ont = "ALL", # Or "BP", "MF", "CC"
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.12, # Default 0.05
                          qvalueCutoff = 0.2)  # Default 0.2

# Generate dot plot for GO enrichment
Go_enrichment_dotplot <- dotplot(Go_enrichment, 
                                 showCategory = 10,
                                 split="ONTOLOGY") +
  facet_grid(ONTOLOGY~.,scale='free')

print(Go_enrichment_dotplot)

### 2.2 Save Enrichment Analysis Results
# Convert enrichment results to data frame and save to file
Go_enrichment_dataframe <- as.data.frame(Go_enrichment)

# Generate output file name
Go_enrichment_name <- paste0("002-", "_Cis_GO_enrichment-", matches, ".tsv")

write.table(Go_enrichment_dataframe, 
            file = Go_enrichment_name, 
            sep = "\t", 
            row.names = FALSE, 
            col.names = TRUE, 
            quote = FALSE)

#### 2.3 Plot GO Enrichment Dot Plot

# Calculate enrichment factor
Go_enrichment_tsv = read.table(Go_enrichment_name,
                               header=TRUE,
                               sep="\t",
                               quote = "",
                               fill = TRUE) 

Go_enrichment_tsv <- separate(data=Go_enrichment_tsv, 
                              col=GeneRatio,
                              into = c("GR1", "GR2"), 
                              sep = "/")  # Split GeneRatio into 2 columns (GR1, GR2)

Go_enrichment_tsv <- separate(data=Go_enrichment_tsv, 
                              col=BgRatio, 
                              into = c("BR1", "BR2"), 
                              sep = "/")  # Split BgRatio into 2 columns (BR1, BR2)

Go_enrichment_tsv <- mutate(Go_enrichment_tsv,
                            enrichment_factor = (as.numeric(GR1)/as.numeric(GR2))/(as.numeric(BR1)/as.numeric(BR2))) # Calculate Enrichment Factor 

### Save Enrichment Factor
go_enrichment_factor_filename <- paste0("003-", "_Cis_GO_enrichment-", matches, "_enrichmentFactor.tsv") 

write.table(Go_enrichment_tsv, 
            file = go_enrichment_factor_filename,
            sep = "\t", 
            col.names = TRUE,
            quote = FALSE,
            row.names = FALSE)

### 2.3.2 Plot GO Enrichment (Top 10 terms)

eGoBP <- Go_enrichment_tsv %>% 
  filter(ONTOLOGY=="BP") %>%
  filter(row_number() >= 1, row_number() <= 10)

eGoCC <- Go_enrichment_tsv %>% 
  filter(ONTOLOGY=="CC") %>%
  filter(row_number() >= 1, row_number() <= 10)

eGoMF <- Go_enrichment_tsv %>% 
  filter(ONTOLOGY=="MF") %>%
  filter(row_number() >= 1, row_number() <= 10)

eGo10 <- rbind(eGoBP, eGoMF, eGoCC)

# Plot GO enrichment results
pal = c("#5d8cc2", "#f1bc3a", "#e7673a", "#50a07d", "#0d294e")

go_dotpot <- ggplot(eGo10, aes(enrichment_factor, fct_reorder(factor(Description), enrichment_factor))) + 
  geom_point(aes(size=Count, color=-1*log10(pvalue), shape=ONTOLOGY)) +
  scale_color_gradient(low="green", high ="red") + 
  labs(color=expression(-log[10](p_value)), size="Count", shape="Ontology",
       x="Enrichment Factor", y="GO term", title="GO enrichment") + 
  theme_bw()

print(go_dotpot)

### 2.3.3 Save GO Enrichment Dot Plot

go_dotplot_name <- paste0("004", "_Cis_GO_enrichment-", matches, "_godotplot.pdf")

ggsave(go_dotpot_name, go_dotpot, width = 10, height = 6, units = "in")

#### 2.4 GO Enrichment Bubble Plot

# Create the gene set for GO bubble plot
GOterms = data.frame(category = eGo10$ONTOLOGY,
                     ID = eGo10$ID,
                     term = eGo10$Description, 
                     genes = gsub("/", ",", eGo10$geneID), 
                     adj_pval = eGo10$p.adjust)

genelist <- data.frame(ID = merged_data2$ENTREZID,
                       logFC = merged_data2$log2FoldChange.x)

# Remove rows with "#N/A" in the ID column
genelist <- genelist[!grepl("#N/A", genelist$ID), ]

circ <- circle_dat(GOterms, genelist)  # Generate bubble plot data

go_bubble_plot <- GOBubble(circ, 
                           labels = 3, 
                           table.legend = FALSE, 
                           table.col = TRUE,
                           ID = TRUE, 
                           display = 'single')

print(go_bubble_plot)

### 2.4.2 Save GO Enrichment Bubble Plot

Go_enrichment_gobabble_name <- paste0("005", "_Cis_GO_enrichment-", matches, "_gobabble.pdf")

ggsave(Go_enrichment_gobabble_name,
       go_bubble_plot,
       width = 14, height = 10, 
       units = "in")

#########################################################################

###################################################################
### 3. KEGG Enrichment Analysis                                    ###
###################################################################

### 3.1 Perform KEGG Enrichment Analysis

kegg_enrichment <- enrichKEGG(gene,
                              organism = "osa",  # Use the correct organism (example: "osa" for Oryza sativa)
                              pvalueCutoff = 0.1)

# Check the KEGG enrichment results
kegg_enrichment$ID  # List of KEGG pathways enriched

# 3.1.2 Clean up the KEGG enrichment descriptions
kegg_enrichment@result$Description <- gsub(" - Oryza sativa japonica \\(Japanese rice\\) \\(RefSeq\\)", "", kegg_enrichment@result$Description)

# Convert the KEGG enrichment results to a data frame
kegg_enrichment_dataframe <- as.data.frame(kegg_enrichment)

# Generate the new file name for KEGG enrichment results
kegg_enrichment_name <- paste0("014-", "_KEGG_enrichment-", matches, ".tsv")

# Write the KEGG enrichment results to a file
write.table(kegg_enrichment_dataframe, 
            file = kegg_enrichment_name, 
            sep = "\t", 
            col.names = TRUE, 
            row.names = FALSE, 
            quote = FALSE)

### 3.2 Plot KEGG Enrichment Results (Dotplot)

kegg_enrich_dotplot <- dotplot(kegg_enrichment, 
                               showCategory = 15)

print(kegg_enrich_dotplot)

### 3.2.2 Calculate KEGG Enrichment Factor
kegg_enrich_tsv = read.table(kegg_enrichment_name,
                             header=TRUE,
                             sep="\t",
                             quote = "",
                             fill = TRUE)  # Read KEGG enrichment results

kegg_enrich_tsv <- separate(data=kegg_enrich_tsv, 
                             col=GeneRatio,
                             into = c("GR1", "GR2"), 
                             sep = "/")  # Split GeneRatio into two columns

kegg_enrich_tsv <- separate(data=kegg_enrich_tsv, 
                             col=BgRatio, 
                             into = c("BR1", "BR2"), 
                             sep = "/")  # Split BgRatio into two columns

kegg_enrich_tsv <- mutate(kegg_enrich_tsv,
                          enrichment_factor = (as.numeric(GR1)/as.numeric(GR2))/(as.numeric(BR1)/as.numeric(BR2))) # Calculate Enrichment Factor 

### Save the Enrichment Factor Results
kegg_enrich_factor_filename <- paste0("015-", "_Cis_KEGG_enrichment-", matches, "_enrichmentFactor.tsv")

write.table(kegg_enrich_tsv, 
            file = kegg_enrich_factor_filename,
            sep = "\t", 
            col.names = TRUE,
            quote = FALSE,
            row.names = FALSE)

###################################################################
### 3.3 KEGG Enrichment Circle Plot

kegg_goterms <- data.frame(
  category = kegg_enrich_tsv$category,
  ID = kegg_enrich_tsv$ID,
  term = kegg_enrich_tsv$Description, 
  genes = gsub("/", ",", kegg_enrich_tsv$core_enrichment), 
  adj_pval = kegg_enrich_tsv$p.adjust
)

# Create bubble plot data
kegg_circ <- circle_dat(kegg_goterms, genelist)  # Create bubble plot data

kegg_circular_plot <- GOCircle(kegg_circ, 
                               rad1 = 2,  # Inner radius
                               rad2 = 3,  # Outer radius
                               label.size = 3,  # Label size
                               label.fontface = 'bold',  # Label font
                               nsub = 10,  # Number of terms to display
                               zsc.col = c('red', 'yellow'),  # Z-score color
                               lfc.col = c('red', '#FDC92680'))  # LogFC color

### 3.3.1 Save the KEGG Enrichment Circle Plot
kegg_circle_plot_name <- paste0("017", "_Cis_KEGG_enrichment-", matches, "_keggCirclePlot.pdf")

ggsave(kegg_circle_plot_name, 
       kegg_circular_plot, 
       width = 14, height = 10, 
       units = "in")

###################################################################
### 3.4 KEGG Enrichment Heatmap

# Generate heatmap data for KEGG enrichment
kegg_heat_data <- chord_dat(kegg_circ, genelist, kegg_goterms$term)  # Use previously created data

kegg_heat_plot <- GOHeat(kegg_heat_data, 
                         nlfc = 1,  # Use logFC if available
                         fill.col = c('red', 'grey60', 'green'))  # Custom colors

print(kegg_heat_plot)

### 3.4.1 Save the KEGG Enrichment Heatmap Plot
kegg_heat_plot_name <- paste0("018", "_Cis_KEGG_enrichment-", matches, "_keggHeatmap.pdf")

ggsave(kegg_heat_plot_name,
       kegg_heat_plot,
       width = 16, height = 10, 
       units = "in")

###################################################################
### 3.5 KEGG Enrichment Cnet Plot

kegg_cnet_plot <- cnetplot(kegg_enrichment, 
                           circular = TRUE, 
                           colorEdge = TRUE, 
                           color.params = list(foldChange = GSEA_genelist,
                                               edge = TRUE))

print(kegg_cnet_plot)

### 3.5.1 Save KEGG Enrichment Cnet Plot
kegg_enrich_cnet_name <- paste0("019", "_Cis_KEGG_enrichment-", matches, "_gseaKeggCnet.pdf")

ggsave(kegg_enrich_cnet_name, 
       kegg_cnet_plot, 
       width = 18, height = 10, 
       units = "in")

###################################################################
### 3.6 KEGG Enrichment Pathway View

# Generate pathway data for KEGG
# Example code for generating pathways using pathview

pathway_ids <- kegg_enrichment$ID  # Get KEGG pathway IDs

for (pathway_id in pathway_ids) {
  pathview(gene.data = GSEA_genelist,
           pathway.id = pathway_id,
           species = "osa",
           limit = list(gene = max(abs(GSEA_genelist)), cpd = 1))
}

### 3.6.1 Save KEGG Enrichment Pathview
pathview_file_name <- paste0("020-", matches, "_KEGG_pathview.pdf")

ggsave(pathview_file_name, 
       width = 14, height = 10, 
       units = "in")

###################################################################
### 4. GSEA Enrichment Analysis                                     ###
###################################################################

# Perform GSEA analysis for genes
GSEA_genelist <- GSEA_input$log2FoldChange
names(GSEA_genelist) <- GSEA_input$ENTREZID

# Run GSEA analysis (example with GSEA_GO)
gsea_go <- gseGO(geneList = GSEA_genelist,
                 OrgDb = rice,
                 keyType = "ENTREZID",
                 ont = "ALL",
                 pvalueCutoff = 0.05)

### 4.1 GSEA GO Ridgeplot
gsea_go_plot <- ridgeplot(gsea_go)

print(gsea_go_plot)

### 4.1.1 Save GSEA GO Ridgeplot
gsea_go_file_name <- paste0("025-", "_GSEA_GO_", matches, "_gseaGoRidge.pdf")

ggsave(gsea_go_file_name, 
       gsea_go_plot, 
       width = 10, height = 10)

### 4.2 GSEA GO Cnetplot

gsea_go_cnet_plot <- cnetplot(gsea_go, 
                              showCategory = 10,
                              foldChange = GSEA_genelist, 
                              circular = TRUE, 
                              colorEdge = TRUE)

print(gsea_go_cnet_plot)

### 4.2.1 Save GSEA GO Cnet Plot
gsea_go_cnet_file_name <- paste0("026-", "_GSEA_GO_", matches, "_gseaGOcnet.pdf")

ggsave(gsea_go_cnet_file_name, 
       gsea_go_cnet_plot, 
       width = 14, height = 10)

### 4.3 GSEA GO Circle Plot

gsea_go_circ_data <- circle_dat(gsea_goterms, GSEA_babble_data)

gsea_go_circle_plot <- GOCircle(gsea_go_circ_data, 
                                rad1 = 2, 
                                rad2 = 3, 
                                label.size = 3, 
                                label.fontface = 'bold', 
                                nsub = 10, 
                                zsc.col = c("red", "#BD378680", "#FCFDBFFF"), 
                                lfc.col = c("#7301A880", "#FDC92680"))

### 4.3.1 Save GSEA GO Circle Plot
gsea_go_circle_file_name <- paste0("027-", "_GSEA_GO_", matches, "_gseaGoCircle.pdf")

ggsave(gsea_go_circle_file_name,
       gsea_go_circle_plot,
       width = 14, height = 10, 
       units = "in")

###################################################################
### 5. KEGG GSEA Enrichment Analysis                                ###
###################################################################

kegg_gsea <- gseKEGG(GSEA_genelist,
                     organism = "osa",
                     pvalueCutoff = 0.05)

### 5.1 KEGG GSEA Ridgeplot
kegg_gsea_plot <- ridgeplot(kegg_gsea)

print(kegg_gsea_plot)

### 5.1.1 Save KEGG GSEA Ridgeplot
gsea_kegg_file_name <- paste0("030-", "_GSEA_KEGG_", matches, "_gseaKeggRidge.pdf")

ggsave(gsea_kegg_file_name, 
       gsea_kegg_plot, 
       width = 10, height = 10)

### 5.2 Save KEGG GSEA Results
gsea_kegg_file <- paste0("031-", "_GSEA_KEGG_", matches, "_gseaKeggSingle.pdf")

ggsave(gsea_kegg_file, 
       gsea_kegg_plot, 
       width = 15, height = 7)

###################################################################


```

## 9) GSEA Functional Enrichment for Key lncRNAs

- Perform **GSEA** (Gene Set Enrichment Analysis) for key lncRNAs identified through WGCNA or DESeq2 to identify enriched pathways and functional categories.
```bash

# Load required libraries
library(enrichplot)
library(future.apply)
library(clusterProfiler)
library(AnnotationHub)
library(AnnotationDbi)
library(riceidconverter)
library(ggplot2)
library(dplyr)

# Set working directory
setwd("./working_directory/")  # Example directory

# Read data
exprSet <- read.table("./example_folder/2024-3-27_normal_count_lncrna_mrna.tsv", 
                      header = TRUE, sep = "\t", check.names = FALSE, 
                      row.names = 1)

# Modified batch_cor function for correlation analysis
batch_cor <- function(target_gene) {
  if (!target_gene %in% rownames(exprSet)) {
    stop("Target gene not found in exprSet row names.")
  }
  
  y <- as.numeric(exprSet[target_gene, ])
  
  # Perform correlation analysis for all genes using parallel processing
  results <- future_lapply(rownames(exprSet), function(gene) {
    x <- as.numeric(exprSet[gene, ])
    if (gene != target_gene) {
      dd <- tryCatch({
        cor.test(x, y, method = "spearman")
      }, error = function(e) {
        return(list(estimate = NA, p.value = NA))
      })
      data.frame(gene = target_gene, correlated_gene = gene, cor = dd$estimate, p.value = dd$p.value)
    }
  })
  
  do.call(rbind, results)
}

# Set up parallel processing
plan(multisession)

# Perform correlation analysis for a specific lncRNA
target_lncrna <- "LncRNA.9562.1"  # Example lncRNA ID
# Perform correlation analysis and measure time taken
system.time(correlation_results <- batch_cor(target_lncrna))

# Sort the correlation results by absolute correlation value (descending)
correlation_results <- correlation_results %>%
  arrange(desc(abs(cor)))

# Save results to a file
write.csv(correlation_results, file = paste0(target_lncrna, "_correlation_results.csv"), row.names = FALSE)

# Print the top 10 correlated genes
print(head(correlation_results, 10))

# Optionally, continue with GSEA analysis using the top correlated genes
# Select top 1000 correlated genes for further analysis
top_genes <- head(correlation_results, 1000)

# Prepare gene list for GSEA
gene_df <- data.frame(logFC = top_genes$cor, SYMBOL = top_genes$correlated_gene)

# Load annotation database
rice <- loadDb("./example_folder/ricedb")  # Example database path

# Convert gene IDs from SYMBOL to ENTREZID
gene <- bitr(gene_df$SYMBOL, fromType="SYMBOL", toType="ENTREZID", OrgDb=rice)

# Merge gene conversion results with correlation data
gene_df <- merge(gene_df, gene, by="SYMBOL")

# Prepare geneList for GSEA
geneList <- gene_df$logFC
names(geneList) <- gene_df$ENTREZID
geneList <- sort(geneList, decreasing = TRUE)

# Perform GO GSEA
go_gsea <- gseGO(geneList = geneList,
                 OrgDb = rice,
                 keyType = "ENTREZID",
                 ont = "ALL",
                 minGSSize = 10,
                 maxGSSize = 500,
                 pvalueCutoff = 0.05)

# Plot GO GSEA results
ridge_go_plot <- ridgeplot(go_gsea, label_format = 100)
print(ridge_go_plot)

# Save the GO GSEA plot
go_gsea_name <- paste0(target_lncrna, "_GO_GSEA_results_ridgeplot.pdf")
ggsave(go_gsea_name, ridge_go_plot, width = 10, height = 12, dpi = 300)

# Perform KEGG GSEA
kegg_gsea <- gseKEGG(geneList,
                     organism = "osa",
                     pvalueCutoff = 0.05,
                     pAdjustMethod = "BH")

# Plot KEGG GSEA results
kegg_ridge_plot <- ridgeplot(kegg_gsea)
print(kegg_ridge_plot)

# Save the KEGG GSEA plot
kegg_gsea_name <- paste0(target_lncrna, "_KEGG_GSEA_results_ridgeplot.pdf")
ggsave(kegg_gsea_name, kegg_ridge_plot, width = 10, height = 12, dpi = 300)

# Save GSEA results
write.csv(as.data.frame(go_gsea), file = paste0(target_lncrna, "_GO_GSEA_results.csv"), row.names = FALSE)
write.csv(as.data.frame(kegg_gsea), file = paste0(target_lncrna, "_KEGG_GSEA_results.csv"), row.names = FALSE)

# Check available SYMBOL keys in rice database
keys(rice, keytype = "SYMBOL")

```

---

### Final Output

The functional analysis pipeline identifies key lncRNAs that regulate biological processes and stress responses in rice. Through **differential expression**, **target gene prediction**, and **enrichment analysis**, we gain insight into the roles of lncRNAs in disease resistance, hormone signaling, and metabolic processes.

This phase of the pipeline helps refine the understanding of **lncRNA function** and provides targets for further experimental validation.

