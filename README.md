**Understanding Differential Gene Expression Analysis with DESeq2**

In this guide, we will explore how DESeq2 identifies differentially expressed genes and how to use this package in R to analyze gene expression differences between two conditions.

We will follow the DESeq2 standard workflow vignette and work with a publicly available dataset, which is also accessible as a Bioconductor package.

The dataset consists of RNA sequencing data from four airway smooth muscle cell lines, each with one treated and one untreated sample. The treatment involves dexamethasone, a medication commonly used by asthma patients to reduce airway inflammation. Our goal is to analyze transcriptional changes induced by dexamethasone treatment.

To perform this analysis, we will use three key R packages:

DESeq2 – for differential gene expression analysis

tidyverse – for data manipulation

airway – to access the dataset from Bioconductor

By the end of this guide, you will understand how to process RNA-seq data, perform differential expression analysis, and interpret the results using DESeq2.



# Script-for-Differential-Gene-Expression-Analysis-using-DESeq2-N

This script uses the DESeq2 package in R to perform differential gene expression (DGE) analysis. The goal is to identify genes that are significantly differentially expressed between two conditions (e.g., treated vs untreated) by analyzing count data from RNA sequencing (RNA-seq) experiments.

**Step 1: Setting Up and Loading Required Libraries**

library(DESeq2)
library(tidyverse)
library(airway)

**DESeq2:** The core package for performing differential gene expression analysis.

**tidyverse:** A collection of R packages for data manipulation and visualization.

**airway:** A dataset package often used for tutorials in RNA-seq analysis.


**Step 2: Preparing the Count Data**
**Loading the Count Data**

counts_data <- read.csv('counts_data.csv')
head(counts_data)


Reads the gene expression count data from a CSV file (counts_data.csv).

This file should have:

  Rows: Genes (Gene IDs)

  Columns: Samples (Sample IDs)

  Each cell represents the raw RNA-seq read counts for a given gene in a particular sample.

  **Loading the Sample Information**

  colData <- read.csv('sample_info.csv')

Reads the metadata file (sample_info.csv), which contains information about each sample.

This file typically includes:

  Sample IDs (matching the column names in counts_data).

  Experimental conditions (e.g., treatment vs. control).

  Other relevant metadata (e.g., batch information, sequencing run details).

  **Ensuring Count Data Matches Sample Information**

all(colnames(counts_data) %in% rownames(colData))
all(colnames(counts_data) == rownames(colData))

The column names in counts_data should match the row names in colData.

The second check ensures that the order of samples in both datasets is identical.

**Step 3: Constructing a DESeqDataSet Object**

dds <- DESeqDataSetFromMatrix(countData = counts_data,
                              colData = colData,
                              design = ~ dexamethasone)


Creates a DESeqDataSet object, which is required for DESeq2.

Arguments:

  countData: The raw count matrix.

  colData: The metadata table.

  design: The experimental design formula, specifying the variable (e.g., "dexamethasone") to compare conditions.

  **Pre-filtering: Removing Low-Count Genes**


keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

Filters out genes with low read counts (less than 10 reads across all samples).

This improves statistical power and reduces noise.

**Setting Factor Levels for Comparison**

dds$dexamethasone <- relevel(dds$dexamethasone, ref = "untreated")

Defines "untreated" as the reference group.

This ensures that differential expression is computed relative to the untreated samples.

**Step 4: Running DESeq Analysis**

dds <- DESeq(dds)
res <- results(dds)


Runs the differential expression analysis using DESeq().

Extracts the results table, which includes:

  log2FoldChange (expression difference between conditions).

  p-value (statistical significance of the difference).

  adjusted p-value (corrected for multiple testing).
  

**Step 5: Exploring the Results**

summary(res)

Provides a summary of the differential expression results.

Reports:

  The number of significantly upregulated and downregulated genes.

  The total number of genes tested.

  **Applying a More Stringent Significance Threshold**

  res0.01 <- results(dds, alpha = 0.01)
summary(res0.01)

Extracts results using a more stringent threshold (alpha = 0.01 instead of the default 0.05).

This reduces false positives by considering only highly significant genes.

**Retrieving Names of Available Comparisons**

resultsNames(dds)

Lists the available conditions in the dataset (e.g., "treated_4hrs", "treated_8hrs", "untreated").

**Performing a Specific Comparison**

results(dds, contrast = c("dexamethasone", "treated_4hrs", "untreated"))

Performs a direct comparison between "treated_4hrs" and "untreated".

Returns the differentially expressed genes between these conditions.

**Step 6: Visualizing the Results**

plotMA(res)

Generates an MA plot, which:

  Plots log2 fold changes (y-axis) against mean expression levels (x-axis).

  Highlights significantly differentially expressed genes in red.

  


