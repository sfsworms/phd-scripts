## This script aims to take counts of peptides and do analysis of them.

library(tidyverse)  #Needed for data wrangling
library(DESeq2)  #Use for differential expression

# This folder should contain the CSV files with a column for sequence and a column for counts

directory = choose.dir()

# Get the count data

fileList <- list.files(directory)

#Keep only the data csv

fileList <- fileList[grepl("Gen", fileList)]

varNames <- fileList %>%
  gsub("Cytoplasmic-", "", .) %>%
  gsub("_001_peptide3_count.csv", "", .) 

for (i in seq_along(fileList)) {

  counts = read.csv(file.path(directory, fileList[i]), header = FALSE)

  colnames(counts) = c("seq", "count")

  assign(varNames[i], counts)

}

## Merge by peptide to form one big dataset the exact format of that dataset, with one column per
## condition and one row per gene, with the row named after genes, is necessary for the DESeq
## (Differential Expression) package

DESeqMergedSet = get(varNames[1])

for (i in 2:length(fileList)) {
  DESeqMergedSet = full_join(DESeqMergedSet, get(varNames[i]), by = "seq")
}

columnNames = c("seq", varNames)
colnames(DESeqMergedSet) = columnNames

# Clean the DESeqMergedSet from stop codons. As of jan 23, I'm not doing that

# dropThose = DESeqMergedSet$seq %>% grepl('\\*' , . )

# DESeqMergedSet <- DESeqMergedSet[!dropThose,]

# Renames the rows of DESeqMergedSet with the gene names

rownames(DESeqMergedSet) <- DESeqMergedSet[, 1]

DESeqMergedSet <- DESeqMergedSet %>%
  dplyr::select(2:ncol(DESeqMergedSet))

DESeqMergedSet[is.na(DESeqMergedSet)] <- 0

# Get the coldata excel sheet

coldata <- readxl::read_xlsx(choose.files()) %>%
  as.data.frame()

rownames(coldata) = coldata[, 1]

coldata <- coldata %>%
  dplyr::select(condition:type)

dds <- DESeqDataSetFromMatrix(countData = DESeqMergedSet, colData = coldata, design = ~condition)
dds <- DESeq(dds)
resultsNames(dds)  # lists the coefficients
res <- results(dds, name = "condition_gen5_induction_vs_gen1", independentFiltering = FALSE)

x <- (res$pvalue < 0.01) %>%
  sum()
x/length((res$pvalue))

alpha <- 0.05  # Threshold on the adjusted p-value
cols <- densCols(res$log2FoldChange, -log10(res$pvalue))
plot(res$log2FoldChange, -log10(res$padj), col = cols, panel.first = grid(), main = "Volcano plot",
  xlab = "Effect size: log2(fold-change)", ylab = "-log10(adjusted p-value)", pch = 20, cex = 0.6)
abline(v = 0)
abline(v = c(-1, 1), col = "brown")
abline(h = -log10(alpha), col = "brown")

gn.selected <- abs(res$log2FoldChange) > 2.5 & res$padj < alpha
text(res$log2FoldChange[gn.selected], -log10(res$padj)[gn.selected], lab = rownames(res)[gn.selected],
  cex = 0.4)

# => If FPR is 0.1, I have a low true rate since half of pos would be FP

