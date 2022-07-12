## This script aims to take likes of peptides and do DESeq analysis of them.

library(tidyverse)
library(DESeq2)

directory = "D:/2022.06.07_drift_seq/90-666155004b/00_fastq/NNK/NNK3/aa_counts_csv"

# Get the count data

fileList <- list.files(directory)

varNames =fileList %>%
  gsub("Cytoplasmic-","",.) %>%
  gsub("_001_peptide3_count_aa.csv","",.)

for(i in seq_along(fileList)){
  
  counts = read.csv(file.path(directory, fileList[i]),
           header = FALSE)
  
  colnames(counts) = c("seq", "count")
  
  assign(varNames[i], 
         counts)

}

## Merge by peptide
mergedSet = get(varNames[1])

for(i in 2:length(fileList)){
  mergedSet = full_join(mergedSet, get(varNames[i]), by = "seq")
}

columnNames = c("seq", varNames)
colnames(mergedSet) = columnNames

#Clean the mergedSet from stop codons

dropThose = mergedSet$seq %>% grepl("\\*" , . ) 
  
mergedSet <- mergedSet[!dropThose,]$seq


dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design= ~ batch + condition)
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
res <- results(dds, name="condition_trt_vs_untrt")
# or to shrink log fold changes association with condition:
res <- lfcShrink(dds, coef="condition_trt_vs_untrt", type="apeglm")