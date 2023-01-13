## This script aims to take counts of peptides and do DESeq analysis of them.

library(tidyverse) #Needed for data wrangling
library(DESeq2) #Use for differential expression

directory = "D:/Destination/peptide_24_count/count_csv"

# Get the count data

fileList <- list.files(directory)

varNames <- fileList %>%
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
  
mergedSet <- mergedSet[!dropThose,]

#Renames the rows of mergedSet with the gene names

rownames(mergedSet) <- mergedSet[,1]

mergedSet <- mergedSet %>%
  dplyr::select(2:ncol(mergedSet))

mergedSet[is.na(mergedSet)] <- 0

# Get the coldata excel sheet

coldata <- readxl::read_xlsx("D:/Destination/coldata.xlsx") %>%
  as.data.frame()

rownames(coldata) = coldata[,1]

coldata <- coldata %>%
  dplyr::select(condition:type) 

dds <- DESeqDataSetFromMatrix(countData = mergedSet,
                              colData = coldata,
                              design= ~ condition)
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
res <- results(dds, name="condition_gen5_induction_vs_gen1", independentFiltering = FALSE)

x <- (res$pvalue < 0.01) %>% sum()
x/ length((res$pvalue))

alpha <- 0.05 # Threshold on the adjusted p-value
cols <- densCols(res$log2FoldChange, -log10(res$pvalue))
plot(res$log2FoldChange, -log10(res$padj), col=cols, panel.first=grid(),
     main="Volcano plot", xlab="Effect size: log2(fold-change)", ylab="-log10(adjusted p-value)",
     pch=20, cex=0.6)
abline(v=0)
abline(v=c(-1,1), col="brown")
abline(h=-log10(alpha), col="brown")

gn.selected <- abs(res$log2FoldChange) > 2.5 & res$padj < alpha 
text(res$log2FoldChange[gn.selected],
     -log10(res$padj)[gn.selected],
     lab=rownames(res)[gn.selected ], cex=0.4)

# => If FPR is 0.1, I have a low true rate since half of pos would be FP

