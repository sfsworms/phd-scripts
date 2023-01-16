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
  gsub("Cytoplasmic-NNK-", "", .) %>%
  gsub("_001_peptide3_count.csv", "", .) %>%
  gsub("-", "_", .)

for (i in seq_along(fileList)) {

  counts = read.csv(file.path(directory, fileList[i]), header = FALSE)

  colnames(counts) = c("seq", "count")

  assign(varNames[i], counts)

}

rm(counts)

## Merge by peptide to form one big dataset.

mergedSet = get(varNames[1])

for (i in 2:length(fileList)) {
  mergedSet = full_join(mergedSet, get(varNames[i]), by = "seq")
}

columnNames = c("seq", varNames)
colnames(mergedSet) = columnNames

# Replace all NAs by 0

mergedSet <- mergedSet %>%
  mutate_at(c(2:ncol(mergedSet)),
           ~replace_na(.,0))

rm(list = varNames)

# Clean the mergedSet from stop codons. As of jan 23, I'm not doing that

# dropThose = mergedSet$seq %>% grepl('\\*' , . )

# mergedSet <- mergedSet[!dropThose,]

## Compute sum of counts for each condition and then compute ratios

mergedSet <- mergedSet %>% 
  mutate(Gen_1_sum = Gen_1_LB_R1 + Gen_1_LB_R2) %>%
  mutate(Gen_5_Ara_sum = Gen_5_Ara_R1 + Gen_5_Ara_R2) %>%
  mutate(Gen_5_Glu_sum = Gen_5_Glu_R1 + Gen_5_Glu_R2) 

Gen1_read_total <- sum(mergedSet$Gen_1_sum)
Gen5_ara_total <- sum(mergedSet$Gen_5_Ara_sum) 
Gen5_glu_total <- sum(mergedSet$Gen_5_Glu_sum)

mergedSet <- mergedSet %>%
  mutate(Gen_1_ratio = Gen_1_sum / Gen1_read_total) %>%
  mutate(Gen_5_ara_ratio = Gen_5_Ara_sum / Gen5_ara_total) %>%
  mutate(Gen_5_glu_ratio = Gen_5_Glu_sum / Gen5_glu_total)

## Compute enrichment ratio by comparing with Gen 1

mergedSet <- mergedSet %>%
  mutate(enrichment_ara = log2(Gen_5_ara_ratio/Gen_1_ratio)) %>%
  mutate(enrichment_glu = log2(Gen_5_glu_ratio/Gen_1_ratio))

#Have an issue with zero ratio for Gen1, need to add a tiiiiiiiny amount to all I think. Or actually just
# refuse to compute ratio for those? With NAs or something