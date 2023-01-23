## This script aims to take counts of peptides and do analysis of them.

library(tidyverse)  #Needed for data wrangling
library(DESeq2)  #Use for differential expression
library(bioseq) # Used to translate sequence

setwd(dir = "C:/Users/worms/Dropbox/PhD/PhD-Scripts/CyclicPeptidePipeline")
source("library_installation.R")  #Install and/or load needed libraries
source("function.R")  #Functions used in the script

# This folder should contain the CSV files with a column for sequence and a column for counts
directory <- "C:/Users/worms/NGS Data/2022.06.07_drift_seq/90-666155004b/00_fastq/NNK/NNK3/counts_csv"

# Get the names of the counts .csv files

file_list <- list.files(directory) 
file_list <- file_list[grepl(".csv", file_list)]

# Trim the file names to use as variable names. 

run_names <- file_list %>%
  gsub("Cytoplasmic-NNK-", "", .) %>%
  gsub("_001_peptide3_count.csv", "", .) %>%
  tolower() %>%
  gsub("-", "_", .) %>%
  gsub("gen_", "gen", .)

# Import the sets
merged_set <- createCountSet(directory, file_list, run_names)

# Going to split the set for induced and repressed conditions. Drop the lines containing only zeros

induced_set <- merged_set %>% 
  select(!contains("glu")) %>%
  filter(rowSums(.[-c(1,2)]) != 0)

# Don't need medium info anymore
colnames(induced_set) <- gsub("ara_|lb_",
                              "",
                              colnames(induced_set))

repressed_set <- merged_set %>%
  select(!contains("ara")) %>%
  filter(rowSums(.[-c(1,2)]) != 0)

colnames(repressed_set) <- gsub("ara_|lb_|glu_",
                              "",
                              colnames(repressed_set))

## Compute sum of counts for each condition and then compute ratios of each read as a percentage of total
## compute an enrichment ratio based on the final ratio over the initial one, and removes NAs caused by divide by zero error

induced_set <- induced_set %>% computeRatios()
repressed_set <- repressed_set %>% computeRatios()

induced_set$enrichment_ratio %>%
  summary()

repressed_set$enrichment_ratio %>%
  summary()
# 
# # This is nice, it hints that some were highly selected. Let's look at them. 
# 
# induced_set %>%
#   arrange(desc(enrichment_ratio)) %>% 
#   head(n = 20)
# 
# # Really cool, look like few stop codons (good) and a common CRSX motif
# 
# induced_set_enriched <- induced_set %>%
#   filter(enrichment_ratio > 4)
# 
# induced_set_depleted <- induced_set %>% 
#   filter(enrichment_ratio < 0.25)
# 
# # Let's see about the subset wdepleted
# 
# 
# induced_set %>%
#   arrange(enrichment_ratio) %>% 
#   head(n = 20)
# 
# #Really strange, they all have the same number of reads in all conditions!
# 
# induced_set %>%
#   arrange(enrichment_ratio) %>% 
#   filter(enrichment_ratio > 0.62)
#   head(n = 20)
# 
# ## Few depleted peptides: we can't detect toxic peptides?
#   
#   
#   repressed_set$enrichment_ratio %>%
#     summary()
#   
#   # This is nice, it hints that some were highly selected. Let's look at them. 
#   
#   repressed_set %>%
#     arrange(desc(enrichment_ratio)) %>% 
#     head(n = 20)
#   
#   # Really cool, look like few stop codons (good) and a common CRS* motif
#   
#   repressed_set_enriched <- repressed_set %>%
#     filter(enrichment_ratio > 4)
#   
#   repressed_set_depleted <- repressed_set %>% 
#     filter(enrichment_ratio < 0.25)
#   
#   # Let's see about the subset wdepleted
#   
#   
#   repressed_set %>%
#     arrange(enrichment_ratio) %>% 
#     head(n = 20)
#   
#   #Really strange, they all have the same number of reads in all conditions!
#   
#   repressed_set %>%
#     arrange(enrichment_ratio) %>% 
#     filter(enrichment_ratio > 0.62)
