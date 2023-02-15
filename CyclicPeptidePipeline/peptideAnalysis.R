# This script is an ETL aiming to take counts of peptides from the csv of counts produced by the
# python scripts and make them into one .csv file containing a dataframe ready for exploratory
# analysis. That data frame has counts for all conditions as well as ratio of total read per
# conditions and enrichment ratios.

setwd(dir = "C:/Users/worms/Dropbox/PhD/PhD-Scripts/CyclicPeptidePipeline")  # Directory containing the scripts
source("function.R")  #Functions used in the script

# Load necessary libraries
library(tidyverse)  #Needed for data wrangling
library(bioseq)  # Used to translate sequence

# This folder should contain the CSV files with a column for sequence and a column for counts
directory <- "C:/Users/worms/NGS Data/2022.06.07_drift_seq/90-666155004b/00_fastq/NNK/NNK7/counts_csv"

# Get the names of the counts .csv files

file_list <- list.files(directory)
file_list <- file_list[grepl(".csv", file_list)]

# Trim the file names to use as variable names.

run_names <- file_list %>%
  gsub("Cytoplasmic-NNK-", "", .) %>%
  gsub("_001_peptide24_count_aa.csv", "", .) %>%
  tolower() %>%
  gsub("-", "_", .) %>%
  gsub("gen_", "gen", .)

# Import the sets
merged_set <- createCountSet(directory, file_list, run_names)

# Going to split the set for induced and repressed conditions. Drop the lines containing only
# zeros

induced_set <- merged_set %>%
  select(!contains("glu")) %>%
  filter(rowSums(.[-c(1, 2)]) != 0)

# Don't need medium info anymore
colnames(induced_set) <- gsub("ara_|lb_|glu_", "", colnames(induced_set))

repressed_set <- merged_set %>%
  select(!contains("ara")) %>%
  filter(rowSums(.[-c(1, 2)]) != 0)

colnames(repressed_set) <- gsub("ara_|lb_|glu_", "", colnames(repressed_set))

## Compute sum of counts for each condition and then compute ratios of each read as a percentage
## of total compute an enrichment ratio based on the final ratio over the initial one, and removes
## NAs caused by divide by zero error

induced_set <- induced_set %>%
  computeRatios()
repressed_set <- repressed_set %>%
  computeRatios()

induced_set$enrichment_ratio %>%
  summary()

repressed_set$enrichment_ratio %>%
  summary()

## Export the two data sets

## Try to generalize the file names
peptide_types <- file_list[1] %>%
  substr(., regexpr("NN", .), regexpr("NN", .) + 2)

peptide_types <- paste(peptide_types, nchar(merged_set[1, 1])/3 - 1, sep = "")

write.csv(induced_set, file = file.path(dirname(directory), 
                                        paste("induced_set_", 
                                              peptide_types, 
                                              ".csv", sep = "")))

write.csv(repressed_set, file = file.path(dirname(directory), paste("repressed_set_", peptide_types, ".csv",
  sep = "")))

## Make a merged long set and write it

count_set <- rbind(induced_set %>% 
                     mutate(induction = "induced"),
                   repressed_set %>%
                     mutate(induction = "repressed")) %>%
              mutate(induction = as.factor(induction))

## Randomize the order of the data frame so I can just grab a part of it.

count_set <- count_set[sample(nrow(count_set)), ]

write.csv(count_set, file = file.path(dirname(directory), paste("count_set_", peptide_types, ".csv",
                                                           sep = "")))
