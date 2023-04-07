# This script is an ETL aiming to take counts of peptides from the csv of counts produced by the
# python scripts and make them into one .csv file containing a dataframe ready for exploratory
# analysis. That data frame has counts for all conditions as well as ratio of total read per
# conditions and enrichment ratios.

setwd(dir = "C:/Users/worms/phd-scripts/cyclic_peptide_pipeline")  # Directory containing the scripts
source("functions/function_analysis.R")  #Functions used in the script

# Load necessary libraries
library(tidyverse)  #Needed for data wrangling
library(bioseq)  # Used to translate sequence

# This folder should contain the CSV files with a column for sequence and a column for counts
directory <- "C:/Users/worms/NGS Data/2022.06.07_drift_seq/90-666155004b/00_fastq/all_files/peptide_count_csv"

# Get the names of the counts .csv files

file_list <- list.files(directory)
file_list <- file_list[grepl(".csv", file_list)]

# Trim the file names to use as variable names.

run_names <- file_list %>%
  gsub("peptide", "", .) %>%
  gsub("_count.csv", "", .) %>%
  tolower() %>%
  gsub("gen_", "gen", .)

# Split the file names and run names based on the peptide size

run_names_short <- run_names[grepl("12", run_names)]
run_names_long <- run_names[grepl("24", run_names)]

file_list_short <- file_list[grepl("12", file_list)]
file_list_long <- file_list[grepl("24", file_list)]

# Import the sets
merged_set_short <- create_count_set(directory, file_list_short, run_names_short)
merged_set_long <- create_count_set(directory, file_list_long, run_names_long)

# Save the merged set
write.csv2(merged_set_short, file = file.path(directory,"merged_set_short.csv"), row.names = FALSE)
write.csv2(merged_set_long, file = file.path(directory,"merged_set_long.csv"), row.names = FALSE)

# Now going to process the short set. 

# Going to split the set for induced and repressed conditions. Drop the lines containing only
# zeros

induced_set_short <- merged_set_short %>%
  select(!contains("glu")) %>%
  filter(rowSums(.[-c(1, 2)]) != 0)

# Don't need medium info anymore
colnames(induced_set_short) <- gsub("ara_|lb_|glu_", "", colnames(induced_set_short))

repressed_set_short <- merged_set_short%>%
  select(!contains("ara")) %>%
  filter(rowSums(.[-c(1, 2)]) != 0)

colnames(repressed_set_short) <- gsub("ara_|lb_|glu_", "", colnames(repressed_set_short))

## Compute sum of counts for each condition and then compute ratios of each read as a percentage
## of total compute an enrichment ratio based on the final ratio over the initial one, and removes
## NAs caused by divide by zero error

induced_set_short <- induced_set_short %>%
  compute_ratios()

repressed_set_short <- repressed_set_short %>%
  compute_ratios()

## Export the two data sets

peptide_types <- "short"

write.csv2(induced_set_short, 
           file = file.path(dirname(directory), 
                                        paste0("induced_set_", 
                                              peptide_types, 
                                              ".csv")),
           row.names = FALSE)

write.csv2(repressed_set_short, file = file.path(dirname(directory), paste("repressed_set", peptide_types, ".csv",
  sep = "")),row.names = FALSE)

## Make a merged long set and write it

count_set <- rbind(induced_set_short %>% 
                     mutate(induction = "induced"),
                   repressed_set_short %>%
                     mutate(induction = "repressed")) %>%
              mutate(induction = as.factor(induction))

## Randomize the order of the data frame so I can just grab a part of it.

count_set <- count_set[sample(nrow(count_set)), ]

write.csv2(count_set, 
          file = file.path(dirname(directory), paste("count_set_", peptide_types, ".csv",
                                                           sep = "")),
          row.names = FALSE)

# Now the same for the long file

merged_set_long <- read.csv2(file.path(directory,"merged_set_long.csv"))

induced_set_long <- merged_set_long %>%
  select(!contains("glu")) %>%
  filter(rowSums(.[-c(1, 2)]) != 0)

# Don't need medium info anymore
colnames(induced_set_long) <- gsub("ara_|lb_|glu_", "", colnames(induced_set_long))

repressed_set_long <- merged_set_long%>%
  select(!contains("ara")) %>%
  filter(rowSums(.[-c(1, 2)]) != 0)

colnames(repressed_set_long) <- gsub("ara_|lb_|glu_", "", colnames(repressed_set_long))

## Compute sum of counts for each condition and then compute ratios of each read as a percentage
## of total compute an enrichment ratio based on the final ratio over the initial one, and removes
## NAs caused by divide by zero error

induced_set_long <- induced_set_long %>%
  compute_ratios()

repressed_set_long <- repressed_set_long %>%
  compute_ratios()

## Export the two data sets

peptide_types <- "long"

write.csv2(induced_set_long, 
           file = file.path(dirname(directory), 
                            paste0("induced_set_", 
                                   peptide_types, 
                                   ".csv")),
           row.names = FALSE)

write.csv2(repressed_set_long, file = file.path(dirname(directory), paste("repressed_set", peptide_types, ".csv",
                                                                           sep = "")),row.names = FALSE)

## Make a merged long set and write it

count_set <- rbind(induced_set_long %>% 
                     mutate(induction = "induced"),
                   repressed_set_long %>%
                     mutate(induction = "repressed")) %>%
  mutate(induction = as.factor(induction))

## Randomize the order of the data frame so I can just grab a part of it.

count_set <- count_set[sample(nrow(count_set)), ]

write.csv2(count_set, 
           file = file.path(dirname(directory), paste("count_set_", peptide_types, ".csv",
                                                      sep = "")),
           row.names = FALSE)
