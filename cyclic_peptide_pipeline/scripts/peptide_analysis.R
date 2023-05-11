# This script is an ETL aiming to take counts of peptides from the csv of counts produced by the
# python scripts and make them into one .csv file containing a dataframe ready for exploratory
# analysis. That data frame has counts for all conditions as well as ratio of total read per
# conditions and enrichment ratios.

# Setup

setwd(dir = "C:/Users/worms/phd-scripts/cyclic_peptide_pipeline")  # Directory containing the scripts
source("functions/function_analysis.R")  #Functions used in the script

## Load necessary libraries
library(tidyverse)  #Needed for data wrangling
library(bioseq)  # Used to translate sequence

## This folder should contain the CSV files with a column for sequence and a column for counts
directory <- "C:/Users/worms/ngs_data/2022_06_07_drift_seq/90-666155004b/00_fastq/all_files/peptide_count_csv"

## Get the names of the counts .csv files
file_list <- list.files(directory)
file_list <- file_list[grepl(".csv", file_list)]

## Trim the file names to use as variable names.
run_names <- file_list %>%
  gsub("peptide", "", .) %>%
  gsub("_count.csv", "", .) %>%
  tolower() %>%
  gsub("gen_", "gen", .)

## Split the file names and run names based on the peptide size

run_names_short <- run_names[grepl("12", run_names)]
run_names_long <- run_names[grepl("24", run_names)]

file_list_short <- file_list[grepl("12", file_list)]
file_list_long <- file_list[grepl("24", file_list)]

# Short set

## Import the short set
merged_set_short <- create_count_set(directory, file_list_short, run_names_short)

# Pivot it to a longer format with a library and a condition columns
merged_set_short <- merged_set_short %>% 
  pivot_longer(cols = c(3:8),names_to = "run_name", values_to = "count") %>%
  separate_wider_regex(cols = "run_name", c(library = ".*?", "_", experiment = ".*")) %>%
  pivot_wider(id_cols = c("seq", "peptide_seq", "library"), names_from = "experiment", values_from = "count") %>%
  pivot_longer( cols = contains("gen5"), names_to = "condition", values_to = "gen5") %>% 
  relocate("condition", .after = 3) %>%
  mutate(condition = str_extract(condition, "(?<=_)[^_]{3}(?=_)")) %>%
  rename("gen1_lb_12" = "gen1")

## Compute ratios for each conditions 

merged_set_short <- merged_set_short %>%
  compute_ratios(df = .,
                 count.cols = c(5,6),
                 info.cols = c(3,4))


# Rename the induction

merged_set_short <- merged_set_short %>%
  mutate(condition = ifelse(condition == "ara", "induced", "repressed"))


## Add a 'standard seq' column to take into account the cyclisation

merged_set_short <- merged_set_short  %>%
  mutate(standard_seq = standard_sequence(peptide_seq), .after = peptide_seq)

# Randomize the order 

merged_set_short <- merged_set_short[sample(nrow(merged_set_short)), ]

# Write the set

write.csv2(merged_set_short, 
           file = file.path(dirname(directory), "count_set_short.csv"),
           row.names = FALSE)


# Long set

## Import the long set
merged_set_long <- create_count_set(directory, file_list_long, run_names_long)

## Pivot it to a longer format with a library and a condition columns
merged_set_long <- merged_set_long %>% 
  pivot_longer(cols = c(3:8),names_to = "run_name", values_to = "count") %>%
  separate_wider_regex(cols = "run_name", c(library = ".*?", "_", experiment = ".*")) %>%
  pivot_wider(id_cols = c("seq", "peptide_seq", "library"), names_from = "experiment", values_from = "count") %>%
  pivot_longer( cols = contains("gen5"), names_to = "condition", values_to = "gen5") %>% 
  relocate("condition", .after = 3) %>%
  mutate(condition = str_extract(condition, "(?<=_)[^_]{3}(?=_)")) %>% # Extrace the glu or ara in the condition name
  rename("gen1_lb_24" = "gen1") %>% # shorten that column
  filter(gen1 > 0 | gen5 > 0) # Remove lines that have no reads

## Compute ratios for each conditions 
merged_set_long <- merged_set_long %>%
  compute_ratios(df = .,
                 count.cols = c(5,6),
                 info.cols = c(3,4))


## Rename the induction
merged_set_long <- merged_set_long %>%
  mutate(condition = ifelse(condition == "ara", "induced", "repressed"))


## Add a 'standard seq' column to take into account the cyclisation
merged_set_long <- merged_set_long  %>%
  mutate(standard_seq = standard_sequence(peptide_seq), .after = peptide_seq)

## Randomize the order 
merged_set_long <- merged_set_long[sample(nrow(merged_set_long)), ]

## Write the set
write.csv2(merged_set_long, 
           file = file.path(dirname(directory), "count_set_long.csv"),
           row.names = FALSE)
