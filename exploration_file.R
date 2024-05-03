# This is an exploration file

## Library loading
if (!require("tidyverse")) {
  install.packages("tidyverse")
  library(tidyverse)
}

if (!require("ShortRead")) {
  install.packages("ShortRead")
  library(ShortRead)
}

## Setting the folder to explore
input_directory <- file.path("C:/Users/Sebastian Worms/ngs_data/2023 NGS Ale/90-933598625/sample_100k/filtered")

file_list <- list.files(input_directory) %>%  #Get files from directory
  grep(pattern = ".fastq.gz$", ., value = TRUE) #only the .Fastq

## Grabs all the reads 
# TODO: Make it to take in chunks
fastq_fw <- readFastq(file.path(input_directory,file_list[1]))
fastq_rv <- readFastq(file.path(input_directory,file_list[2]))

## Get the reads from the fastQ
reads_fw <- sread(fastq_fw)
reads_rv <- sread(fastq_rv)

## Exploring the merged data

merged_directory <- file.path("C:/Users/Sebastian Worms/ngs_data/2023 NGS Ale/90-933598625/sample_100k/filtered/merged_reads")
merged_file_list <- list.files(merged_directory) %>%  #Get files from directory
  grep(pattern = ".fastq.gz$", ., value = TRUE) #only the .Fastq

fastq_merged <- readFastq(file.path(merged_directory,merged_file_list[1]))

reads_merged <- sread(fastq_merged)

width(reads_merged) %>% hist()
