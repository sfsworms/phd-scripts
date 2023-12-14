# The goal of this script is to sample bits of a very large fastq.gz file and save it for use with other script when running optimisation.

## Library loading
if (!require("tidyverse")) {
  install.packages("tidyverse")
  library(tidyverse)
}

if (!require("ShortRead")) {
  install.packages("ShortRead")
  library(ShortRead)
}

## Link to the data
input_directory <- file.path("") #This should link to the folder with the .fastq.gz
output_directory <- file.path("C:/Users/worms/ngs_data/2022_06_07_drift_seq/90-666155004b/00_fastq/NNB")

## Script

### Get the files

#Get a list of the merged .fastq files
file_list <- list.files(directory) %>%  #Get files from directory
  grep(pattern = ".fastq.gz$", ., value = TRUE) #only the .Fastq