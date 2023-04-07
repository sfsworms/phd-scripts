## This file serves to extract the peptide sequences from filtered and merged files.

setwd(dir = "C:/Users/worms/Dropbox/phd_script/cyclic_peptide_pipeline")
source("functions/peptide_extraction_NNK_functions.R")  #Functions used in the script

# Libraries need
library(ShortRead)
library(tidyverse)

# Directory store the directory with all the files.
directory <- "C:/Users/worms/NGS Data/2022.06.07_drift_seq/90-666155004b/00_fastq/NNK/filtered_merged"

#Get a list of the merged .fastq files
file_list <- list.files(directory) %>%  #Get files from directory
  grep(pattern = ".fastq.gz$", ., value = TRUE) #only the .Fastq
 

destination <- file.path("C:/Users/worms/NGS Data/2022.06.07_drift_seq/90-666155004b/00_fastq/NNK") #Set where to put the fasta with peptides files.

tracking_df <- data.frame(n_reads = integer(),
                          n_short_peptide = integer(),
                          n_long_peptide = integer(),
                          file = character())

for(i in seq_along(file_list)){ #Run the extraction script on all the fastq NNK files
progress_df <-  extract.peptides.fastq(file_name = file.path(directory,file_list[i]),
                         destination = file.path(destination,"peptide_sequence"))

tracking_df <- rbind(tracking_df, progress_df)
}

write.csv2(tracking_df, file.path(destination,"tracking_NNB.csv"))
