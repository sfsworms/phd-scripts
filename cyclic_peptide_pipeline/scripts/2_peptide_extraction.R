## This file serves to extract the peptide sequences from filtered and merged files.

# Set working directory to the script's location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Source external functions
source(file.path(
   dirname(getwd()),
  "functions",
  "peptide_extraction_functions.R"
    )
  ) 

# List of required packages
required_packages <- c("ShortRead", "tidyverse", "rstudioapi")

# Install and load required packages
lapply(required_packages, install_if_missing)

#Source directory
directory <- "C:/Users/sebas/2023_OXA_NGS/90-957741147/demultiplex/Npu-final"

#Get a list of the merged .fastq files
file_list <- list.files(directory) %>%  #Get files from directory
  grep(pattern = ".fastq.gz$", ., value = TRUE) #only the .Fastq
 

destination <- file.path("C:/Users/sebas/2023_OXA_NGS/90-957741147/npu_pep") #Set where to put the fasta with peptides files.

tracking_df <- data.frame(n_reads = integer(),
                          n_peptide = integer(),
                          file = character())

for(i in seq_along(file_list)){ #Run the extraction script on all the fastq NNK files
      progress_df <-  extract.peptides.fastq(file_name = file.path(directory,file_list[i]),
                             destination = file.path(destination,"peptide_sequence"),
                             front_pattern = "CTTCATTGCGAGCAAT",
                             back_pattern = "TGTCTGTCTTACGACACCG",
                             short_peptide_size = 2,
                             long_peptide_size = 18)

tracking_df <- rbind(tracking_df, progress_df)
}

write.csv2(tracking_df, file.path(destination,"tracking_NNB.csv"))
