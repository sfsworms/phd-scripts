#TODO: Turn the sampling into a function

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

if (!require("fs")) {
  install.packages("fs")
  library(fs)
}

## Inputs
### Folder to the data
input_directory <- file.path("C:/Users/Sebastian Worms/ngs_data/NGS Ale/90-933598625/00_fastq") #This should link to the folder with the .fastq.gz
output_directory <- file.path("C:/Users/Sebastian Worms/ngs_data/NGS Ale/90-933598625/00_fastq","sampled_fastq") #Directory where the samples will be stored

### Number of sequences to sample
n = 10e5

# Check if the subfolder exists
if (!dir_exists(output_directory )) {
  # Create the subfolder if it does not exist
  dir_create(output_directory )
}

## Script

### Get the files

#Get a list of the .fastq files.
file_list <- list.files(input_directory) %>%  #Get files from directory
  grep(pattern = ".fastq.gz$", ., value = TRUE) #only the .Fastq

for(i in file_list){
  stream <- FastqStreamer(file.path(input_directory,i), n = n)
  on.exit(close(stream))
  fq <- yield(stream)
  writeFastq(fq, file.path(output_directory,
                           paste0(gsub(".fastq.gz",
                                       "_sample.fastq.gz",
                                       i)
                                  )
                           )
             )
}

