## Analysis project

#Note to self: investigate if garbage collection 'gc()' can speed up the loop. Probably a good practice to insert it after function taking > 30s.

setwd(dir = "C:/Users/worms/Dropbox/PhD/PhD-Scripts/CyclicPeptidePipeline")
source("function.R")  #Functions used in the script

# Directory store the directory with all the files.

directory <- "C:/Users/worms/NGS Data/2022.06.07_drift_seq/90-666155004b/00_fastq/NNK"
#Get a list of the NNK .fastq files
file_list <- list.files(directory) %>%  #Get files from directory
  grep(pattern = ".fastq.gz$", ., value = TRUE) #only the .Fastq
 

destination <- file.path("C:/Users/worms/Dropbox/PhD/PhD-Scripts/CyclicPeptidePipeline/test") #Set where to put the fasta with peptides files.

for(i in seq_along(file_list)){ #Run the extraction script on all the fastq NNK files
  extract.peptides.fastq(file_name = file.path(directory,file_list[i]),
                         destination = file.path(destination,"peptide_sequence"))
}
