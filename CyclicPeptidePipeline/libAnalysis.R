## Analysis project

setwd(dir = "C:/Users/worms/Dropbox/PhD/PhD-Scripts/CyclicPeptidePipeline")
source("library.R")  #Install and/or load needed libraries
source("function.R")  #Functions used in the script

# Directory store the directory with all the files.

directory <- "D:/2022.06.07 Drift Seq/90-666155004b/00_fastq/NNK"
#Get a list of the NNK .fastq files
fileList <- list.files(directory) %>%  #Get files from directory
  grep(pattern = ".fastq$", ., value = TRUE) %>% #only the .Fastq
  grep(pattern = "NNK", ., value = TRUE)  #that have NNK in the name

destination <- file.path("C:/Users/worms/Desktop/Test peptide2") #Set where to put the fasta with peptides files.


for(i in seq_along(fileList)){ #Run the extraction script on all the fastq NNK files
  extract.peptides.fastq(fileName = file.path(directory,fileList[i]), destination = file.path(directory,"Destination"))
}

### Get sequence counts from fasta file of peptides
peptideFileList <- list.files(file.path(directory,"Destination"))


#At this point I think I'll use biopython to edit my script
stream <- FastqStreamer(file.path(directory,"Destination",peptideFileList[1]))


x <- readAAStringSet(file.path(directory,
                               "Destination",
                               peptideFileList[1]))

