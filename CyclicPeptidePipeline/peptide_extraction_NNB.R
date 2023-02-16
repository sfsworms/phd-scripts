## Analysis project

#Note to self: investigate if garbage collection 'gc()' can speed up the loop. Probably a good practice to insert it after function taking > 30s.

setwd(dir = "C:/Users/worms/Dropbox/PhD/PhD-Scripts/CyclicPeptidePipeline")
source("library.R")  #Install and/or load needed libraries
source("function.R")  #Functions used in the script

# Directory store the directory with all the files.

directory <- "C:/Users/worms/NGS Data/2022.06.07_drift_seq/90-666155004b/00_fastq/NNB"
#Get a list of the NNK .fastq files
file_list <- list.files(directory) %>%  #Get files from directory
  grep(pattern = ".fastq.gz$", ., value = TRUE) #only the .Fastq
 

destination <- file.path(directory,"NNB2") #Set where to put the fasta with peptides files.

for(i in seq_along(file_list)){ #Run the extraction script on all the fastq NNB files
  extract.peptides.fastq(fileName = file.path(directory,file_list[i]), 
                         destination = file.path(destination,"peptide_sequence"),
                         shortPeptideSize = 12,
                         largePeptideSize = 24)
}

### Get sequence counts from fasta file of peptides
peptidefile_list <- list.files(file.path(directory,"Destination"))


##Investigate the quality of NNB files

test <- FastqSampler(file.path(directory,file_list[1]))
seq <- yield(test)


x <- readAAStringSet(file.path(directory,
                               "Destination",
                               peptidefile_list[1]))

fileName = "D:/2022.06.07 Drift Seq/90-666155004b/00_fastq/Periplasmic-NNB-Gen-1-LB_R1_001.fastq.gz"

libCharNNK = list(frontPattern = "TGGCTTCATTGCGAGCAAT",
               backPattern = "TGTCTGTCTTACG",
               shortPeptideSize = 12,
               largePeptideSize = 24)

libCharNNB = list(frontPattern = "ATCATTGTCCATAAC",
                  backPattern = "TGCATCAGTGGAGAT",
                  shortPeptideSize = 18,
                  largePeptideSize = 30)

directory ="D:/2022.06.07 Drift Seq/90-666155004b/00_fastq"

file_list <- list.files("D:/2022.06.07 Drift Seq/90-666155004b/00_fastq/") %>%  #Get files from directory
  grep(pattern = ".fastq.gz$", ., value = TRUE) %>% #only the .Fastq
  grep(pattern = "NNB", ., value = TRUE)  #that have NNK in the name

for(i in seq_along(file_list)){ #Run the extraction script on all the fastq NNK files
  extract.peptides.fastq(fileName = file.path(directory,file_list[i]), destination = file.path(directory,"Destination"), libChar =  libCharNNB)
}


extract.peptides.fastq(fileName, 
                       destination = "D:/2022.06.07 Drift Seq/90-666155004b/00_fastq/test",
                       libChar = libCharNNB)


#Weird.

test <- readFasta(file.path(directory,"Destination","Periplasmic-NNB-Gen-5-Glu_R1_001_peptide.fa")) %>% sread()

alphabetFrequency(test)
alphabetByCycle(test)

#So the peptides look legit. There just aren't many of them'
#Let's look at the quality ranking.


qaSummary <- qa(file.path(directory,file_list[6]))



?freqSequence
?.freqSequences(qa, "read")

fls <- dir(directory, "*fastq.gz$", full=TRUE)
qaSummary <- qa(fls, type="fastq")
directory
