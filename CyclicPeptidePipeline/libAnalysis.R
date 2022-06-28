## Analysis project

setwd(dir = "C:/Users/worms/Dropbox/PhD/PhD-Scripts/CyclicPeptidePipeline")
source("library.R")  #Install and/or load needed libraries
source("function.R")  #Functions used in the script

# Directory store the directory with all the files.




extract.peptides.fastq <- function(fileName, destination = sprintf("%s_subset", fl)) {

  frontPattern <- "TGGCTTCATTGCGAGCAAT"
  backPattern <- "TGTCTGTCTTACG"
  revPattern <- "GTCGTAAGACAGACA"
  # Load a fastQ file with the peptides. The files are too big to be used
  # entirely, so I'll have to use FastqStremer
  
  # I first need to get a connection established.
  # open the connection
  stream <- FastqStreamer(fileName)
  on.exit(close(stream))
  
  i = 1
  ptm = proc.time()
  
  dnaSeq <- yield(stream)
  
  while(length(dnaSeq) > 0){
    print(i)
    print(proc.time() - ptm)
    i = i+1
    ### Fix that
    
    # Filter the one that are way too small (<145 bp out of 150, about 0.5% of seq for
    # NNK7)
    dnaSeq <- yield(stream)
    dnaSeq = dnaSeq[width(dnaSeq) > 145]

    # Convert to a DNAStringSet
    
    dnaSeq = sread(dnaSeq)
    
    # Get all the sequence that match the fwd pattern for intein

    fwdDnaSeq <- dnaSeq[grepl(pattern = frontPattern, dnaSeq %>%
      as.character())]

    # I will subset the NNK7 and NNK7

    fwdNNK3seq = fwdDnaSeq[grepl(pattern = "^TGCA(.)*", fwdDnaSeq %>%
      as.character())]

    fwdNNK7seq = fwdDnaSeq[grepl(pattern = "^AAAA(.)*", fwdDnaSeq %>%
      as.character())]

    # And get the sequence of the peptides

    fwdNNK3peptide = extract.peptide(dnaSeq = fwdNNK3seq, 
                                     regexPattern = frontPattern,
                                     pepSize = 12)
    
    fwdNNK7peptide = extract.peptide(dnaSeq = fwdNNK7seq, 
                                     regexPattern = frontPattern,
                                     pepSize = 24)

    # Get all the sequence that match the rev pattern for intein and take their reverse
    # complement

    revDnaSeq = dnaSeq[grepl(pattern = revPattern, dnaSeq %>%
      as.character())] %>%
      reverseComplement()
    
    #Remove those that do not contain fronPattern or backPattern
    revDnaSeq = revDnaSeq[grepl(pattern = frontPattern, revDnaSeq %>%
                                  as.character())]
    revDnaSeq = revDnaSeq[grepl(pattern = backPattern, revDnaSeq %>%
                                  as.character())]
    
    # I do not have a barcode but I can get the size of the peptide by looking at the
    # intein sequence anyway#Take their reverse complement: problem, we don't get
    # barcode on those reads
    
    frontPos = str_locate(revDnaSeq, frontPattern)[,2]+1

    backPos = str_locate(revDnaSeq, backPattern)[,1]-1

    makeSense <- frontPos < backPos #Check the back is after the front
    
    revPeptide = subseq(x = revDnaSeq[makeSense], 
                        start =  frontPos[makeSense], 
                        end = backPos[makeSense])
    
    revNNK3peptide = revPeptide[width(revPeptide) == 12]  
      
    revNNK7peptide = revPeptide[width(revPeptide) == 24]  

    # I'll store all the sequence in a ShortRead object

    fwdNNK7peptide = ShortRead(sread = fwdNNK7peptide, id = BStringSet(rep("NNK7 fw",
      length(fwdNNK7peptide))))

    revNNK7peptide = ShortRead(sread = revNNK7peptide, id = BStringSet(rep("NNK7 rv",
      length(revNNK7peptide))))

    NNK7peptide <- append(fwdNNK7peptide, revNNK7peptide)


    fwdNNK3peptide = ShortRead(sread = fwdNNK3peptide, id = BStringSet(rep("NNK3 fw",
      length(fwdNNK3peptide))))

    revNNK3peptide = ShortRead(sread = revNNK3peptide, id = BStringSet(rep("NNK3 rv",
      length(revNNK3peptide))))

    NNK3peptide <- append(fwdNNK3peptide, revNNK3peptide)

    writeFasta(NNK7peptide, file.path(destination, gsub(basename(fileName), pattern = ".fastq", replacement = "_peptide3.fa")), mode = "a")  # The mode append it to a file if existing
    writeFasta(NNK3peptide, file.path(destination, gsub(basename(fileName), pattern = ".fastq", replacement = "_peptide7.fa")), mode = "a")
  }
}

directory <- "D:/2022.06.07 Drift Seq/90-666155004b/00_fastq/NNK"

#Get a list of the .fastq files
fileList <- list.files(directory) %>% 
            grep(pattern = ".fastq$", ., value = TRUE)

for(i in seq_along(fileList)){
  extract.peptides.fastq(fileName = file.path(directory,fileList[i]), destination = file.path(directory,"Destination"))
}

fileName <- file.path(directory,fileList[1])

destination <- file.path("C:/Users/worms/Desktop/Test peptide2")

list.files(directory)

grep(pattern = ".fastq$", fileList, value = TRUE)


