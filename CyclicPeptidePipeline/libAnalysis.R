## Analysis project

setwd(dir = "C:/Users/worms/Dropbox/PhD/PhD-Scripts/CyclicPeptidePipeline")
source("library.R")  #Install and/or load needed libraries
source("function.R")  #Functions used in the script

# Load a fastQ file with the peptides. The files are too big to be used
# entirely, so I'll have to use FastqStremer

# I first need to get a connection established.

fileName <- file.path("C:/Users/worms/Documents/2022.06.07 Drift Seq/90-666155004b/Cytoplasmic-NNK-Gen-5-Ara_R1_001.fastq")

destination <- file.path("C:/Users/worms/Desktop/Test peptide2")


extract.peptides.fastq <- function(fileName, destination = sprintf("%s_subset", fl)) {

  frontPattern <- "TGGCTTCATTGCGAGCAAT"
  backPattern <- "TGTCTGTCTTACG"
  
  # open the connection
  stream <- open(FastqStreamer(fileName, verbose = TRUE))
  on.exit(close(stream))
  
  i = 1
  ptm = proc.time()
  repeat{
    
    print(i)
    print(proc.time() - ptm)
    i = i+1
    
    dnaSeq = yield(stream) 

    # Filter the one that are way too small (<145 bp out of 150, about 0.5% of seq for
    # NNK7)

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

    fwdNNK3peptide = extract.peptide(dnaSeq = fwdNNK3seq , pattern = frontPattern,
      pepSize = 12)
    fwdNNK7peptide = extract.peptide(dnaSeq = fwdNNK7seq , pattern = frontPattern,
      pepSize = 24)

    # Get all the sequence that match the rev pattern for intein and take their reverse
    # complement

    revDnaSeq = dnaSeq[grepl(pattern = "GTCGTAAGACAGACA", dnaSeq %>%
      as.character())] %>%
      reverseComplement()

    # I do not have a barcode but I can get the size of the peptide by looking at the
    # intein sequence anyway#Take their reverse complement: problem, we don't get
    # barcode on those reads

    #I'm running alignment multiple times here, this ain't effficienc

    frontPos = pairwiseAlignment(pattern = revDnaSeq, subject = frontPattern, type = "local") %>%
      pattern() %>%
      start() + nchar(frontPattern)

    backPos = pairwiseAlignment(pattern = revDnaSeq, subject = backPattern, type = "local") %>%
      pattern() %>%
      start()

    pepLength = backPos - frontPos

    revNNK3peptide = extract.peptide(dnaSeq = revDnaSeq[pepLength == 12], pattern = "TGGCTTCATTGCGAGCAAT",
      pepSize = 12)

    revNNK7peptide = extract.peptide(dnaSeq = revDnaSeq[pepLength == 24], pattern = "TGGCTTCATTGCGAGCAAT",
      pepSize = 24)

    rm(backPos, frontPos, pepLength)

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

    writeFasta(NNK7peptide, paste(destination, "NNK7peptide.fa", sep = "/"), mode = "a")  # The mode append it to a file if existing
    writeFasta(NNK3peptide, paste(destination, "NNK3peptide.fa", sep = "/"), mode = "a")
  }
}

extract.peptides.fastq(fileName = fileName, destination = destination)
