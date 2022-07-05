#Functions for the peptide analysis project

#Load libraries
library(tidyverse)
library(Biostrings)
library(universalmotif)

codonAlphabet <- names(GENETIC_CODE)
codonNNK <- codonAlphabet[grep("(.)(.)[G,T]$",codonAlphabet)]
codonNNB <- codonAlphabet[grep("(.)(.)[C,G,T]$",codonAlphabet)]

#Generate a peptide from a table of codons. Might need to multithread to make it faster

generate.sequence <- function(seqLength = 4, alphabet = c("XXX","YYY")){
  x <- runif(seqLength, min = 1, max = length(alphabet)) %>% round()
  seq <- alphabet[x] %>% toString() %>% gsub(", ","",.)
  return(seq)
}

#Iterate the above to create a DNAStringSet library

create.lib <- function(libSize = 10, pepLength = 4, alphabet = c("XXX","YYY")){
  lib <- c("1","2")
  for(i in 1:libSize){
    lib[i] <- generateSequence(pepLength, alphabet)
  }
  lib <- DNAStringSet(lib)
  return(lib)
}

#An attempt to implement extract.peptide faster 

extract.peptide = function(dnaSeq, regexPattern = "TGGCTTCATTGCGAGCAAT", pepSize = 24){
  x <- str_locate(dnaSeq, regexPattern)[,2]
  sensiblePos <- x+pepSize < 150 #Added this for weird case where the pattern is present near the end
  peptideList <- subseq(x = dnaSeq[sensiblePos], start = x[sensiblePos]+1, end = x[sensiblePos]+pepSize)
  return(peptideList)
}


# extract.peptides.NNK.fastq <- function(fileName, destination = sprintf("%s_subset", fl)) {
#   
#   frontPattern <- "TGGCTTCATTGCGAGCAAT"
#   backPattern <- "TGTCTGTCTTACG"
#   revPattern <- "GTCGTAAGACAGACA"
#   # Load a fastQ file with the peptides. The files are too big to be used
#   # entirely, so I'll have to use FastqStremer
#   
#   # I first need to get a connection established.
#   # open the connection
#   stream <- FastqStreamer(fileName)
#   on.exit(close(stream))
#   
#   i = 1
#   ptm = proc.time()
#   
#   dnaSeq <- yield(stream)
#   
#   while(length(dnaSeq) > 0){
#     print(i)
#     print(proc.time() - ptm)
#     i = i+1
#     ### Fix that
#     
#     # Filter the one that are way too small (<145 bp out of 150, about 0.5% of seq for
#     # NNK7)
#     dnaSeq <- yield(stream)
#     dnaSeq = dnaSeq[width(dnaSeq) > 145]
#     
#     # Convert to a DNAStringSet
#     
#     dnaSeq = sread(dnaSeq)
#     
#     # Get all the sequence that match the fwd pattern for intein
#     
#     fwdDnaSeq <- dnaSeq[grepl(pattern = frontPattern, dnaSeq %>%
#                                 as.character())]
#     
#     # I will subset the NNK7 and NNK7
#     
#     fwdNNK3seq = fwdDnaSeq[grepl(pattern = "^TGCA(.)*", fwdDnaSeq %>%
#                                    as.character())]
#     
#     fwdNNK7seq = fwdDnaSeq[grepl(pattern = "^AAAA(.)*", fwdDnaSeq %>%
#                                    as.character())]
#     
#     # And get the sequence of the peptides
#     
#     fwdNNK3peptide = extract.peptide(dnaSeq = fwdNNK3seq, 
#                                      regexPattern = frontPattern,
#                                      pepSize = 12)
#     
#     fwdNNK7peptide = extract.peptide(dnaSeq = fwdNNK7seq, 
#                                      regexPattern = frontPattern,
#                                      pepSize = 24)
#     
#     # Get all the sequence that match the rev pattern for intein and take their reverse
#     # complement
#     
#     revDnaSeq = dnaSeq[grepl(pattern = revPattern, dnaSeq %>%
#                                as.character())] %>%
#       reverseComplement()
#     
#     #Remove those that do not contain fronPattern or backPattern
#     revDnaSeq = revDnaSeq[grepl(pattern = frontPattern, revDnaSeq %>%
#                                   as.character())]
#     revDnaSeq = revDnaSeq[grepl(pattern = backPattern, revDnaSeq %>%
#                                   as.character())]
#     
#     # I do not have a barcode but I can get the size of the peptide by looking at the
#     # intein sequence anyway#Take their reverse complement: problem, we don't get
#     # barcode on those reads
#     
#     frontPos = str_locate(revDnaSeq, frontPattern)[,2]+1
#     
#     backPos = str_locate(revDnaSeq, backPattern)[,1]-1
#     
#     makeSense <- frontPos < backPos #Check the back is after the front
#     
#     revPeptide = subseq(x = revDnaSeq[makeSense], 
#                         start =  frontPos[makeSense], 
#                         end = backPos[makeSense])
#     
#     revNNK3peptide = revPeptide[width(revPeptide) == 12]  
#     
#     revNNK7peptide = revPeptide[width(revPeptide) == 24]  
#     
#     # I'll store all the sequence in a ShortRead object
#     
#     fwdNNK7peptide = ShortRead(sread = fwdNNK7peptide, id = BStringSet(rep("NNK7 fw",
#                                                                            length(fwdNNK7peptide))))
#     
#     revNNK7peptide = ShortRead(sread = revNNK7peptide, id = BStringSet(rep("NNK7 rv",
#                                                                            length(revNNK7peptide))))
#     
#     NNK7peptide <- append(fwdNNK7peptide, revNNK7peptide)
#     
#     
#     fwdNNK3peptide = ShortRead(sread = fwdNNK3peptide, id = BStringSet(rep("NNK3 fw",
#                                                                            length(fwdNNK3peptide))))
#     
#     revNNK3peptide = ShortRead(sread = revNNK3peptide, id = BStringSet(rep("NNK3 rv",
#                                                                            length(revNNK3peptide))))
#     
#     NNK3peptide <- append(fwdNNK3peptide, revNNK3peptide)
#     
#     writeFasta(NNK3peptide, file.path(destination, gsub(basename(fileName), pattern = ".fastq", replacement = "_peptide3.fa")), mode = "a")  # The mode append it to a file if existing
#     writeFasta(NNK7peptide, file.path(destination, gsub(basename(fileName), pattern = ".fastq", replacement = "_peptide7.fa")), mode = "a")
#   }
# }


extract.peptides.fastq <- function(fileName, destination, libChar) {
  
  frontPattern = libChar$frontPattern
  backPattern = libChar$backPattern

  # Load a fastQ.gz file with the peptides. The files are too big to be used
  # entirely, so I'll have to use FastqStremer. .gz compressed files actually process faster than .fastq files.
  
  # I first need to get a connection established.
  # open the connection
  stream <- FastqStreamer(fileName)
  on.exit(close(stream))
  
  i = 1
  ptm = proc.time()
  
  dnaSeq <- yield(stream)
  
  while(length(dnaSeq) > 0){
    print(c(i,fileName))
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
    
    fwdShortSeq = fwdDnaSeq[grepl(pattern = "^TGCA(.)*", fwdDnaSeq %>%
                                   as.character())]
    
    fwdLongSeq = fwdDnaSeq[grepl(pattern = "^AAAA(.)*", fwdDnaSeq %>%
                                   as.character())]
    
    # And get the sequence of the peptides
    
    fwdShortPep = extract.peptide(dnaSeq = fwdShortSeq, 
                                     regexPattern = frontPattern,
                                     pepSize = 12)
    
    fwdLongPep = extract.peptide(dnaSeq = fwdLongSeq, 
                                     regexPattern = frontPattern,
                                     pepSize = 24)
    
    # Get all the sequence that match the rev pattern for intein and take their reverse
    # complement
    
    revDnaSeq = dnaSeq[grepl(pattern = libChar$backPattern %>% 
                               DNAString() %>%
                               reverseComplement(), 
                             dnaSeq %>%
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
    
    revShortPep = revPeptide[width(revPeptide) == libChar$shortPeptideSize]  
    
    revLongPep = revPeptide[width(revPeptide) == libChar$largePeptideSize]  
    
    # I'll store all the sequence in a ShortRead object
    
    fwdLongPep = ShortRead(sread = fwdLongPep, id = BStringSet(rep(paste(libChar$largePeptideSize,"AA peptide fw"),
                                                                           length(fwdLongPep))))
    
    revLongPep = ShortRead(sread = revLongPep, id = BStringSet(rep(paste(libChar$largePeptideSize,"AA peptide rv"),
                                                                           length(revLongPep))))
    
    longPep <- append(fwdLongPep, revLongPep)
    
    
    fwdShortPep = ShortRead(sread = fwdShortPep, id = BStringSet(rep(paste(libChar$shortPeptideSize,"AA peptide fw"),
                                                                           length(fwdShortPep))))
    
    revShortPep = ShortRead(sread = revShortPep, id = BStringSet(rep(paste(libChar$shortPeptideSize,"AA peptide rv"),
                                                                           length(revShortPep))))
    
    shortPep <- append(fwdShortPep, revShortPep)
    
    writeFasta(shortPep, 
               file.path(destination, 
                         gsub(basename(fileName), 
                              pattern = ".fastq.gz", 
                              replacement = paste("_peptide", 
                                                  libChar$shortPeptideSize, 
                                                  ".fa", sep=""))), 
               mode = "a")  # The mode append it to a file if existing
    
    writeFasta(longPep, 
               file.path(destination, 
                         gsub(basename(fileName), 
                              pattern = ".fastq.gz", 
                              replacement = paste("_peptide", 
                                                  libChar$longPeptideSize, 
                                                  ".fa", sep=""))), 
               mode = "a")
  }
}
