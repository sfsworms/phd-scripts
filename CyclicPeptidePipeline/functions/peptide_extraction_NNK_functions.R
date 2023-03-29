#Extract the peptide from a dnaSeq sequence. regexPattern should be the sequence right before the peptides. 

extract.peptide <- function(dnaSeq, regexPattern = "TGGCTTCATTGCGAGCAAT", pepSize = 24){
  pepPosition <- str_locate(dnaSeq %>% as.character(), regexPattern)[,2]
  sensiblePos <- pepPosition+pepSize < 145 #Added this for weird case where the pattern is present near the end
  peptideList <- subseq(x = dnaSeq[sensiblePos], start = pepPosition[sensiblePos]+1, end = pepPosition[sensiblePos]+pepSize)
  rm(pepPosition, sensiblePos)
  return(peptideList)
}

#Extract all peptides from the fastq of a run. Takes the name of the file, the destination where tos tore the peptide,
#the sequence just in front and behing the peptide and the size of the peptides.

extract.peptides.fastq = function(file_name, 
                                  destination, 
                                  frontPattern = "TTCATTGCGAGCAAT",
                                  backPattern = "TGTCTGTCTTACGACA",
                                  shortPeptideSize = 12,
                                  largePeptideSize = 24
) {
  # Load a fastQ.gz file with the peptides. The files are too big to be used
  # entirely, so I'll have to use FastqStremer. .gz compressed files actually process faster than .fastq files.
  
  # I first need to get a connection established.
  # open the connection
  reads_per_batch <- 10^6
  stream <- FastqStreamer(file_name, n = reads_per_batch)
  on.exit(close(stream))
  
  #This is just to track the speed of processing of the file.
  i = 1
  ptm = proc.time()
  
  dnaSeq <- yield(stream)
  
  # This DF will be used to write out the total number of peptides produced at the end.
  tracking_df <- data.frame(n_reads = 0,
                            n_short_peptide = 0,
                            n_long_peptide = 0)
  
  # This bit below checks before the while loop that output doesn't exist to avoid appending to an existing file
  
  if (file.exists(file.path(destination, 
                        gsub(basename(file_name), 
                             pattern = ".fastq.gz", 
                             replacement = paste("_peptide", 
                                                 largePeptideSize, 
                                                 ".fa", sep=""))))){
    stop("The output file for large peptides already exists.")
  }
  
  if (file.exists(file.path(destination, 
                            gsub(basename(file_name), 
                                 pattern = ".fastq.gz", 
                                 replacement = paste("_peptide", 
                                                     shortPeptideSize, 
                                                     ".fa", sep=""))))){
    stop("The output file for short peptides already exists.")
  }
  
  while(length(dnaSeq) > 0){
    # Filter the one that are way too small (<145 bp out of 150, about 0.5% of seq for
    # NNK7)
    length_dnaSeq <- length(dnaSeq)
    tracking_df <- tracking_df %>% mutate(n_reads = n_reads+length_dnaSeq) # Tracking purpose
    
    dnaSeq = dnaSeq[width(dnaSeq) > 145]
    
    # Convert to a DNAStringSet
    
    dnaSeq = sread(dnaSeq)
    
    # Get all the sequence that match the fwd pattern for intein
    
    fwdDnaSeq <- dnaSeq[grepl(pattern = frontPattern, dnaSeq %>%
                                as.character())]
    
    # I will subset the NNK7 and NNK3
    
    fwdShortSeq = fwdDnaSeq[grepl(pattern = "^TGCA(.)*", fwdDnaSeq %>%
                                    as.character())]
    
    fwdLongSeq = fwdDnaSeq[grepl(pattern = "^AAAA(.)*", fwdDnaSeq %>%
                                   as.character())]
    
    # And get the sequence of the peptides
    
    fwdShortPep = extract.peptide(dnaSeq = fwdShortSeq, 
                                  regexPattern = frontPattern,
                                  pepSize = shortPeptideSize)
    
    fwdLongPep = extract.peptide(dnaSeq = fwdLongSeq, 
                                 regexPattern = frontPattern,
                                 pepSize = largePeptideSize)
    
    # Get all the sequence that match the rev pattern for intein and take their reverse
    # complement
    
    revDnaSeq = dnaSeq[grepl(pattern = backPattern %>% 
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
    # intein sequence anyway
    
    frontPos = str_locate(revDnaSeq %>% as.character(), frontPattern)[,2]+1
    
    backPos = str_locate(revDnaSeq %>% as.character(), backPattern)[,1]-1
    
    makeSense <- frontPos < backPos #Check the back is after the front
    
    revPeptide = subseq(x = revDnaSeq[makeSense], 
                        start =  frontPos[makeSense], 
                        end = backPos[makeSense])
    
    revShortPep = revPeptide[width(revPeptide) == shortPeptideSize]  
    
    revLongPep = revPeptide[width(revPeptide) == largePeptideSize]  
    
    # I'll store all the sequence in a ShortRead object
    
    fwdLongPep = ShortRead(sread = fwdLongPep, id = BStringSet(rep(paste(largePeptideSize,"AA peptide fw"),
                                                                   length(fwdLongPep))))
    
    revLongPep = ShortRead(sread = revLongPep, id = BStringSet(rep(paste(largePeptideSize,"AA peptide rv"),
                                                                   length(revLongPep))))
    
    longPep <- append(fwdLongPep, revLongPep)
    
    
    fwdShortPep = ShortRead(sread = fwdShortPep, id = BStringSet(rep(paste(shortPeptideSize,"AA peptide fw"),
                                                                     length(fwdShortPep))))
    
    revShortPep = ShortRead(sread = revShortPep, id = BStringSet(rep(paste(shortPeptideSize,"AA peptide rv"),
                                                                     length(revShortPep))))
    
    shortPep <- append(fwdShortPep, revShortPep)
    
    tracking_df <- tracking_df %>% mutate(n_short_peptide = n_short_peptide + length(shortPep),
                                          n_long_peptide = n_long_peptide + length(longPep),)
    
    
    writeFasta(shortPep, 
               file.path(destination, 
                         gsub(basename(file_name), 
                              pattern = ".fastq.gz", 
                              replacement = paste("_peptide", 
                                                  shortPeptideSize, 
                                                  ".fa", sep=""))), 
               mode = "a")  # The mode append it to a file if existing
    
    writeFasta(longPep, 
               file.path(destination, 
                         gsub(basename(file_name), 
                              pattern = ".fastq.gz", 
                              replacement = paste("_peptide", 
                                                  largePeptideSize, 
                                                  ".fa", sep=""))), 
               mode = "a")
    dnaSeq <- yield(stream)
    
    # Give progress info
    
    paste("Out of ", length_dnaSeq," reads in batch ", i, " of file ",basename(file_name)," ",length(shortPep)+length(longPep)," encoded for peptides.") %>% print()
    print(proc.time() - ptm)
    i = i+1
    gc()
  }
  #Load the next batch of reads
  
  tracking_df %>% mutate(file = basename(file_name))
  
  return(tracking_df)
  
}