#Functions for the peptide analysis project

#Load libraries needed
library(tidyverse)
library(Biostrings)
library(universalmotif) #Do generation of peptides

#Those are useful values
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

#Extract the peptide from a dnaSeq sequence. regexPattern should be the sequence right before the peptides. 

extract.peptide = function(dnaSeq, regexPattern = "TGGCTTCATTGCGAGCAAT", pepSize = 24){
  pepPosition <- str_locate(dnaSeq, regexPattern)[,2]
  sensiblePos <- pepPosition+pepSize < 145 #Added this for weird case where the pattern is present near the end
  peptideList <- subseq(x = dnaSeq[sensiblePos], start = pepPosition[sensiblePos]+1, end = pepPosition[sensiblePos]+pepSize)
  return(peptideList)
}

#Extract all peptides from the fastq of a run. Takes the name of the file, the destination where tos tore the peptide,
#the sequence just in front and behing the peptide and the size of the peptides.
extract.peptides.fastq = function(fileName, 
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
  stream <- FastqStreamer(fileName)
  on.exit(close(stream))
  
  #This is just to track the speed of processing of the file.
  i = 1
  ptm = proc.time()
  
  dnaSeq <- yield(stream)
  
  while(length(dnaSeq) > 0){
    print(c(i,fileName))
    print(proc.time() - ptm)
    i = i+1
    
    # Filter the one that are way too small (<145 bp out of 150, about 0.5% of seq for
    # NNK7)
    dnaSeq <- yield(stream)
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
    # intein sequence anyway#Take their reverse complement: problem, we don't get
    # barcode on those reads
    
    frontPos = str_locate(revDnaSeq, frontPattern)[,2]+1
    
    backPos = str_locate(revDnaSeq, backPattern)[,1]-1
    
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
    
    writeFasta(shortPep, 
               file.path(destination, 
                         gsub(basename(fileName), 
                              pattern = ".fastq.gz", 
                              replacement = paste("_peptide", 
                                                  shortPeptideSize, 
                                                  ".fa", sep=""))), 
               mode = "a")  # The mode append it to a file if existing
    
    writeFasta(longPep, 
               file.path(destination, 
                         gsub(basename(fileName), 
                              pattern = ".fastq.gz", 
                              replacement = paste("_peptide", 
                                                  largePeptideSize, 
                                                  ".fa", sep=""))), 
               mode = "a")
  }
}

## Compute sum of counts for each condition and then compute ratios of each read as a percentage of total
## compute an enrichment ratio based on the final ratio over the initial one, and removes NAs caused by divide by zero error

computeRatios <- function(df){
  
  df <- df %>%
    mutate(gen1_sum = gen1_r1 + gen1_r2) %>%
    mutate(gen5_sum = gen5_r1 + gen5_r2)
  
  df <- df %>%
    mutate(gen1_ratio = gen1_sum / sum(df$gen1_sum)) %>%
    mutate(gen5_ratio = gen5_sum / sum(df$gen5_sum)) %>%
    mutate(enrichment_ratio = gen5_ratio / gen1_ratio)
  
  # Drop the ratios that are infinite (where there is nothing originally) into NA
  df <- df %>%
    mutate(enrichment_ratio = ifelse(is.infinite(enrichment_ratio), NA, enrichment_ratio)) 
  # Add a log version
  df <- df %>%
    mutate(enrichment_ratio_log = log2(enrichment_ratio))
  
  return(df)
}

  
#Create a function that drop all the peptides with stop codons from a data frame
remove_stop_codons <- function(peptide_data_frame){
  output_df <- peptide_data_frame[!grepl("\\*", peptide_data_frame$peptide_seq),]
  return(output_df)
}
  
  

## Get all the sequence counts from within a directory into one set. Requires a directory, a list of files to load and the names of the variables in the set
createCountSet <- function(directory, file_list, run_names){
  for (i in seq_along(file_list)) {
    assign(run_names[i],
           read.csv(file.path(directory, file_list[i]), header = FALSE) %>%
             setNames(., c("seq", run_names[i]))
    )
  }
  ## Merge by sequence to form one big dataset.
  
  merged_set = get(run_names[1])
  
  for (i in 2:length(file_list)) {
    merged_set = full_join(merged_set, get(run_names[i]), by = "seq")
  }
  
  rm(run_names)
  #Replace NA by zeros
  merged_set <- merged_set %>%
    mutate_at(c(2:ncol(merged_set)),
              ~replace_na(.,0))
  
  # Add a peptide_seq that contains the translation
  merged_set <- merged_set %>%
    mutate(peptide_seq = seq %>%
             dna() %>%
             seq_translate() %>%
             as.character()) %>%
    relocate(peptide_seq, .after = seq)
  gc()
  return(merged_set)
}