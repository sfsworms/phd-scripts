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

#Iterate the abose to create a DNAStringSet library


create.lib <- function(libSize = 10, pepLength = 4, alphabet = c("XXX","YYY")){
  lib <- c("1","2")
  for(i in 1:libSize){
    lib[i] <- generateSequence(pepLength, alphabet)
  }
  lib <- DNAStringSet(lib)
  return(lib)
}

#This function shoud take a ShortReadQ element, find a pattern in each of the reads, then return X nucleotides after the pattern

extract.peptide = function(dnaSeq, pattern = "TGGCTTCATTGCGAGCAAT", pepSize = 24) {
  targetAlignment = pairwiseAlignment(pattern = dnaSeq, subject = pattern,
    type = "local")  #Align the target

  pos = targetAlignment %>%
    pattern() %>%
    start() + nchar(pattern)  #This gives me the start of the peptide in each position
  sensiblePos = pos + pepSize < 148 #Added this for weird case where the pattern is present near the end
  
  targetlist <- subseq(x = dnaSeq[sensiblePos], start = pos[sensiblePos], width = pepSize)
  
  return(targetlist)
}
