#Functions for the peptide analysis project

#Load libraries
library(tidyverse)
library(Biostrings)
library(universalmotif)

codonAlphabet <- names(GENETIC_CODE)
codonNNK <- codonAlphabet[grep("(.)(.)[G,T]$",codonAlphabet)]
codonNNB <- codonAlphabet[grep("(.)(.)[C,G,T]$",codonAlphabet)]

#Generate a peptide from a table of codons. Might need to multithread to make it faster

generateSequence <- function(seqLength = 4, alphabet = c("XXX","YYY")){
  x <- runif(seqLength, min = 1, max = length(alphabet)) %>% round()
  seq <- alphabet[x] %>% toString() %>% gsub(", ","",.)
  return(seq)
}

#Iterate the abose to create a DNAStringSet library


createLib <- function(libSize = 10, pepLength = 4, alphabet = c("XXX","YYY")){
  lib <- c("1","2")
  for(i in 1:libSize){
    lib[i] <- generateSequence(pepLength, alphabet)
  }
  lib <- DNAStringSet(lib)
  return(lib)
}
