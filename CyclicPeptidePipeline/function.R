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
