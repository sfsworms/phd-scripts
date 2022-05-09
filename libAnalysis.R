#So, let me play with fake libs.

#Load libraries
library(tidyverse)
library(Biostrings)

#Create functions I want to create a random peptide
createPeptide <- function(size = 8){
  aaList <- AA_ALPHABET[c(1:20,27)]
  peptide <- aaList[sample(1:length(aaList), size)] %>% 
    unlist %>%
    paste(.,collapse="")
  return(peptide)
}

createLibrary <- function(n = 1, size = 8){
  library <- vector()
  for(i in 1:n){
  library <- c(library,createPeptide(size))
  } 
  library <- AAStringSet(x = library)
  return(library)
}


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
  
BiocManager::install("universalmotif")


