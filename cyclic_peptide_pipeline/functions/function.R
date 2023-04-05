#Functions for the peptide analysis project

#Load libraries needed
library(tidyverse)
library(Biostrings)
library(universalmotif) #Do generation of peptides
library(ShortRead) #Needed to read the large fastq files

#Those are useful values
# codonAlphabet <- names(GENETIC_CODE)
# codonNNK <- codonAlphabet[grep("(.)(.)[G,T]$",codonAlphabet)]
# codonNNB <- codonAlphabet[grep("(.)(.)[C,G,T]$",codonAlphabet)]

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
    gc()
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
