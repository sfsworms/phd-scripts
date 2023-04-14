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



# compute_ratios
## compute an enrichment ratio based on the final ratio over the initial one, and removes NAs caused by divide by zero error
## Assumed a merged count set where each row has count and is named library_generation_other

compute_ratios <- function(df){
  # Get all the colnames and decide new colnames
  col_names <- colnames(df)
  new_col_names <- paste0(col_names[-c(1:3)],"_ratio")
  
  # For all columns apart from the seq and peptide_seq, compute a ratio and add it to the df
  
  for(name in col_names[-c(1:3)]){
    column <- df %>% select(all_of(name))
    new_column <- column/sum(column)
    df <- cbind(df, new_column)
  }
  
  colnames(df) <- c(col_names, new_col_names)
  
  # Compute a gen5 over gen1 ratio
  
  gen1 <- df %>% select(contains("gen1") & contains("ratio"))
  gen5 <- df %>% select(contains("gen5") & contains("ratio"))
  enrichment_ratio <- gen5/gen1
  colnames(enrichment_ratio) <- "enrichment_ratio"
  
  df <- df %>%
    cbind(enrichment_ratio) %>% 
    mutate(enrichment_ratio_log = log2(enrichment_ratio))
    
  return(df)
}

# create_count_set()
## Get all the sequence counts from within a directory into one set. Requires a directory, a list of files to load and the names of the variables in the set
create_count_set <- function(directory, file_list, run_names){
  for (i in seq_along(file_list)) {
    print(paste0("Loading file ",file_list[i]))
    assign(run_names[i],
           read.csv(file.path(directory, file_list[i]), header = FALSE) %>%
             setNames(., c("seq", run_names[i])))
    if (i == 1){
      merged_set = get(run_names[1])
      print("Loading first file")
    }
    else {
      print(paste0("Merging file ",file_list[i]))
      merged_set = full_join(merged_set, get(run_names[i]), by = "seq")
    }
    rm(list = run_names[i])
  }

  rm(run_names)
  #Replace NA by zeros
  merged_set <- merged_set %>%
    mutate_at(c(2:ncol(merged_set)),
              ~replace_na(.,0))
  
  print("Removed NAs")
  
  # Add a peptide_seq that contains the translation
  merged_set <- merged_set %>%
    mutate(peptide_seq = seq %>%
             dna() %>%
             seq_translate() %>%
             as.character()) %>%
    relocate(peptide_seq, .after = seq)
  
  print("Added translation.")
  gc()
  return(merged_set)
}

# remove_stop_codons()
##Create a function that drop all the peptides with stop codons from a data frame
remove_stop_codons <- function(peptide_data_frame){
  output_df <- peptide_data_frame[!grepl("\\*", peptide_data_frame$peptide_seq),]
  return(output_df)
}
