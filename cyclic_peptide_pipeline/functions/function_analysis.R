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
## Assumed a merged count set where each row has counts and is named library_generation_other

#' compute_ratios
#' compute an enrichment ratio based on the final ratio over the initial one, and removes NAs caused by divide by zero error
#' This function assumes a count set where each columns has counts. still need to make the enrichment work if the columns aren't name gen1 and gen5
#' 
#' @param df The count_set with the counts column to be ratio'ed. 
#' @param count.cols Columns containing the counts to be ratio'ed
#' @param compute_enrichment Should an enrichment be computed. Default to TRUE, remove if not applicable. Only works if you have two columns in count.cols()
#' @param info.cols Columns containing information by which to split the count_set if columns contains counts from multiple sequencing
#'
#' @return
#' @export
#'
#' @examples
compute_ratios <- function(df, count.cols = c(5:6), info.cols = c(3:4), compute_enrichment = TRUE){
  # Get all the colnames and decide new colnames
  col_names <- colnames(df)
  new_col_names <- paste0(col_names[count.cols],"_ratio")
  
  df <- df %>% 
      mutate(gen1 = gen1 + 0.1, gen5=gen5+0.1) %>% #adds a pseudocount
       group_by(across(all_of(col_names[info.cols]))) %>% # Group them by the columns used 
      mutate(across(all_of(col_names[count.cols]), function(x){x/sum(x)}, .names = "{.col}_ratio")) %>% # Calculate a ratio for each column specified in count.cols
      ungroup() # Remove the grouping
      
  # Compute a gen5 over gen1 ratio if the option is set
  
  if(compute_enrichment){
    if(length(count.cols) != 2){
      break
    }
    
    df <- df %>%
      mutate(enrichment_ratio = gen5_ratio/gen1_ratio) %>%
      mutate(enrichment_ratio_log = log2(enrichment_ratio))
  }
    
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
## Create a function that drop all the peptides with stop codons from a data frame 

remove_stop_codons <- function(peptide_data_frame){
  output_df <- peptide_data_frame[!grepl("\\*", peptide_data_frame$peptide_seq),]
  return(output_df)
}

## This function goes through a character list representing cyclic peptides. It looks at all possible Cs and defined a "standard" by
## comparing ASCII value so that the sequence starting with C with the highest value is taken

## Note on 09/05/2023 this was renamed standard_sequence from standard_sequence2, will need to rename in the functions.

#' Standard sequence
#' This function goes through a character list representing cyclic peptides. It looks at all possible Cs and defined a "standard" by
#' Note on 09/05/2023 this was renamed standard_sequence from standard_sequence2, will need to rename in the functions.
#' @param aa_seq 
#'
#' @return A character vector containing the sequence starting with C with the highest ASCII value.
#' @export
#'
#' @examples

standard_sequence <- function(aa_seq) {
  
  # Calculate the number of C's in each sequence
  num_c <- str_count(aa_seq, "C")
  
  # Find sequences with only one C
  single_c_idx <- num_c == 1
  single_c_seq <- aa_seq[single_c_idx]
  
  # Find sequences with more than one C
  multi_c_idx <- num_c > 1
  multi_c_seq <- aa_seq[multi_c_idx]
  
  # Extract positions of C's in multi-C sequences
  multi_c_pos <- str_locate_all(multi_c_seq, "C")
  
  max_seq <- vector()
  # Generate all possible sequences for multi-C sequences and select the best
  for(j in seq_along(multi_c_pos)){
    max_seq[j] <- sapply(multi_c_pos[[j]][,1], function(i) {
      paste0(substring(multi_c_seq[j], i, nchar(multi_c_seq[j])),
             substring(multi_c_seq[j], 1, i-1)) 
    }) %>% max()
    if(j %% 1000 ==0){print(paste0("Treated multi-c sequences: ",j))}
  }
  
  # Combine the sequences with one C and the chosen sequences for multi-C sequences
  standard_seq <- rep(NA, length(aa_seq))
  standard_seq[single_c_idx] <- single_c_seq
  standard_seq[multi_c_idx] <- max_seq
  
  return(standard_seq)
}

# The below function takes a count_set. For each peptide sequence produced (including circular homonyms), it looks at all the pairwise
# combination of gene producing those peptides, and look at the correlation coefficient of enrichment ratios. 
# As of now, it just ignore infinite enrichment rations (where no reads are seen in one of the generation)

get_enrichment_list <- function(df){
  ratio_df <- df %>%
    filter(!is.na(enrichment_ratio_log)) %>%
    split(., .$standard_seq) %>%
    Filter(function(x) nrow(x) >1, .) %>%
    lapply(FUN= function(df) df$enrichment_ratio_log) %>%
    lapply(FUN= function(list) t(combn(list, m=2))) %>%
    lapply(FUN = data.frame) %>%
    Reduce(rbind,.)
  
  return(ratio_df)
}


