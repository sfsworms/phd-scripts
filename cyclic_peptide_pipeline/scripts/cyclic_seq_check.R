setwd("C:/Users/worms/phd-scripts/cyclic_peptide_pipeline")

source("functions/function_analysis.R")

# Get a sorted list of the 20 AA, and add a "*" for stop codons
sorted_AA <- c(sort(AA_STANDARD),"*")

# Let's create a function creating the "standard sequence"

standard_sequence <- function(aa_seq = "CABCDC"){
  # If there is one C, return its sequence
  if(str_count(aa_seq,"C") == 1) return(aa_seq)
  
  # If there isn't a C, this is a mistake
  if(str_count(aa_seq,"C") == 0) return(NA)
  
  else{
  # Get the positions of all the Cs
    positions <- str_locate_all(aa_seq,"C")[[1]][,1] 
    possible_sequence <- vector()
    
    for(i in positions){
      rotated_seq <- paste0(substr(aa_seq, i, nchar(aa_seq)),
                            substr(aa_seq, 0, i-1))
      
      possible_sequence <- c(possible_sequence, rotated_seq)

    }
  possible_sequence %>% max() %>% return()
  }
}

# Let's add that to my short peptide merged set.

short_pep_set <- read.csv2(file.choose())

short_pep_set_test <- short_pep_set %>%
  mutate(std_seq = standard_sequence(peptide_seq))

