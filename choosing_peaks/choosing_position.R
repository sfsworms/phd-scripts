library(tidyverse)

#This code will look at positions where I have single nucleotides

path <- file.choose()

peptide <- read.csv(path)

output_df <- data.frame(pep_name = character(),
                        pos = integer(),
                        seq = character(),
                        nt= character()
                        )

for(i in seq_len(nchar(peptide$seq[1]))){
  nucleotide <- peptide %>% pull(seq) %>% substring(i,i)
  freq <- table(unlist(nucleotide))
  unique_nucleotide <- names(freq)[freq == 1]
  
  for(j in unique_nucleotide){
    new_line <- data.frame(pep_name = peptide[j == nucleotide,] %>% pull(pep),
                           pos = i,
                           seq = peptide[j == nucleotide,] %>% pull(seq),
                           nt = j
                           )
    output_df <- rbind(output_df, new_line)
  }
}

if(!(output_df %>% pull(pep_name) %>% unique() %>% length()) == nrow(peptide)){
  print("Not all peptides have a position where they're uniques.")
  break()
  }

write.csv(output_df,
          file = file.path(dirname(path), "suggested_positions.csv"))
