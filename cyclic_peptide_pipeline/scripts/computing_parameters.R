# The goal of this script is to sort the mess of my analysis notebook by separating the computations from the discussions

directory <- "C:/Users/worms/ngs_data/2022_06_07_drift_seq/90-666155004b/00_fastq/all_files"

library(Peptides)
library(tidyverse)

count_set <- read.csv2(file = file.path(directory, "count_set_long.csv"))  %>%
  filter(library == "nnk" & condition == "induced") %>% # Take only NNK, induced peptides
  mutate(enrichment_ratio = ifelse(is.infinite(enrichment_ratio),NA,enrichment_ratio)) %>% # Going to drop all the infinite enrichment ratio to NA
  filter(grepl("^TGC([ATGC][ATGC][GT]){7}$", seq)) 

count_set <- count_set %>%
  group_by(standard_seq) %>%
  mutate(homonym_seq = n()) %>%
  mutate(average_enrich_ratio = psych::geometric.mean(enrichment_ratio)) %>%
  ungroup()

count_set <- count_set %>%
  mutate(hydrophob = hydrophobicity(standard_seq)) %>%
  mutate(charge = charge(standard_seq)) %>%
  mutate(pi = pI(standard_seq)) %>%
  mutate(aindex = aIndex(standard_seq)) %>%
  mutate(boman = boman(standard_seq)) %>%
  mutate(RRcount = str_count(standard_seq, "RR")) %>%
  mutate(RRcount = as.factor(Rcount))

write.csv(count_set, file = file.path(directory,"count_set_params2.csv"), row.names = FALSE)