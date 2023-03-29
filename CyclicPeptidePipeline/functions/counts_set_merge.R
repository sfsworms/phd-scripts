# Goal: Merge the NNB3 and NNK3 files

library(tidyverse)

file_NNK <- file.choose()
file_NNB <- file.choose()

NNB <- read.csv(file_NNB)
NNK <- read.csv2(file_NNK)

NNB <- NNB %>% #There is that column of x
  select(-1)

NNB <- NNB %>%
  mutate(library = 'NNB')

NNB <-rename(NNB, condition = induction)

NNK <- NNK %>%
  mutate(library = 'NNK')

all_counts <- rbind(NNB,NNK) # Merge the data sets, then randomize them to make sampling later easier.

all_counts <- all_counts[all_counts %>% nrow() %>% sample(),]

destination <- dirname(file_NNK) %>% dirname() %>% dirname()

write.csv2(all_counts, 
           file = file.path(destination,"all_NNK3.csv"),
           row.names = FALSE,
           append = FALSE)

sample_output <- read.csv2(file.path(destination,"all_NNK3.csv"), #Just checking it looks ok.
                           nrows = 1000)

## Going to just creat my correct sequence dataset

all_correct_counts <- all_counts %>%
  filter(grepl("^TGC.*", seq))

write.csv2(all_correct_counts, 
           file = file.path(destination,"all_correct_NNK3.csv"),
           row.names = FALSE,
           append = FALSE)
