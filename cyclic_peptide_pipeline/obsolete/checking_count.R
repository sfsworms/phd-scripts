library(tidyverse)

path <- file.path("C:/Users/worms/ngs_data/2022_06_07_drift_seq/90-666155004b/00_fastq/all_files/peptide_count_csv")

file_list <- list.files(path)
file_list <- file_list[grepl("_count.csv", file_list)]

df <- data.frame(file = character(),
                 sequence = integer(),
                 count = integer())

for(i in seq_along(file_list)){
  x <- file.path(path,file_list[i]) %>%
    read.csv(header = FALSE)
  new_line <- data.frame(file = file_list[i],
                         sequence = nrow(x),
                         count = x[,2] %>% sum())
  print(new_line)
  df <- rbind(df, new_line)
}

write.csv2(df, file = file.path(path,"count_check.csv"), row.names = FALSE)

## Checking counts in the merged set
path_merged_set_short <- file.choose()

merged_set_short <- read.csv2(path_merged_set_short)

merged_set_short %>% 
  filter(library == "nnk") %>%
  filter(gen1_lb_12 > 0) %>% 
  pull(4)  %>%
  sum()

# The above code did it manually, let's automate that a bit for the filtered set

# Only select the ones with an NNK sequence
merged_set_short <- merged_set_short %>%
  filter(grepl("^TGC([ATGC][ATGC][GT]){3}$", seq))

out_df <- data.frame(library = character(),
                     gen = character(),
                     counts = integer(),
                     sequences = integer())

for(i in 4:6){
  for(j in c("nnb","nnk")){
    newline <- data.frame(library = j,
                          gen = (merged_set_short %>% colnames)[i],
                          counts = merged_set_short %>% filter(library == j) %>% pull(i) %>% sum(),
                          sequences = (((merged_set_short %>% filter(library == j) %>% pull(i)) > 0 )%>% sum())
    )
    out_df <- rbind(out_df, newline)
  }
}

write.csv(out_df, file = file.path(dirname(path_merged_set_long),
                                   "counts_merged_short_TGC.csv"),
          row.names = FALSE)

# Now for the long set

path_merged_set_long <- file.choose()

merged_set_long <- read.csv2(path_merged_set_long)

colnames <- colnames(merged_set_long)

merged_set_long %>% 
  filter(nnb_gen1_lb_24 > 0) %>% 
  nrow()

test_df <- data.frame(column = character(),
                      sequences = integer(),
                      counts = integer())

for(i in 4:9){
  column <- colnames[i]
  sequences <- ((merged_set_long %>% pull(i)) > 0 )%>% sum()
  counts <- merged_set_long %>% pull(i) %>% sum()
  new_line <- data.frame(column, sequences, counts)
  test_df <- rbind(test_df, new_line)
}

write.csv2(test_df, file = file.path(dirname(path_merged_set_long),
                                     "counts_long_merged_set.csv"))

merged_set_long <- merged_set_long %>%
  filter(grepl("^TGC([ATGC][ATGC][GT]){7}$", seq))

for(i in 4:9){
  column <- colnames[i]
  sequences <- ((merged_set_long %>% pull(i)) > 0 )%>% sum()
  counts <- merged_set_long %>% pull(i) %>% sum()
  new_line <- data.frame(column, sequences, counts)
  test_df <- rbind(test_df, new_line)
}

write.csv2(test_df, 
           file = file.path(dirname(path_merged_set_long),
                                     "counts_long_merged_set_TGC.csv"),
           row.names = FALSE)

