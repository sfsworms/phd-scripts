path <- file.path("C:/Users/worms/ngs_data/2022_06_07_drift_seq/90-666155004b/00_fastq/all_files/peptide_count_csv")

file_list <- list.files(path)
file_list <- file_list[grepl("_count.csv", file_list)]

df <- data.frame(file = character(),
                 count = numeric())

for(i in seq_along(file_list)){
  x <- file.path(path,file_list[i]) %>%
    read.csv(header = FALSE)
  new_line <- data.frame(file = file_list[i],
                         count = x[,2] %>% sum())
  print(new_line)
  rbind(df, new_line)
}


count_set <- read.csv(path)
