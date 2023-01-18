## This script aims to take counts of peptides and do analysis of them.

library(tidyverse)  #Needed for data wrangling
library(DESeq2)  #Use for differential expression

# This folder should contain the CSV files with a column for sequence and a column for counts

directory <- "C:/Users/worms/NGS Data/2022.06.07_drift_seq/90-666155004b/00_fastq/NNK/NNK3/counts_csv"

#directory = choose.dir()

# Get the count data

file_list <- list.files(directory)

#Keep only the data csv, trim the name a bit and then import all the csvs.

file_list <- file_list[grepl("Gen", file_list)]

run_names <- file_list %>%
  gsub("Cytoplasmic-NNK-", "", .) %>%
  gsub("_001_peptide3_count.csv", "", .) %>%
  gsub("-", "_", .) %>%
  tolower()


for (i in seq_along(file_list)) {
  assign(run_names[i],
         read.csv(file.path(directory, file_list[i]), header = FALSE) %>%
                        setNames(., c("seq", run_names[i]))
  )
}

## Merge by sequence to form one big dataset.

merged_set = get(run_names[1])

for (i in 2:length(file_list)) {
  merged_set = full_join(merged_set, get(run_names[i]), by = "seq")
}
rm(list = run_names)

#Gotta get all my columns lowercase. Down with capital!
 
# Replace all NAs by 0

merged_set <- merged_set %>%
  mutate_at(c(2:ncol(merged_set)),
           ~replace_na(.,0))

# Going to split the set for induced and repressed conditions. 

induced_set <- merged_set %>% 
  select(!contains("glu"))

repressed_set <- merged_set %>%
  select(!contains("ara"))
  


# Clean the merged_set from stop codons. As of jan 23, I'm not doing that

# dropThose = merged_set$seq %>% grepl('\\*' , . )
# merged_set <- merged_set[!dropThose,]

## Compute sum of counts for each condition and then compute ratios

merged_set <- merged_set %>% 
  mutate(Gen_1_sum = Gen_1_LB_R1 + Gen_1_LB_R2) %>%
  mutate(Gen_5_Ara_sum = Gen_5_Ara_R1 + Gen_5_Ara_R2) %>%
  mutate(Gen_5_Glu_sum = Gen_5_Glu_R1 + Gen_5_Glu_R2) 

Gen1_read_total <- sum(merged_set$Gen_1_sum)
Gen5_ara_total <- sum(merged_set$Gen_5_Ara_sum) 
Gen5_glu_total <- sum(merged_set$Gen_5_Glu_sum)

merged_set <- merged_set %>%
  mutate(Gen_1_ratio = Gen_1_sum / Gen1_read_total) %>%
  mutate(Gen_5_ara_ratio = Gen_5_Ara_sum / Gen5_ara_total) %>%
  mutate(Gen_5_glu_ratio = Gen_5_Glu_sum / Gen5_glu_total)

## Compute enrichment ratio by comparing with Gen 1

merged_set <- merged_set %>%
  mutate(enrichment_ara = log2(Gen_5_ara_ratio/Gen_1_ratio)) %>%
  mutate(enrichment_glu = log2(Gen_5_glu_ratio/Gen_1_ratio))

#Have an issue with zero ratio for Gen1, need to add a tiiiiiiiny amount to all I think. Or actually just
# refuse to compute ratio for those? With NAs or something

#So the proper way to do it is to split the dataset into two and run the analysis separately for now.





#How many peptides don't have any zeroes?
x <- list()
nRow <- seq(nrow(merged_set))
nCol <-  8:10

for (i in nRow){
  y = FALSE
  for (j in nCol){
    if (merged_set[[j]][i] == 0){
      y = TRUE
    }
  }
  x <- append(x, y)  
  
  print(i)
  }
    


hasZeroes <- merged_set %>%
  select(8:10) %>%
  rowSums(merged_set == 0)


rm(hasZeroes)  

oriSet <- merged_set %>%
  select(8:10)

x <- rowSums(filterSet == 0)

filterSet <- merged_set[!rowSums(oriSet == 0)]
