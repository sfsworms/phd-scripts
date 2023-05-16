# Set the directory path where the data is located
directory <- "C:/Users/worms/ngs_data/2022_06_07_drift_seq/90-666155004b/00_fastq/all_files"

# Load required packages
library(tidyr)
library(DESeq2)

# Read the count data CSV file located in the specified directory
count_set <- read.csv2(file = file.path(directory, "count_set_long_pseudocount.csv"))

# Filter the count_set dataframe to retain only rows where 'library' is 'nnk' and 'condition' is 'induced'
count_set <- count_set %>%
  filter(library == "nnk" & condition == "induced")

# Replace all infinite 'enrichment_ratio' values with NA
count_set <- count_set %>%
  mutate(enrichment_ratio = ifelse(is.infinite(enrichment_ratio),NA,enrichment_ratio)) 

# Filter the count_set dataframe to retain only rows where 'seq' matches the specific regular expression pattern
count_set <- count_set %>%
  filter(grepl("^TGC([ATGC][ATGC][GT]){7}$", seq)) 

# Select the first, sixth, and seventh columns of the count_set dataframe
count_set_deseq <- count_set %>% select(1,6,7)

# Convert the selected dataframe into a standard dataframe
count_set_deseq <- count_set_deseq %>% data.frame()

# Assign row names to the dataframe based on the 'seq' column
row.names(count_set_deseq) <- count_set_deseq %>% pull(seq)

# Select the second and third columns of the dataframe
count_set_deseq <- count_set_deseq %>% select(c(2,3))

# Convert the 'gen1' and 'gen5' columns to integer type
count_set_deseq <- count_set_deseq %>% 
  mutate(gen1 = as.integer(gen1), gen5 = as.integer(gen5))

# Create a new dataframe for condition information
condition_df <- data.frame(gen = factor(c(1,5)),
                           row.names = c("gen1", "gen5"))

# Construct a DESeq dataset object from the matrix of count data and the condition dataframe, with the experimental design formula ~ gen
dds <- DESeqDataSetFromMatrix(countData = count_set_deseq, colData = condition_df, design = ~ gen)

# Run the DESeq2 differential expression pipeline on the dataset
dds <- DESeq(dds)

# Extract the results from the DESeq dataset
res <- results(dds)
