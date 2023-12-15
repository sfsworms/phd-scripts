library(ggplot2)
library(tibble)
library(Peptides)
library(scatterplot3d)
library(stringr)

# Read in the count_set
directory <- "C:/Users/worms/ngs_data/2022_06_07_drift_seq/90-666155004b/00_fastq/all_files"

# Read a sample of the set
count_set <- read.csv2(file = file.path(directory, "count_set_params.csv"))  %>%
  filter(library == "nnk" & condition == "induced") %>% # Take only NNK, induced peptides
  mutate(enrichment_ratio = ifelse(is.infinite(enrichment_ratio),NA,enrichment_ratio)) %>% # Going to drop all the infinite enrichment ratio to NA
  filter(grepl("^TGC([ATGC][ATGC][GT]){7}$", seq)) %>%
  filter(!grepl("\\*", standard_seq)) #Remove sequences with stop codons

# Subsetting the count_set for sanity reasons
set.seed(123)
#count_set_sample <- count_set[sample(c(1:nrow(count_set)), 10^4),] # Subsett some rows
# Subset count with higher homonym_seq

count_set_sample <- count_set %>% filter(homonym_seq >= 4)

## Add insta index, aliphatic moment, RR counts, H counts, P counts
count_set_sample <- count_set_sample %>%
  mutate(insta_index = instaIndex(standard_seq)) %>%
  mutate(hydrophobic_moment = hmoment(standard_seq)) %>%
  mutate(rr_count = str_count(standard_seq, pattern = "RR")) %>%
  mutate(h_count = str_count(standard_seq, pattern = "H")) %>%
  mutate(P_count = str_count(standard_seq, pattern = "P")) 

# Add a category variable (enriched, neutral, depleted)
enrichment_cutoff <- 6

count_set_sample <- count_set_sample %>% mutate(category = ifelse(enrichment_ratio_log < -enrichment_cutoff, "depleted",
                                              ifelse(enrichment_ratio_log < enrichment_cutoff, "neutral", "enriched")))

count_set_sample <- count_set_sample %>% select(-X) 

count_set_sample <- count_set_sample %>% filter(category == "depleted")

row.names(count_set_sample) <- NULL # Remove the name introduced by the sample
count_set_sample <- column_to_rownames(count_set_sample, var = "seq") 

# Rename rows with sequence for the PCA and select the variables for PCA

count_set_pca <- count_set_sample %>% # Rename the rows of the set with seq (to use as identified)
  select(c(hydrophob,charge,aindex,boman, insta_index, hydrophobic_moment, rr_count, h_count, P_count))

#Scaling the variables I'll use prior to PCA
count_set_pca_scaled <- scale(count_set_pca)

# Perform PCA
pca_result <- prcomp(count_set_pca_scaled)

# Print summary of the PCA result
print(summary(pca_result))

# Show the importance of components
screeplot(pca_result, type="lines")

## Okay, gotta keep PC 1 & 2

# Visualizing the PCA result
biplot(pca_result)

# Convert PCA results to a data frame
pca_df <- data.frame(sample = rownames(pca_result$x), pca_result$x)

# Assume 'grouping_variable' is your external grouping variable
pca_df$category <- count_set_sample$category[match(pca_df$sample, rownames(count_set_sample))]

ggplot(pca_df %>% filter(category != "neutral"), aes(x = PC1, y = PC2, color = category)) +
  geom_point(alpha = 0.5) +
  theme_minimal()

# Try to do k-mean clustering with 3 clusters, using the first 2 PCA component
kmeans_result <- kmeans(pca_result$x[,1:2], centers = 3)

# Create a data frame with the PCA results and the cluster assignment
pca_df <- pca_df %>% mutate(cluster = factor(kmeans_result$cluster))

# Create the plot
ggplot(pca_df, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point(alpha = 0.4) +
  theme_minimal()

# Selecting sequence

pca_df %>% 
  filter(cluster == 3) 
