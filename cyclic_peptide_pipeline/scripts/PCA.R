library(ggplot2)
library(tibble)

# Subsetting the count_set for sanity reasons
set.seed(123)
count_set_sample <- count_set[sample(c(1:nrow(count_set)), 10^4),] # Subsett 10^5 rows


# Create a "enrich" categories to split enriched, purged and neutral-ish seq

row.names(count_set_sample) <- NULL # Remove the name introduced by the sample

count_set_pca <- column_to_rownames(count_set_sample, var = "seq") %>% # Rename the rows of the set with seq (to use as identified)
  select(c(hydrophob,charge,pi,aindex,boman))

#Scaling the variables I'll use prior to PCA
count_set_pca_scaled <- scale(count_set_pca)

# Perform PCA
pca_result <- prcomp(count_set_pca_scaled)

# Print summary of the PCA result
print(summary(pca_result))

# Show the importance of components
screeplot(pca_result, type="lines")

# Visualizing the PCA result
biplot(pca_result)

# Convert PCA results to a data frame
pca_df <- data.frame(Sample = rownames(pca_result$x), pca_result$x)

# Assume 'grouping_variable' is your external grouping variable
pca_df$grouping_variable <- your_data$grouping_variable[match(rownames(pca_result$x), rownames(your_data))]

ggplot(pca_df, aes(x = PC1, y = PC2, color = grouping_variable)) +
  geom_point() +
  theme_minimal()
