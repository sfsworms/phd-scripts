# Load libraries and functions
source("functions.R")


# Simulate a library of 100 different clones
X <- 1e2
original_library <- simulate_library(X)
#original_library[1,2] <- 10000
plot(original_library)

# Subset 200 bacteria from the original library
Y <- 10000
subsetted_library <- subset_library(original_library, Y)
plot(subsetted_library)


# Plot the distributions
df <- data.frame(rbind(original_library %>% 
                         mutate(type = "original"),
                        subsetted_library %>%
                          mutate(type = "subsetted")
                       )
                 )

ggplot(df, aes(x = clone_id, y = clone_count, fill = type)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Distribution of Bacterial Clones",
       x = "Clone ID",
       y = "Count")

# Add normalized value
df <- df %>% 
  group_by(type) %>%
  mutate(ratio = clone_count/sum(clone_count)) %>%
  ungroup() %>% 
  pivot_wider(names_from = type, values_from = c(ratio,clone_count))

# Plot normalized values against each other
ggplot(df, aes(x = ratio_original, y = ratio_subsetted)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +  # Optional: Add a linear fit
  labs(title = "Comparison of Clone Ratio: Original vs Subsetted",
       x = "Original Library Counts",
       y = "Subsetted Library Counts")

# Compute correlation coefficient

cor(df$clone_count_original, df$clone_count_subsetted)

# Do it at the function level
subsetting_correlation <- function(sampling_ratio, X = 100){
  ori_lib <- simulate_library(X)
  sub_lib <- subset_library(ori_lib, X * sampling_ratio)
  correlation <- cor(ori_lib$clone_count, sub_lib$clone_count)
  return(correlation)
  
}

ratio = c(1:10)

cor_df <- data.frame(sampling_ratio = c(1:200)) %>%
  mutate(correlation = sapply(sampling_ratio, subsetting_correlation))

ggplot(cor_df, aes(x = sampling_ratio, y = correlation)) +
  geom_point() +
  labs(title = "Correlation of ratio after subsetting depending on sampling ratios",
       x = "Sampling ratios",
       y = "Correlation coefficient") + 
  geom_smooth() +
  theme_bw()
## Add line at 10, 50
## Add trendline
