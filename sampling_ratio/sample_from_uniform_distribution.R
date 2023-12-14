source("functions.R")

# Making an initial, uniformly distributed library

clone_number <- 40
sampling_ratio <- 10

uniform_lib <- data.frame(clone_id = c(1:clone_number),
                          clone_count = rep(100, clone_number))

sampled_lib <- subset_library(uniform_lib, clone_number*sampling_ratio, replacement = TRUE)

# Plot the distributions
df <- data.frame(rbind(uniform_lib %>% 
                         mutate(type = "original"),
                       sampled_lib %>%
                         mutate(type = "sampled")
                       )
                 ) %>%
  group_by(type) %>%
  mutate(ratio = clone_count/sum(clone_count)) %>%
  ungroup()


plot <- ggplot(df, aes(x = clone_id, y = ratio, fill = type)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Clone",
       y = "Proportion",
       fill = "Library") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_discrete(labels = c("Original library", "Post sampling library")) +
  theme_bw()

plot

ggsave(filename = file.path(getwd(),
                            "sample_from_uniform.png"),
       plot = plot,
       width = 12,
       height = 8,
       units = 'cm')
