# Load required libraries
library(ggplot2)
library(tidyverse)


# Function to simulate bacterial library using a poisson distribution.

simulate_library <- function(X, mean_count = 10) {
  data.frame(clone_id = c(1:X),
             clone_count = rpois(X, mean_count)*100) %>%
    return()
}

# Function to subset bacterial library
subset_library <- function(lib, Y, replacement = FALSE) {
  if (sum(lib$clone_count) < Y) stop("Not enough bacteria to subset")
  count_array <- lib %>% pull(clone_count)
  
  new_lib <- rep(0, nrow(lib))
  
  for (i in 1:Y) {
    # Randomly select a clone based on its proportion in the library
    selected_clone <- sample(1:length(count_array), size = 1, prob = count_array / sum(count_array))
    new_lib[selected_clone] <- new_lib[selected_clone] + 1
    if (!replacement) lib[selected_clone,2] <- lib[selected_clone,2] - 1
  }
  
  return(data.frame(clone_id = lib$clone_id,
                    clone_count = new_lib))
}

