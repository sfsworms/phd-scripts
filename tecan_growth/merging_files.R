# Merging the datasets
 ## This function is used to merge datasets representing multiple repetitions of the same experiment.

library(tidyverse) # tidyverse

script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)# Get the path to this script.

source(file = file.path(script_dir, "utils.R")) # Load the functions for analysis of the tecan data

folders_names <- c("2023.06.28_test_select_ara", "2023.06.29_test_select_ara" )
output_name <- "merged_ara.csv"

directory <- "C:/Users/worms/Dropbox/data/tecan"

tecan_data <- data.frame()

for(i in seq_along(folders_names)){
  folder_loc <- file.path(directory,folders_names[i])
  new_data <- get_data_tecan(folder_loc = folder_loc) %>% # Load the tecan data
    add_well_name(., folder_loc = folder_loc) %>%
    add_experiment_details(., folder_loc = folder_loc) %>%
    mutate(experiment = folders_names[i]) %>%
    mutate(experiment_number = i)
  
  tecan_data <- rbind(tecan_data, new_data)
}

write.csv(tecan_data, file = file.path(directory, output_name))
