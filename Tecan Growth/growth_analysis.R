## This link to the function to extract the data from the tecan

source(file = "C:/Users/worms/Dropbox/PhD/PhD-Scripts/Tecan Growth/utils.R")

folder.loc <- "C:/Users/worms/Dropbox/PhD/Data/Tecan/2022.05.27 Test peptide marker"

## This extract the data from the tecan file, add well names to the well ID and add details contained in 'experiment.xlsx'

data_frame <- get_data_tecan(folder_loc = folder.loc) %>%
  add_well_name(., folder_loc = folder.loc) %>%
  add_experiment_details(., folder_loc = folder.loc)


# Now I can graph things.

