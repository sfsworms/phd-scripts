# Script to get extract and sort data from Tecan excel files

library(tidyverse)  #General tidy stuff
library(readxl)  #Use to read the excel file
library(openxlsx)  #USed to write the tidy data to
library(platetools)  #Needed to play with plates a bit more easily


get_data_tecan <- function(folder_loc = folder.loc, tecan_file_name = "tecan.xlsx") {
  
  # This does the data import. Note that it assumes that the data file is name
  # tecan.xlsx. If that is not the case edit the code below
  
  tecanDataLoc <- paste(folder_loc, tecan_file_name, sep = "/")
  
  tecan_data <- read_excel(tecanDataLoc)
  
  # This detect how much fluff there is in front of the data and then re-import
  # dropping that fluff.  Works for TECAN file because they have a 'Cycle Nr.'
  # column. It then detect the last cycle and cut the bottom off.
  
  for (i in 1:dim(tecan_data)[1]) {
    if (isTRUE((tecan_data[i, 1] == "Cycle Nr.")[1])) {
      break
    } else {
      i <- i + 1
    }
  }
  tecan_data <- read_excel(tecanDataLoc, skip = i)
  rm(i)
  
  names(tecan_data)[1] <- "cycle"
  names(tecan_data)[2] <- "time"
  tecan_data <- tecan_data %>%
    select(-3) %>%
    mutate(cycle = as.numeric(cycle)) %>%
    mutate(time = as.numeric(time))
  
  
  for (i in 1:nrow(tecan_data)) {
    if (tecan_data[i, 1] %>%
        is.na()) {
      break
    }
  }
  tecan_data <- tecan_data[1:i - 1, ]
  rm(i)
  
  ## This convert the data to to a 'long' format with the tidyverse package.
  tecan_data <- gather(tecan_data, wellID, OD, -c(cycle, time), factor_key = TRUE)
  
  return(tecan_data)
}



## This create a data frame containing the well numbers (wellID) as well as the well names (wellName) from a
## plate.xlsx file in folder_loc.It assumes that well marked with 'x' are empty and discard them.  It then merge
## that info with tecan_data

add_well_name <- function(data = tecan_data, folder_loc = folder.loc){
  
  well_name <- data.frame(wellID = num_to_well(1:96, plate = 96) %>%
                            sub("^([A-H])0+", "\\1", .), 
                          wellName = paste(folder_loc, "plate.xlsx", sep = "/") %>%
                            read.xlsx(., colNames = FALSE) %>%
                            as.matrix() %>%
                            t() %>%
                            as.list() %>%
                            unlist())
  
  well_name <- well_name %>% # 'X' indicate an empty well
    filter(toupper(wellName) != "X")
  
  data <- inner_join(data, well_name, by = "wellID")
  
  return(data)
}

## This function just go and get the experimental details from an experiment.xlsx file

add_experiment_details <- function(data= tecan_data, folder_loc = folder.loc){
  experiment.details <- paste(folder_loc, "experiment.xlsx", sep="/") %>%
    read.xlsx()
  
  data <- left_join(data, experiment.details, by = "wellName")
  return(data)
}

blanking_tecan_data <- function(tecan_data){
  # Check if the data frame is a data frame
  if(!is.data.frame(tecan_data)){
    stop("The argument of the function needs to be a data frame.")
  }
  
  # First I need to check that the blank wells aren't growing
  blank_wells <- tecan_data %>%
    filter(wellName == "blank")
  
  # If there are no blanck wells, throw an error
  if(!(blank_wells %>% nrow())){
    stop("No blank well found. Are you sure they are labelled 'blank'")
  }
  
  # # Check that the blank wells aren't growing. I put a cutoff at 20%
  # min_blank <- blank_wells %>% pull(OD) %>% min()
  # max_blank <- blank_wells %>% pull(OD) %>% max()
  # if((max_blank/min_blank) > 1.2){
  #   stop("The blank wells are growing more than 20% over the time frame.")
  #}
  
  # I'm going to average the blanks, then substract the appropriate blank from each of the OD.

  # Keep only the relevant data for blank wells
  blank_wells <- blank_wells %>%
    select(c("cycle", "time", "wellID", "OD"))
  
  # Shift to a wider format
  blank_wells <- blank_wells %>%
    pivot_wider(names_from = wellID, values_from = OD)
  
  # Compute the average OD and drop the unneeded columns
  
  blank_wells <- blank_wells %>%
    rowwise() %>%
    mutate(avg_blank = mean(c_across(3:ncol(blank_wells)))) %>%
    select(cycle, avg_blank)
  
  # Add the column to the tecan_data, correct OD by substracting blank, drop the blank wells and column.
  
  tecan_data <- full_join(tecan_data, blank_wells, by = "cycle") 

  tecan_data <- tecan_data %>%
    mutate(OD = OD - avg_blank) %>%
    select(-avg_blank) %>%
    filter(wellName != "blank")
  
  
  return(tecan_data)
}