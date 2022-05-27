#Script to get mumax values for Gol

library(tidyverse)
library(growthrates)
library(growthcurver)
library(lattice)

library(readxl) #Use to read the excel file
library(openxlsx) #USed to write the tidy data to excel

get_data_tecan <- function(tecan_file_name="tecan.xlsx", folder_loc = "C:/Users/worms/Desktop/Gol_ Data for R analysis/TECAN"){
  
  #This does the data import. Note that it assumes that the data file is name tecan.xlsx. If that is not the case edit the code below
  
  tecanDataLoc <- paste(folder_loc,tecan_file_name, sep="/")
  
  growth_data <- read_excel(tecanDataLoc)
  
  #This detect how much fluff there is in front of the data and then re-import dropping that fluff. 
  #Works for TECAN file because they have a "Cycle Nr." column. It then detect the last cycle and cut the bottom off.
  
  for (i in 1:dim(growth_data)[1]){
    if(isTRUE((growth_data[i,1] == "Cycle Nr.")[1])){
      break
    }
    else{
      i <- i+1
    }
  }
  
  growth_data <- read_excel(tecanDataLoc, skip = i)
  rm(i)
  
  #Rename first two columns;, convert them to num and drop temperature to make things easier and trim the datestamp
  
  names(growth_data)[1] <-"cycle"
  names(growth_data)[2] <-"time"
  growth_data <- growth_data %>%
    select(-3) %>% 
    mutate(cycle = as.numeric(cycle)) %>%
    mutate(time = as.numeric(time))
  
  
  for (i in 1:dim(growth_data)[1]){
    if(isTRUE((growth_data[i,1] > 0)[1])){
      i <- i+1
    }
    else{
      break
    }
  }
  growth_data <- growth_data[1:i-1,]
  rm(i)
  
  #This convert the data to to a "long" format with the tidyverse package. 
  growth_data <- gather(growth_data, wellID, OD,-c(cycle,time), factor_key=TRUE)
  
  return(growth_data)
}

get_data_hidex <- function(tecan_file_name=tecan.xlsx, folder_loc = "C:/Users/worms/Desktop/Gol_ Data for R analysis/Hidex"){
 
  hidexDataLoc <- paste(folder_loc,tecan_file_name, sep="/")
  
  #Drop the crap before the data
  growth_data <- read_excel(hidexDataLoc, sheet = "Raw OD(600)")
  
  for (i in 1:dim(growth_data)[1]){
    if(isTRUE((growth_data[i,2] == "Well")[1])){
      break
    }
    else{
      i <- i+1
    }
  }
  
  growth_data <- read_excel(hidexDataLoc, sheet = "Raw OD(600)", skip = i, col_names = T)
  rm(i)
  
  #Drop all the times but the first one and the well column then reorganize
  
  growth_data <- growth_data %>%
    select(-1, -3) %>%
    pivot_longer(!Well, names_to = "time", values_to = "OD") %>%
    dplyr::rename(wellID = 1) 
  
 return(growth_data)
}



#Get the TECAN files, run the data function on each one, get the mumax for each then export to xls

folder_loc = "C:/Users/worms/Dropbox/PhD/Data/Tecan/2022.05.27 Test peptide marker"

file_list <- list.files(folder_loc)

for(i in 1:(file_list %>% length)){
  growth_data <- gol.get_data_tecan(file_list[i])
  
  splitted.data <- multisplit(data, c("wellID"))
  
  mumax <- data.frame(strain = character(),
                      mumax = numeric())
  
  for(j in 1:length(splitted.data)){
    
    fit <- fit_easylinear(splitted.data[[j]]$time, splitted.data[[j]]$OD %>% as.numeric())
    mumax_value <- coef(fit)[3]
    strain_value <- splitted.data[[j]]$wellID[1] %>% as.character()
    new_line <- c(strain_value, mumax_value)
    mumax <- rbind(mumax, new_line)
  }
  
  write.xlsx(mumax, paste("C:/Users/worms/Desktop/Gol_ Data for R analysis/mumax/tecan",file_list[i], sep="/"))
}

data <- get_data_tecan(folder_loc = folder_loc)


#Automate the growthcurver
folder_loc = "C:/Users/worms/Desktop/Gol_ Data for R analysis/TECAN"

file_list <- list.files(folder_loc)

for(i in 1:(file_list %>% length)){
  gc_out <- gol.get_data_tecan(file_list[i]) %>%
    select(-1) %>%
    pivot_wider(names_from = wellID, values_from = OD)%>%
    SummarizeGrowthByPlate()
  
  
  write.xlsx(gc_out, 
             paste("C:/Users/worms/Desktop/Gol_ Data for R analysis/GrowthCurver/TECAN",
                   gsub(".xlsx", "_growthcurver_output.xlsx", file_list[i]),
                   sep="/"))
}

folder_loc = "C:/Users/worms/Desktop/Gol_ Data for R analysis/Hidex"

file_list <- list.files(folder_loc)


for(i in 1:(file_list %>% length)){
  gc_out <- gol.get_data_hidex(file_list[i]) %>%
    select(2,1,3) %>%
    pivot_wider(names_from = wellID, values_from = OD)%>% 
    mutate(time = as.numeric(time)) %>%
    SummarizeGrowthByPlate()
  
  
  write.xlsx(gc_out, 
             paste("C:/Users/worms/Desktop/Gol_ Data for R analysis/GrowthCurver/Hidex",
                   gsub(".xlsx", "_growthcurver_output.xlsx", file_list[i]),
                   sep="/"))
}



folder_loc = "C:/Users/worms/Desktop/Gol_ Data for R analysis/Triplicates"

file_list <- list.files(folder_loc)

for(i in 1:(file_list %>% length)){
  growth_data <- read.csv(paste(folder_loc,file_list[i], sep = "/"))
  
  growth_data <- growth_data %>% 
    pivot_longer(!1, names_to = "strain", values_to = "OD") %>%
    dplyr::rename(time = 1) %>% 
    pivot_wider(names_from = strain, values_from = OD) %>%
    SummarizeGrowthByPlate()

  write.xlsx(gc_out, 
             paste("C:/Users/worms/Desktop/Gol_ Data for R analysis/GrowthCurver/Other",
                   gsub(".csv", "_growthcurver_output.xlsx", file_list[i]),
                   sep="/"))
  
}
