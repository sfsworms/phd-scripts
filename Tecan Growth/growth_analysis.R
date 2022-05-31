## This link to the function to extract the data from the tecan

source(file = "C:/Users/worms/Dropbox/PhD/PhD-Scripts/Tecan Growth/utils.R")

folder.loc <- "C:/Users/worms/Dropbox/PhD/Data/Tecan/2022.05.27 Test peptide marker"

## This extract the data from the tecan file, add well names to the well ID and add details contained in 'experiment.xlsx'

data_frame <- get_data_tecan(folder_loc = folder.loc) %>%
  add_well_name(., folder_loc = folder.loc) %>%
  add_experiment_details(., folder_loc = folder.loc)

write.csv(data_frame,paste(folder.loc,"tidy.csv", sep="/"))

# Now I can graph things.

ggplot(data_frame %>% 
         filter(strain != "NA") %>% 
         filter(wellName != "bSW069") %>%
         mutate(time = time/60) %>%
         filter(time < 600), aes(x = time, y = OD, color = wellName)) +
  labs(title = "Growth of E.coli with the various barcoding peptides") +
  xlab("Time (min)") +
  ylab("OD600") +
  geom_smooth()

ggsave(paste(folder.loc,"growth.png", sep="/"))  



## This part of the script is meant to compute the various growth parameter

splitted.data <- multisplit(data_frame, c("wellID")) #Split the data in a list

growth.param <- data.frame(wellID = character(),
                    wellName = character(),
                    mumax = numeric(),
                    lag = numeric())
## TODO: Replace this with apply() function 
for(j in 1:length(splitted.data)){
    fit <- fit_easylinear(splitted.data[[j]]$time, splitted.data[[j]]$OD)
    growth.param[j,"mumax"] <- coef(fit)[3]
    growth.param[j,"lag"] <- coef(fit)[4]
    growth.param[j,"wellID"] <- splitted.data[[j]]$wellID[1]
    growth.param[j,"wellName"] <- splitted.data[[j]]$wellName[1]
}
 
growth.param <- growth.param %>%
  mutate(wellID = as.factor(wellID)) %>%
  mutate(wellName = as.factor(wellName))

ggplot(growth.param %>%
                    mutate(wellName = as.character(wellName)) %>%
                    filter(wellName != "blank") %>%
  #                  filter(wellName != "bSW069") %>%
                    filter(wellName != "LB") %>%
                    mutate(wellName = as.factor(wellName)) , aes(x = wellName, y = mumax)) +
  geom_boxplot()

ggsave(paste(folder.loc,"mumax.png", sep="/"))

ggplot(growth.param %>%
                    mutate(wellName = as.character(wellName)) %>%
                    filter(wellName != "blank") %>%
                    filter(wellName != "bSW069") %>%
                    filter(wellName != "LB") %>%
                    mutate(wellName = as.factor(wellName)) , aes(x = wellName, y = lag)) +
  geom_boxplot()

ggsave(paste(folder.loc,"lag.png", sep="/"))



write.xlsx(mumax, paste("C:/Users/worms/Desktop/Gol_ Data for R analysis/mumax/tecan",file_list[i], sep="/"))
}

