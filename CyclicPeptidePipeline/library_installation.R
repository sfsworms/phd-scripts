#This is just to install of the libraries I need for my script.


# Specify your packages
myPackages <- c("tidyverse", #General data wrangling
                 "Biostrings", #Managing biological sequences
                 "ShortRead", #Read and process large fastq of short file
                 "edgeR", 
                 "formatR", #used to clean code
                "stringi") #used to subset code

notInstalled <- myPackages[!(myPackages %in% installed.packages()[ , "Package"])]    # Extract not installed packages

if(length(notInstalled)) install.packages(notInstalled)  
# Install not installed packages

lapply(myPackages, require, character.only = TRUE)

rm(myPackages, notInstalled)
