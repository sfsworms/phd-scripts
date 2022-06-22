#This is just to install of the libraries I need for my script.



my_packages <- c("tidyverse", "Biostrings", "ShortRead", "edgeR", "formatR")                                        # Specify your packages
not_installed <- my_packages[!(my_packages %in% installed.packages()[ , "Package"])]    # Extract not installed packages
if(length(not_installed)) install.packages(not_installed)                               # Install not installed packages
rm(my_packages, not_installed)

library(tidyverse)
library(Biostrings)
library(ShortRead)
library(formatR) #Automatically clean code
