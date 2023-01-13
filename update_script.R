#Generic update script

if (!requireNamespace("installr", quietly = TRUE))
  install.packages("installr")
library(installr)
updateR()

update.packages()

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.16")

Sys.setenv(LANG = "en")
