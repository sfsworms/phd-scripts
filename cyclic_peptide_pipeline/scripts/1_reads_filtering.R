# Dada2 pair-end merging

## This script will serve to filter the reads prior to merging with FLASH

library(tidyverse)
library(dada2)

path <- "C:/Users/Sebastian Worms/ngs_data/2023 OXA NGS/90-957741147"
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq.gz and SAMPLENAME_R2_001.fastq.gz
# The sort() function is there to make sure they're in the same order

fnFs <- sort(list.files(path, pattern = "_R1_001.fastq.gz$", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2_001.fastq.gz$", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)


# In gray-scale is a heat map of the frequency of each quality score at each base position. The mean
# quality score at each position is shown by the green line, and the quartiles of the quality score
# distribution by the orange lines. The red line shows the scaled proportion of reads that extend to
# at least that position (this is more useful for other sequencing technologies, as Illumina reads are
# typically all the same length, hence the flat red line). Here the quality is mostly good. Also
# quality is chunky here, with values at 13, 25 and 37.


png(file = file.path(path,"quality_plot.png"), bg = "transparent",  width = 600)
plotQualityProfile(c(fnFs,fnRs))
dev.off()

#Quality is pretty good here. I commented out the file because it takes time.


#Reverse read are usually not as good and need to be trimmed more aggressively. In my case, they're good.

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names

names(filtRs) <- sample.names

#The below truncate the fw and reverse reads. If assigned, "out" is a list of the files and reads obtained.

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(100,100), #Those are the cutoff for fw and rv
                     maxN=0, maxEE=c(3,5), truncQ=c(10,10), rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE

write.csv(x = out, file = file.path(path,"filter_report.csv"))


