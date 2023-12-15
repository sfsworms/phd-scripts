# Dada2 pair-end merging

## This script will serve to merge the forward and reverse reads or the files. Should be run after filtering.

library(tidyverse)
library(dada2)


# File parsing
filtpath <- "path/to/FWD/filtered" # CHANGE ME to the directory containing your filtered fastq files
filts <- list.files(filtpath, pattern="fastq.gz", full.names=TRUE) # CHANGE if different file extensions
sample.names <- sapply(strsplit(basename(filts), "_"), `[`, 1) # Assumes filename = sample_XXX.fastq.gz
names(filts) <- sample.names
# Learn error rates
set.seed(100)
err <- learnErrors(filts, nbases = 1e8, multithread=TRUE, randomize=TRUE)
# Infer sequence variants
dds <- vector("list", length(sample.names))
names(dds) <- sample.names
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derep <- derepFastq(filts[[sam]])
  dds[[sam]] <- dada(derep, err=err, multithread=TRUE)
}
# Construct sequence table and write to disk
seqtab <- makeSequenceTable(dds)
saveRDS(seqtab, "path/to/run1/output/seqtab.rds") # CHANGE ME to where you want sequence table saved