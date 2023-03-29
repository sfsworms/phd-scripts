# Dada2 pair-end merging

## This script will serve to merge the forward and reverse reads or the files.

library(dada2)

path <- "C:/Users/worms/NGS Data/2022.06.07_drift_seq/90-666155004b/00_fastq/NNK"
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern = "_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2_001.fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

plotQualityProfile(fnFs[1])
# In gray-scale is a heat map of the frequency of each quality score at each base position. The mean
# quality score at each position is shown by the green line, and the quartiles of the quality score
# distribution by the orange lines. The red line shows the scaled proportion of reads that extend to
# at least that position (this is more useful for other sequencing technologies, as Illumina reads are
# typically all the same length, hence the flat red line).  The forward reads are good quality. We
# generally advise trimming the last few nucleotides to avoid less well-controlled errors that can
# arise there. These quality profiles do not suggest that any additional trimming is needed. We will
# truncate the forward reads at position 240 (trimming the last 10 nucleotides).

plotQualityProfile(fnRs[1:2])
plotQualityProfile(fnRs[1:3])
#Reverse read are usually not as good and need to be trimmed more aggressively.

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names


#The below truncate the fw and reverse reads
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(150,150), #Those are the cutoff for fw and rv
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE
head(out)

# The DADA2 algorithm makes use of a parametric error model (err) and every amplicon dataset has a different set of error rates. The learnErrors method learns this error model from the data

errF <- learnErrors(filtFs, multithread=TRUE)

errR <- learnErrors(filtRs, multithread=TRUE)

# We are now ready to apply the core sample inference algorithm to the filtered and trimmed sequence data.

dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
# We now merge the forward and reverse reads together to obtain the full denoised sequences. Merging is performed by aligning the denoised forward reads with the reverse-complement of the corresponding denoised reverse reads, and then constructing the merged “contig” sequences. By default, merged sequences are only output if the forward and reverse reads overlap by at least 12 bases, and are identical to each other in the overlap region (but these conditions can be changed via function arguments).

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])
