## Analysis project

setwd(dir = "C:/Users/worms/Dropbox/PhD/PhD-Scripts/CylicPeptidePipeline")
source("library.R")  #Install and/or load needed libraries
source("function.R")  #Functions used in the script

## This is still me playing with the count paper.

# count <- read.csv(choose.files()) From the count file of the count
# paper I can go and do the rest of hte graphes.

# Load a fastq file with the peptides. The files are too big to be used
# entirely, so I'll have to use FastqStremer

# I first need to get a connection established.

con <- file.path("C:/Users/worms/Documents/2022.06.07 Drift Seq/90-666155004b/test",
  "Cytoplasmic-NNK-Gen-1-LB_R1_001.fastq")

# Get the quality report

qaSummary <- qa(con, type="fastq")

qaSummary %>% report() %>%
  browseURL()


# I can then use the connection to read bits of the fastQ file

seq = FastqStreamer(con, 10000)

dnaSeq = yield(seq)

# Filter the one that are way too small (<145 bp out of 150, about 0.5% of seq for NNK3)

dnaSeq = dnaSeq[width(dnaSeq) > 145]

# Get all the sequence that match the fwd pattern for intein
fwdPat = "TGGCTTCATTGCGAGCAAT"

fwd = grepl(pattern = fwdPat, dnaSeq %>%
              sread() %>%
              as.character())

fwdDnaSeq <- dnaSeq[fwd]

# Get all the sequence that match the rev pattern for intein

revPat = "GTCGTAAGACAGACA"

rev = grepl(pattern = revPat, dnaSeq %>%
              sread() %>%
              as.character())

revDnaSeq <- dnaSeq[rev]



# Get all seq containing that

rev = grepl(pattern = revpat, DNAString %>%
  as.character())
revString <- DNAString[rev]

# get all seq not containing either pat or revpat
restString <- DNAString[(!(fwd | rev))]

