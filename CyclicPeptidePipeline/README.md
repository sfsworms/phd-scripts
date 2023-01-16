# CylicPeptidePipeline

My pipeline for analysis. Started in april 2022 as part of my PhD project. Re updated in january 23. 

# Peptide analysis pipeline

## Peptide sequence extraction

Done with peptideExtraction.R. This saves all peptides from the sequencing files to a fasta file.in a folder called "peptide_sequence"

Output has this format: 
>NNK3 fw
TGCGCTTCTGGT
>NNK3 fw
TGCGCGTGGAGT
>NNK3 fw
TGCAATGGGCTT
>NNK3 fw
TGCCGTATAGTT

## Turning into counts

Done with CountPeptides.py. Produce a list of sequences and counts in a "pickled" file
, stored in a folder called peptide_counts, or peptide_counts_aa for the AA version.

## Unpickling the files

Done with pickle_to_csv.py, produces csv files that I can import in R. Export them to a folder called either counts_csv or aa_counts_csv.

## Count analysis

Done in R. Mostly in the peptideAnalysis.R file for now.

