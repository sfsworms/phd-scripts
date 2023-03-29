# CylicPeptidePipeline

My pipeline for analysis. Started in april 2022 as part of my PhD project. Updated in 2023. 

# Peptide analysis pipeline

## Peptide sequence extraction

Done with peptide_extraction.R. This saves all peptides from the sequencing files to a fasta file.in a folder called "peptide_sequence". There are versions depending on the type of intein.

Output has this format: 
>NNK3 fw
TGCGCTTCTGGT
>NNK3 fw
TGCGCGTGGAGT
>NNK3 fw
TGCAATGGGCTT
>NNK3 fw
TGCCGTATAGTT

It also outputs a csv with a summary of counts.

## Turning into counts

Done with count_peptides.py. Produce a list of sequences and counts in a "pickled" file, stored in a folder called peptide_counts, or peptide_counts_aa for the AA version.

## Unpickling the files

Done with pickle_to_csv.py, produces csv files that I can import in R. Export them to a folder called either counts_csv or aa_counts_csv.

## Count analysis

Done in R. Mostly in the peptide_analysis.R file for now. This create a file with enrichment ratios.

## Reports

Done in R Markdown files in the documents/ folder.

