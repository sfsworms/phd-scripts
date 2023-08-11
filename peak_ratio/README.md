#####################################
Peak ratio analysis
#####################################

## Intro

This is the python script written by Steve and Guillaume Tong. The aim of this script is to take sequencing of mixed sequences and give a ratio of different sequences

## Input

This script needs three types of file:
	-The ab1 file with the mixed sequences
		-Should be in a "abi" subfolder
	-ab1 files of the various "pure" sequences
		-Should be in a "standards" subfolder
		-Those are needed as the relative intensity of peaks vary depending on the local sequence context
		-Because of that use, it's important that they align in the same direction as the reads .ab1
	-fasta files of the various pure sequences
		-This should be in a multipart .fasta called "default.fasta" in the "fasta" subfolder
		-All standards should be of the same length and cover the same region
		-These are used to align the various ab1 files. Without the alignment, we run into issues as ab1 don't always start at the same position
		
## Output

This script outputs the following:
	-Zoomed in plots of all the mutation regions of the various samples chromatogram
	-Three excel files with peak intensity information
		-One PIC_SAMPLES.xlsx file containing, for each samples, the relative intensity of the four chromatograms at each position mutation as well as the base for the standards.at the mutation position as well as the reference.
		-Two PIC_STANDARD_X.xlsx files containing, for each sample, the relative peaks intensity of the two standards at the position at which they differ. 
		
## Use 

### Module dependencies
The script relies on a number of modules. Prior to running it, you should create an environment and install the following modules:
	-biopython
	-pandas
	-matplotlib
	-xlswriter

### Running the script
In the folder containing the script, you should create the subfolders needed for the various inputs and deposit the relevant .ab1 and .fasta files.

Once that is done, you can run the program by navigating to the folder in the terminal, activating your environment and typing 
python main.py

The script will create an excel file named "samples.xlsx" if it doesn't exist already, with the following columns: SAMPLE_ABI_FILES	FASTA_1	FASTA_2	STANDARD_1	STANDARD_2

The sample column will be already filled. You should fill the FASTA_1 and FASTA_2 columns with the name of the reference sequence within the .fasta as they appear in the fasta subfolder. The STANDARD_1 and STANDARD_2 columns should be filled with the name of the .ab1 files corresponding to those standards. 

### TODO:
In addition to idea commented in the code, it seems to struggle with reverse reads if the template is fw. 