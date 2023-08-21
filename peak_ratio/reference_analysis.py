# -*- coding: utf-8 -*-
"""
Created on Fri Aug 11 15:57:27 2023
GPT code to take my reference reads and give me intensities for each position-nucleotide pair present in the references. 
@author: worms
"""
import sys
import pandas as pd

# Check if input and output file arguments were provided, otherwise set default values
if len(sys.argv) > 1:
    input_file_path = sys.argv[1]
else:
    input_file_path = "data\ref.xlsx"

if len(sys.argv) > 2:
    output_file_path = sys.argv[2]
else:
    output_file_path = "reference_averaged.xlsx"

def process_data(input_file, output_file):
    # Load the data from an Excel file
    data = pd.read_excel(input_file)

    # Extract relevant relative intensity based on the nucleotide present at that position for each peptide
    data['rel_intensity'] = data.apply(lambda row: row[f"{row['nucleotide'].lower()}_rel"], axis=1)

    # Group by position and nucleotide to compute average relative intensity, standard deviation, and count unique peptides
    grouped_data = data.groupby(['position', 'nucleotide'])
    result = pd.DataFrame({
        'average_intensity': grouped_data.rel_intensity.mean(),
        'sd_intensity': grouped_data.rel_intensity.std(),
        'unique_peptides': grouped_data.peptide.nunique()
    }).reset_index()

    # Export the results to a new Excel file
    result.to_excel(output_file, index=False)

print("Processing file "+input_file_path)
    
process_data(input_file_path, output_file_path)

print("Done.")
