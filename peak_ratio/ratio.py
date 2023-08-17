# -*- coding: utf-8 -*-
"""
Created on Mon Aug 14 11:54:10 2023

This is the follow on of the script to get %age of each peptide.

@author: worms
"""

import sys
import pandas as pd

# Check if input and output file arguments were provided, otherwise set default values
if len(sys.argv) > 1:
    reference_file_path = sys.argv[1]
else:
    reference_file_path = "data/reference_values.xlsx"

if len(sys.argv) > 2:
    sample_file_path = sys.argv[2]
else:
    sample_file_path = "data\sample.xlsx"
    
if len(sys.argv) > 3:
    output_file_path = sys.argv[2]
else:
    output_file_path = "data\percentage.xlsx"

# Load the reference values from the Excel file
reference_df = pd.read_excel(reference_file_path)

# Convert the reference dataframe to a dictionary for faster lookup.
# The key is a tuple of (position, nucleotide) and the value is the average intensity.
reference_dict = {(row['position'], row['nucleotide']): row['average_intensity'] for _, row in reference_df.iterrows()}

# Load the sample values from the Excel file
sample_df = pd.read_excel(sample_file_path)

# Function to compute the percentage for a given base and position
def compute_percentage(position, nucleotide, sample_intensity):
    # Get the reference intensity for the given position and nucleotide
    ref_intensity = reference_dict.get((position, nucleotide))
    
    # If a reference intensity exists, compute the percentage, otherwise return 0
    if ref_intensity:
        return sample_intensity / ref_intensity * 100
    else:
        return 0

# Compute the percentages for each base in the sample dataframe using the 'apply' method
sample_df['a_perc'] = sample_df.apply(lambda row: compute_percentage(row['pos'], 'A', row['a_rel']), axis=1)
sample_df['c_perc'] = sample_df.apply(lambda row: compute_percentage(row['pos'], 'C', row['c_rel']), axis=1)
sample_df['g_perc'] = sample_df.apply(lambda row: compute_percentage(row['pos'], 'G', row['g_rel']), axis=1)
sample_df['t_perc'] = sample_df.apply(lambda row: compute_percentage(row['pos'], 'T', row['t_rel']), axis=1)

# Extract only the relevant columns for final output
output_df = sample_df[['peptide', 'pos', 'a_perc', 'c_perc', 'g_perc', 't_perc']]

# Save the output dataframe to an Excel file
output_df.to_excel(output_file_path, index=False)