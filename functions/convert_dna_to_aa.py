# -*- coding: utf-8 -*-
"""
Created on Fri Aug 18 14:04:57 2023

@author: worms
"""

import csv
from Bio.Seq import Seq
import sys


# Check if input and output file arguments were provided, otherwise set default values
if len(sys.argv) > 1:
    input_file = sys.argv[1]
else:
    input_file = "dna.csv"

if len(sys.argv) > 2:
    output_file  = sys.argv[2]
else:
    output_file = "aa.csv"


def translate_csv_to_aa(input_file, output_file):
    with open(input_file, 'r') as csv_file:
        csv_reader = csv.reader(csv_file)
        
        # Assuming the DNA sequences are in the first column of the CSV file.
        dna_sequences = [row[0] for row in csv_reader]
        
        aa_sequences = [Seq(dna).translate() for dna in dna_sequences]
        
    with open(output_file, 'w', newline='') as out_file:
        csv_writer = csv.writer(out_file)
        for aa in aa_sequences:
            csv_writer.writerow([aa])


translate_csv_to_aa(input_file, output_file)