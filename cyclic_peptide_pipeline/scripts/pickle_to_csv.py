# -*- coding: utf-8 -*-
"""
Created on Mon Jul 11 11:52:55 2022

@author: worms

Goal is to be able to take the pickled dictionarries created by count_peptides.py and turn them into files that can be imported in
R for use with DESEq or dowstream analysis.

Input is a folder with pickled count files (ending in '_pickle')
Output is a folder containing the unpickled csv
"""

import pickle # To unpickle things
import os # To get list of files
import csv # To export to csv
import sys # To get the arguments

if len(sys.argv) == 1:
    print("Error: you need to provide at least a source folder")
    sys.exit(1)
elif len(sys.argv) == 2:
    SOURCE_FOLDER = sys.argv[1] # Nom du dossier avec les pickles
    DESTINATION_FOLDER = os.path.join(os.path.dirname(SOURCE_FOLDER), "peptide_csv") # Stocke les pickles dans un dossier "count_csv" a coté du dossier source
elif len(sys.argv) == 3:
    SOURCE_FOLDER = sys.argv[1] # Nom du dossier avec les pickles
    DESTINATION_FOLDER = sys.argv[2] # Nom du dossier de destination

file_list = []

for file in os.listdir(SOURCE_FOLDER):
    if file.endswith("_pickle"):
        file_list.append(file)
    
if len(file_list) == 0:
    print("Source folder doesn't contain files ending in '_pickle'")
    sys.exit(1)

def write_csv(file_to_write, nom):
    with open(nom, 'w') as f:
        writer = csv.writer(f)
        for k, v in file_to_write.items():
            writer.writerow([k, v])

for file in file_list:
    print("Currently treating "+file)
    readable_file = open(SOURCE_FOLDER+"\\"+file, mode="rb")
    peptides = pickle.load(readable_file)
    write_csv(peptides, DESTINATION_FOLDER+"\\"+file.replace("_pickle",".csv"))
    
    



