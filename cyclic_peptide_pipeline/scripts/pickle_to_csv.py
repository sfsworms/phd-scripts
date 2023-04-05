# -*- coding: utf-8 -*-
"""
Created on Mon Jul 11 11:52:55 2022

@author: worms

Goal is to be able to take the pickled dictionarries created by count_peptides.py and turn them into files that can be imported in
R for use with DESEq or dowstream analysis
"""

import pickle # To unpickle things
import os # To get list of files
import csv # To export to csv
import sys # To get the arguments

SOURCE_FOLDER = sys.argv[1] # Nom du dossier avec les pickles
DESTINATION_FOLDER = sys.argv[2]

file_list = os.listdir(SOURCE_FOLDER)

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
    
    



