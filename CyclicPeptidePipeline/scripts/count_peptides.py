from Bio import SeqIO
import pickle
import sys
import os
import csv

SOURCE_FOLDER = sys.argv[1] # Nom du dossier contenant les fichiers qu'on veut split, qui doit être dans le même répertoire. Importé comme un argument du programme.
DESTINATION_FOLDER = sys.argv[2]
FILE_TYPE = "fasta"

# Vérification de la longueur des peptides
PEPTIDE_LENGTH = 12
NUCLEOTIDES = list("ACGT")

def check_size(FILE, size=PEPTIDE_LENGTH):
    """ Parcours l'intégralité du fichier FILE pour déterminer si les peptides ont la même longueur"""
    for record in SeqIO.parse(FILE, FILE_TYPE):
        if len(record) != size:
            return
    return True

# Comptage


def get_count_dna(fasta_file):
    """Compte l'occurence de chaque gène encodant un peptide.
    Sortie: un dictionnaire dont les clés sont des séquences (objet seq)"""
    peptides = dict()
    count = 0
    indiv_count = 0
    for record in SeqIO.parse(fasta_file, FILE_TYPE):
        count = count + 1
        if count % 1000000 == 0:
            print("We have processed " + str(count/1000000) + " millions records.") #Just so I know it is running
        try:
            peptides[record.seq] += 1
        except:
            peptides[record.seq] = 1
            indiv_count = indiv_count + 1
    print("We have processed " + str(count) + " records, containing " + str(indiv_count) + " unique sequences.")
    return peptides


# Pour sauvegarder le fichier
# NB : avec pickle on sauvegarde un seul objet
# This leads to a highly efficient binary file but one that can't easily be processed. 

# Get the list of .fa files to process


for file in os.listdir(SOURCE_FOLDER):
    if file.endswith(".fa"):
        file_list.append(file)

# Treat the files and save the pickled results
results = [] #Used to store the results

csv_filename = os.path.basename(SOURCE_FOLDER)  + "_summary.csv" # Create the tracking csv


for file in file_list:
    output_file_name = file.replace(".fa", "_count_pickle")
    saved_file = open(DESTINATION_FOLDER+"\\"+output_file_name, mode="ab")
    peptides_list = get_count_dna(SOURCE_FOLDER+"\\"+file)
    pickle.dump(peptides_list, saved_file) # Dump the files
    total_seqs = len(peptides_list) # Count the sequences
    total_counts = sum(peptides_list.values())
    
    with open(DESTINATION_FOLDER + "\\" + csv_filename, mode="a", newline='') as output_file:
        if output_file.tell() == 0:  # check if file is empty
            writer = csv.writer(output_file)
            writer.writerow(["file", "total_counts", "total_seqs"])
        writer = csv.writer(output_file)
        writer.writerow([file, total_counts, total_seqs])



# Readable_file = open("Count Saves bis", mode="rb")
# peptides = pickle.load(Readable_file)
# # print(peptides) permet l'affichage du dictionnaire
