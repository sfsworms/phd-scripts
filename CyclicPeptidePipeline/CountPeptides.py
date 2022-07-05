from Bio import SeqIO
import pickle
FILE = "Cytoplasmic-NNK-Gen-1-LB_R1_001_peptide3.fasta" # Nom du fichier qu'on veut split, qui doit être dans le même répertoire
FILE_TYPE = "fasta"

# Vérification de la longueur des peptides
PEPTIDE_LENGTH = 24
NUCLEOTIDES = list("ACGT")

def check_size(size=PEPTIDE_LENGTH):
    """ Parcours l'intégralité du fichier FILE pour déterminer si les peptides ont la même longueur"""
    for record in SeqIO.parse(FILE, FILE_TYPE):
        if len(record) != size:
            return
    return True

# Comptage
def get_count(fasta_file, peptides=dict()):
    """Compte l'occurence de chaque peptide.
    Sortie: un dictionnaire dont les clés sont des séquences (objet seq)"""
    for record in SeqIO.parse(fasta_file, FILE_TYPE):
        try:
            peptides[record.seq] += 1
        except:
            peptides[record.seq] = 1
    return peptides


# Pour sauvegarder le fichier
# NB : avec pickle on sauvegarde un seul objet
SAVED_FILE = open("Count Saves bis", mode="ab")
peptides = get_count(FILE)
pickle.dump(peptides, SAVED_FILE)

Readable_file = open("Count Saves bis", mode="rb")
peptides = pickle.load(Readable_file)
# print(peptides) permet l'affichage du dictionnaire
