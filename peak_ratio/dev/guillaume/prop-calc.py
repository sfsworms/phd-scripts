from Bio import SeqIO
import numpy as np

# Entrer ci-dessous le chemin d'accès au fichier .abi. Dans les faits, seuls DATA9-12, PLOC2 et PBAS2 sont utilisés.
# NB: Voir s'il suffit de garder le fichier .abi dans le même que celui du programme pour fonctionner
FILE_MIX = r"D:\AgroParisTech\2A\Stage\UCLouvain\Lecture de chromatogramme\L107 (NNB)6_81.ab1"
#FILE_CLONE_1 = 
#FILE_CLONE_2 = 

CHANNELS = ["DATA9", "DATA10", "DATA11", "DATA12", "PLOC1","PLOC2","PBAS1","PBAS2"]


def extract_data_from_abi(abi_file=FILE_MIX):
    """
    A partir d'un fichier ab1 issu d'un séquençage, cette fonction permet d'extraire les informations explicitées par CHANNELS
    Les choix du fichier FILE et des cannaux CHANNELS se font en dehors de la fonction.
    
    NB : Les DATA et le PLOC1-2 sont de type tuple mais PBAS1-2 est de type byte
    Suivre byte_to_string pour la conversion en chaîne de caractère et ainsi faciliter l'utilisation
    
    Entrée : fichier FILE (.abi)
    Sortie : Un dictionnaire qui regroupe tous les cannaux séléctionnés
    """
    record = SeqIO.read(abi_file,"abi")
    trace = {}
    for c in CHANNELS:
        trace[c] = record.annotations["abif_raw"][c]
    return trace


FICHIER_EXTRACTED = extract_data_from_abi(abi_file=FILE_MIX)
POSITION = 359
TIMES = 4

# bof bof start
def prop_pic(position = POSITION, data = FICHIER_EXTRACTED):
    """
    Permet de déterminer la proportion de chaque pic (DATA) dans le fichier abi extrait
    Entrée : POSITION à la quelle on veut regarder la prévalence des nucléotides
            FICHIER_EXTRACTED est le fichier abi qui a été extrait pour interprétation python
    Sortie : Un array indiquant les proportions selon les nucléotides
    
    NB: On compte à partir de 0, donc position .py = position .abi -1
    """ 

    peak_high = [["DATA9","DATA10","DATA11","DATA12"],[0,0,0,0]]
    position_in_DATA = data["PLOC2"][position]
    
    peak_high[1][0] = data["DATA9" ][position_in_DATA] # G - Noir
    peak_high[1][1] = data["DATA10"][position_in_DATA] # A - Vert
    peak_high[1][2] = data["DATA11"][position_in_DATA] # T - Rouge
    peak_high[1][3] = data["DATA12"][position_in_DATA] # C - Bleu
    
    # Somme de tous les pics pour ensuite faire la moyenne pondérée
    sum_high = peak_high[1][0] + peak_high[1][1] + peak_high[1][2] + peak_high[1][3]
    
    for k in range(4):
        peak_high[1][k] /= sum_high
        
    return (peak_high) 


def prop_pic_chain(position = POSITION, data = FICHIER_EXTRACTED, times = TIMES):
    """
    Permet de déterminer la proportion de chaque nd sur plusieurs positions à la suite (TIMES)
    """
    chain = []
    for k in range (times):
        chain += prop_pic(position, data)
        position += 1
    return chain


BAR_CODE_1 = "___" #Veuillez indiquer ici le code-barre 1 tout attaché en supprimant les ___.
BAR_CODE_2 = "___" #Veuillez indiquer ici le deuxième code barre

def prop_clone(position = POSITION, data = FICHIER_EXTRACTED, bar_code_1 = BAR_CODE_1, bar_code_2 = BAR_CODE_2):  
    """
    Détermine la proportion du clone 1 dans le mélange. Pour ce faire, on part des barres codes des deux clones.
    On va ensuite chercher les proportions des 4 nucléotides à la position donnée.
    On normalise ensuite le 1er nucléotide du clone 1 entre les deux clones.
    Puis le programme tourne sur les 3 nucléotides suivants
    
    Entrées : position : position à laquelle on commence à regarder // type : int
              data : fichier extrait à partir du .abi  // type : dict
              bar_code_1 : code barre du premier clone  // type : list
              bar_code_2 : code barre du deuxième clone  // type : list
    Sortie : une liste, donnant les proportions normalisées de chaque nucléotide pour le code barre 1
    
    ex : Position 359, data
        BAR_CODE_1 = ["A","T","G","C"]   BAR_CODE_2 = ["G","T","A","C"]
         
         Le premier calcul réalisé est la détermination des proportions de chaque pique à la position 359.
         Ensuite, on regarde la proporion normalisée du premier nucléotide du code-barre 1, ici A
         Le calcul est : A / (A + G)  Car G est le premier nd du code-barre 2
         Ce résultat est stocké dans prop_clone[0]
         
         Et ainsi de suite jusqu'à la fin du code-barre 
         
    NB: les pics ne sont pas normalisés
    """
    
    prop_clone_1 = [0 for k in range(len(bar_code_1))] # Liste dans la quelle sera stocké les proportions normalisées de chaque nd
    
    bar_code_true_1 = list(bar_code_1)
    bar_code_true_2 = list(bar_code_2)
    for k in range(len(bar_code_1)):
        Prop_pic = prop_pic(position = POSITION, data = FICHIER_EXTRACTED) # On récupère les proportions de chaque nucléotides à la position donnée
        
        G = Prop_pic[1][0]
        A = Prop_pic[1][1]
        T = Prop_pic[1][2]
        C = Prop_pic[1][3]
        
        dict_prop = {"G":G, "A": A, "T": T, "C":C}
        
        position += 1
        
        prop_clone_1[k] = dict_prop[bar_code_true_1[k]] / (dict_prop[bar_code_true_1[k]] + dict_prop[bar_code_true_2[k]])

    return (prop_clone_1)
#bof bof end

##New
# position = .... On le retrouve grâce au programme de Steve
def where_in_data(position = POSITION_SEQ, data = FILE):
    """
    Permet d'indiquer la correspondance entre la position de la base lue et la position du pic dans les DATA
    """
    return (data["PLOC2"][position])

def choix_max_pic(data = FICHIER_EXTRACTED, fourchette = 5, position = POSITION):
    max = 0
    DATA = ["DATA9","DATA10","DATA11","DATA12"]
    index = 0
    for k in range(fourchette):
        for piste in DATA:
            position_in_DATA = data["PLOC2"][position-k]
            if  data[piste][position_in_DATA] > max:
                max = data[piste][position_in_DATA]
                index = k
                couleur = piste
    return (max, index, couleur)
        
def nd_to_DATA(nd = NUCLEOTIDE):
    """
    Permet de donner la correspondance entre le nucléotide et le bon fichier DATA(9-12)
    Entrée : un nucléotide sous format string
    Sortie : Un string indiquant le DATA
    """
    if nd == "G":
        return ("DATA9")
    if nd == "A":  
        return ("DATA10")
    if nd == "T":
        return ("DATA11")
    if nd == "C":
        return ("DATA12")

 

def normalisation_pic(position = POSITION, data = FILE_TO_NORMALIZE, nd = NUCLEOTIDE):
    """ 
    Permet de normaliser le pic à la position POSITION avec un pic_ref du même fichier .abi
    """
    position_in_DATA = where_in_data(position,data)
    
    DATA = nd_to_DATA(nd)
    # TODO : dans les fait le nd sera celui du code barre
    pic_ref = choix_max_pic(data = FICHIER_EXTRACTED, fourchette = 5, position = POSITION - 20)
    # On choisit à - 20 nd car c'est là où il semble avoir le moins d'erreur
        
    return(data[DATA][position_in_DATA] / pic_ref)

def normalized_peak_height(position = POSITION, data_to_normalize = FILE_TO_NORMALIZE, data_bruit = FILE_BRUIT, nd = NUCLEOTIDE):
    """
    Normalisation en tenant compte du bruit, c'est à dire le petit signal qu'on retrouve chez l'autre clone pur
    """
    return (normalisation_pic(position, data_to_normalize, nd) - normalisation_pic(position, data_bruit, nd))



def peak_ratio(position = POSITION, data_to_normalize = FILE_TO_NORMALIZE, data_bruit = FILE_BRUIT, data_pure = FILE_PURE, nd = NUCLEOTIDE):
    NPH_final = normalized_peak_height(position, data_to_normalize, data_bruit , nd )
    NPH_expected = normalized_peak_height(position, data_pure, data_bruit, nd)
    return (NPH_final / NPH_expected)


def copy_number_proportion(position, data_to_normalize_1, data_bruit_1, data_pure_1,data_to_normalize_2, data_bruit_2, data_pure_2, nd = NUCLEOTIDE):
    ratio_1 = peak(position, data_to_normalize_1, data_bruit_1, data_pure_1, nd = NUCLEOTIDE)
    ratio_2 = peak(position, data_to_normalize_2, data_bruit_2, data_pure_2, nd = NUCLEOTIDE)
    return ( ratio_1 / ( ratio_1 + ratio_2 ))



 
    
###### Pas utile après

def byte_to_string(b_file): # Utile pour convertir les données de PBAS qui sont des bytes
    """
    Les données dans PBAS1 et 2 sont de type byte, pour pouvoir les traiter on les convertit en string avec cette fonction.
    Entrée : données de type byte
    Sortie : données de type str
    """    
    return b_file.decode() # Retourne un str



def count_nd(position,liste_reads):
    """
    Permet de compter sur plusieurs fichier .abi les nucléotides à une position donnée.
    
    """
    comtpage_ATGC = np.array[["A","T","G","C"],[0,0,0,0]]
    
    nb_read = len(liste_reads)
    
    for k in (range(nb_read)): #nombre de reads à lire
        if position <= len(liste_reads[k]): #nb de nucléotides dans le read k
            if liste_reads[k][position] == "A":
                comtpage_ATGC[1][0] += 1
            if liste_reads[k][position] == "T":
                comtpage_ATGC[1][1] += 1
            if liste_reads[k][position] == "G":
                comtpage_ATGC[1][2] += 1
            if liste_reads[k][position] == "C":
                comtpage_ATGC[1][3] += 1
    return (comptage_ATGC)




GATC={"G":0.4, "A": 0.5, "T": 0.07, "C":0.03}
key_GATC = list(GATC.keys())

prop_clone_1 = GATC[key_GATC[0]] / (GATC[key_GATC[0]]+GATC[key_GATC[1]])


G= 5
CTGA={"C":0,"T":0,"G":0,"A":0 }
A=5
T=6
C=7
G=8
BC=["A","T","C","G"]














        
        