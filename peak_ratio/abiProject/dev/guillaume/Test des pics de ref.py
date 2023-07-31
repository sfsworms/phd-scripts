"""
Le papier de Inferring relative proportions of DNA variants from sequencing electropherograms par I. M. Carr∗, montre une façon de normaliser les pics.
En effet, ces derniers peuvent avoir un intensité variable selon leur environnement. L'aauteur propose de choisir la référence entre 5 et 10 nd avant notre pic d'intérêt.
Nous allons tester ici l'influence du pic, à savoir jusqu'à quel point nous pouvons aller loin dans le choix du pic sans pour autant affecter la corrélation
"""

from Bio import SeqIO
import numpy as np
import statistics as st

CHANNELS = ["DATA9", "DATA10", "DATA11", "DATA12", "PLOC1","PLOC2","PBAS1","PBAS2"]
FILES_TEST = [r"D:\AgroParisTech\2A\Stage\UCLouvain\2021.01.15 Seq of vanco peptide\3.ab1",
              r"D:\AgroParisTech\2A\Stage\UCLouvain\2021.01.15 Seq of vanco peptide\4.ab1",
              r"D:\AgroParisTech\2A\Stage\UCLouvain\2021.01.15 Seq of vanco peptide\5.ab1",
              r"D:\AgroParisTech\2A\Stage\UCLouvain\2021.01.15 Seq of vanco peptide\6.ab1",
              r"D:\AgroParisTech\2A\Stage\UCLouvain\2021.01.15 Seq of vanco peptide\7.ab1",
              r"D:\AgroParisTech\2A\Stage\UCLouvain\2021.01.15 Seq of vanco peptide\8.ab1",
              r"D:\AgroParisTech\2A\Stage\UCLouvain\2021.01.15 Seq of vanco peptide\9.ab1"]



#FILES_TEST est une liste qui contient les chemin d'accès des différents fichiers .ab1 à étudier

def extract_data_from_abi(abi_file):
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

FILES_EXTRACTED = [extract_data_from_abi(abi_file = FILE) for FILE in FILES_TEST]


def prop_pic(position, data):
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


def prop_pic_chain(position, data, times):
    """
    Permet de déterminer la proportion de chaque nd sur plusieurs positions à la suite (TIMES)
    """
    chain = []
    for k in range (times):
        chain += prop_pic(position, data)
        position += 1
    return chain


"""
Pour les essais on utilise le pic A qui est indexé comme suit
[3,4,5,6,7,8,9] Fichier
[185,182,183,182,184,183,184] Position du pic testé dans les différents fichiers
"""
test = []
#list_position_peak = [85,82,83,82,84,83,84]
list_position_peak = [185,182,183,182,184,183,184]
#list_position_peak = [285,282,283,282,284,283,284]
#list_position_peak = [385,382,383,382,384,383,384]
for k in range(len(FILES_EXTRACTED)):
    test = test +[k] +prop_pic_chain(position = list_position_peak[k], data = FILES_EXTRACTED[k], times = 3)


def choix_max_pic(data, debut_fourchette, fin_fourchette, position):
    max = 0
    DATA = ["DATA9","DATA10","DATA11","DATA12"]
    index = 0
    for k in range(debut_fourchette, fin_fourchette):
        for piste in DATA:
            position_in_DATA = data["PLOC2"][position - k]
            if  data[piste][position_in_DATA] > max:
                max = data[piste][position_in_DATA]
                index = k
                couleur = piste
    return ([max, index, couleur]) #Dans les faits, seul max nous intéresse. index et couleur servent à mieux visualiser sur le chromatogramme




def test_pas(list_pas):
    """
    Entrée : une liste de pas, à savoir le décalage de la fenêtre dans laquelle je cherche un pic de référence
    Sortie : Deux listes
    La première liste contient les listes des pas testés avec le rapport ecart-type / moyenne pour un même pic sur plusieurs fichiers
        Elle est de taille le nb de pas 
    La deuxième liste contient les listes des valeurs des pics normalisées par pas
        Elle est de taille le nb de pas, chq élément contenant le nb de fichiers 
        
    [[pas 1, rapport 1 sd/moy], [pas 2, rapport 2 sd/moy], ...],
    [[liste du même pic sur différents fichier normalisés au pas 1],[liste du même pic sur différents fichier normalisés au pas 2], ...]
    """
    compil = []
    compil_list=[]
    for pas in list_pas:
        list_norm_1 = []
        for k in range(len(FILES_EXTRACTED)):
            choix_pic = choix_max_pic(FILES_EXTRACTED[k], pas, 5 + pas, list_position_peak[k])
            # on choisit la taille de la fenêtre comme étant de 5
            pic_norm = FILES_EXTRACTED[k]["DATA9"][list_position_peak[k]] / choix_pic[0]
            list_norm_1 += [pic_norm]
        
        compil += [[pas,st.pstdev(list_norm_1) / st.mean(list_norm_1)]]
        compil_list += [list_norm_1] #récupère la liste des valeurs normalisées
    return (compil)



# Quelques suggestions de test
compil_large = test_pas([5,10,15,20,25,30,35])
compil_large_bis = test_pas([40,45,50,55,65,70])

compil_precise = test_pas([15,16,17,18,19,20])
compil_precise_bis = test_pas([20,21,22,23,24,25])

####################### Résultats

Pour list_position_peak = [85,82,83,82,84,83,84]


compil_large
Out[14]: 
[[5, 0.05678659607669972],
 [10, 0.09580419675070607],
 [15, 0.08770704964219286],
 [20, 0.12350391814693085],
 [25, 0.11639246457112076],
 [30, 0.1120978793415546],
 [35, 0.11200329880544684]]

compil_large_bis
Out[15]: 
[[40, 0.11120897875984616],
 [45, 0.1272794610335809],
 [50, 0.132697607054403],
 [55, 0.0871678758229194],
 [65, 0.0763933943703051],
 [70, 0.06691078469685961]]


Pour list_position_peak = [185,182,183,182,184,183,184]
compil_large
Out[3]: 
[[5, 0.4844799333808041],
 [10, 0.23797733968397952],
 [15, 0.2579435159625268],
 [20, 0.20704373532562093],
 [25, 0.3513754609066071],
 [30, 0.239009761823049],
 [35, 0.22159843000464174]]

compil_large_bis
Out[4]: 
[[40, 0.30599291972156156],
 [45, 0.2763006219802879],
 [50, 0.12721225903686695],
 [55, 0.28265763939963323],
 [65, 0.25395045668746563],
 [70, 0.2588524496178935]]

Pour list_position_peak = [285,282,283,282,284,283,284]
compil_large
Out[7]: 
[[5, 0.2913948597181122],
 [10, 0.28373085516933944],
 [15, 0.24259249024516757],
 [20, 0.2502932401663851],
 [25, 0.26089307739984496],
 [30, 0.24656765862600385],
 [35, 0.2810533799028373]]

compil_large_bis
Out[8]: 
[[40, 0.2767632214146113],
 [45, 0.2590316917830092],
 [50, 0.25030809425826417],
 [55, 0.27636148619277184],
 [65, 0.271171018736835],
 [70, 0.2653475425118425]]


Pour list_position_peak = [385,382,383,382,384,383,384]
compil_large
Out[10]: 
[[5, 0.2554278023143299],
 [10, 0.2734593930366401],
 [15, 0.35274357809355544],
 [20, 0.29385725036573845],
 [25, 0.44647389986449393],
 [30, 0.4405660899838005],
 [35, 0.3584497172738224]]

compil_large_bis
Out[11]: 
[[40, 0.26872516032041616],
 [45, 0.32905578176096806],
 [50, 0.3421892848897187],
 [55, 0.3385917165939083],
 [65, 0.41794197810437533],
 [70, 0.21907015530600382]]


Test sur pic non-normalisés:
    
def erreur_pic():
    compil = []
    list_1 = []
    for k in range(len(FILES_EXTRACTED)):
        pic = FILES_EXTRACTED[k]["DATA9"][list_position_peak[k]]
        list_1 += [pic]
    
    compil += [[st.pstdev(list_1) / st.mean(list_1)]]
    return (compil)

#Inutile plus bas
    # def test_pas(list_pas):
    #     compil = []
    #     for pas in list_pas:
    #         list_ref = [] 
    #         for k in range(len(FILES_EXTRACTED)):
    #             list_ref += [choix_max_pic(data = FILES_EXTRACTED[k], debut_fourchette = pas, fin_fourchette = 5 + pas, position = list_position_peak[k])]
            
    #         list_normal = []
    #         for k in range (len(list_ref)):
    #             list_normal += [FILES_EXTRACTED[k]["DATA9"][list_position_peak[k]] / list_ref[k][1]]
                
    #         compil += [list_normal]
        
    #     # list_same_ref = [[] for k in range (len(list_pas))]
        
    #     # for i in range (len(FILES_EXTRACTED)):
    #     #     for k in range(len(list_pas)):
    #     #         list_same_ref += [compil[i][k]]
            
    #     return (compil)

