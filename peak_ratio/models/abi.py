"""
Edited on Tue Aug  8 22:24:48 2023

@author: worms

This section contains code to make a list of the ab1 files in a folder and read the sequences using SeqIO.read, this leaves all the info that were in the ab1 file.
"""

import os
import glob
import numpy as np

from Bio import SeqIO
from Bio.Seq import Seq


TYPE_AB1 = "abi"
GENERIC_FILE_AB1 = "*.ab1"
EXTENSION_AB1 = ".ab1"
LEN_EXTENSION_AB1 = len(EXTENSION_AB1)

ABIF_RAW = "abif_raw"

DATA9 = "DATA9"
DATA10 = "DATA10"
DATA11 = "DATA11"
DATA12 = "DATA12"
PLOC = "PLOC1"
PBAS = "PBAS1"
ABI_CHANNELS = [DATA9, DATA10, DATA11, DATA12, PLOC, PBAS]

# This class is there to store the info existing in a .ab1 file and some related info. It takes a record produced by seqIO and a file name
class Read:
    def __init__(self, record, filename):
        self.filename = filename
        self.name = self.filename[:-LEN_EXTENSION_AB1]
        self.trace = dict()
        self.extract_data_from_abi(record)
        self.pbas = self.get_sequence() # Why is there a self.pbas and also a self.sequence?
        self.strand = None
        self.sequence = None
        self.chromatograms = None
        self.base_location = None
        self.length_sequence = None
        self.alignment = None
        self.template_aligned = None
        self.read_aligned = None

    def extract_data_from_abi(self, record):
        for channel in ABI_CHANNELS:
            self.trace[channel] = record.annotations[ABIF_RAW][channel]

    def set_strand(self, aligner, template):
        # Ne fais rien si une valeur est déjà présente
        if self.strand is not None:
            return
        self.strand = aligner.get_strand(template=template.sequence, query=self.pbas)

    def extract_trace(self):
        # Oblige de définir le strand avant
        if self.strand is None:
            return
        # Ne fais rien si une valeur est déjà présente
        if self.sequence is not None or self.chromatograms is not None or self.base_location is not None:
            return
        # Pour le sens "+"
        if self.strand == "+":
            self.sequence = self.get_sequence()
            self.chromatograms = self.get_chromatograms()
            self.base_location = self.get_base_location()
        # Pour le sens "-"
        else:
            self.sequence = self.get_sequence_reverse()
            self.chromatograms = self.get_chromatograms_reverse()
            self.base_location = self.get_base_location_reverse()

    def set_length_sequence(self):
        # Oblige de définir la sequence avant
        if self.sequence is None:
            return
        # Ne fais rien si une valeur est déjà présente
        if self.length_sequence is not None:
            return
        self.length_sequence = len(self.sequence)

    # Set the alignment for the read.
    def set_alignment(self, aligner, template):
        # Ne fais rien si une valeur est déjà présente
        if self.alignment is not None:
            return
        self.alignment = aligner.run(template=template.sequence, query=self.sequence) #Align the read on the template, the query is the sequence
        self.template_aligned = self.alignment.aligned[0]
        self.read_aligned = self.alignment.aligned[1]

    # This method extract a sequence from the ab1 trace file
    def get_sequence(self):
        return Seq(self.trace[PBAS].decode())

    def get_sequence_reverse(self):
        return self.get_sequence().reverse_complement()

    def get_chromatograms(self):
        chromatograms = {"A": list(self.trace[DATA10]),
                         "C": list(self.trace[DATA12]),
                         "G": list(self.trace[DATA9]),
                         "T": list(self.trace[DATA11])}
        return chromatograms

    def get_chromatograms_reverse(self):
        chromatograms = {"A": list(reversed(self.trace[DATA11])),
                         "C": list(reversed(self.trace[DATA9])),
                         "G": list(reversed(self.trace[DATA12])),
                         "T": list(reversed(self.trace[DATA10]))}
        return chromatograms

    # Trace[PLOC] returns the position of all the individual base call within the abi file.
    # PLOC can be "PLOC1" or "PLOC2" depending on the machine who made the sequencing
    def get_base_location(self):
        return list(self.trace[PLOC])

    def get_base_location_reverse(self):
        base_location_array = np.array(self.get_base_location())
        chromatogram_index_value = len(self.chromatograms["A"]) - 1
        base_location = chromatogram_index_value - base_location_array
        return list(reversed(base_location))

    def find_tuple_index_in_aligned(self, sample_locations):
        # Trouve l'index du tuple qui concerne les positions de sample_locations
        tuple_index = []
        # Check that sample_locations isn't empty

            
        for location in sample_locations:
            for region in self.template_aligned:
                if location in range(region[0], region[1]):
                    tuple_index.append(self.template_aligned.index(region))
        return tuple_index

    # This method takes a tuple of positions of mutations on the reference sequence, and find the corresponding position on the read
    def get_locations_on_read(self, sample_locations):
        tuple_index = self.find_tuple_index_in_aligned(sample_locations)
        if not tuple_index:
            print("Error : tuple_index is empty. Check if all reads align to their fasta templates.")
            return None
        index = tuple_index[0]
        # Vérifie que toutes les mutations sont contenues dans un tuple identique
        for value in tuple_index:
            if value != index:
                print("ERROR set_read_locations: Gap detected!")
                return None
        # converti l'index mutation_locations en index read_locations
        locations = []
        for location in sample_locations:
            locations.append(location - (self.template_aligned[index][0] - self.read_aligned[index][0])) # This converts the index on the fasta to an index on the ab1, by comparing the alignment of the template with the alignment of the read
        return locations

# The Reads class has for a goal to get all the ab1 files from a folder, read the sequences using SeqIO and make a list with the sequences and filename

class Reads(list): # Inehrit from the list class
    def __init__(self, root_dir=None):
        super().__init__() # Init with the above class
        if root_dir is not None:
            self.directory = root_dir + os.sep # If root dir exists, add a directory attribute with the directory
        else:
            self.directory = ""
        self.files = sorted(glob.glob(pathname=GENERIC_FILE_AB1, root_dir=root_dir)) # Make a sorted list of all the ab1 files in the root_dir. the constant encodes  '*.ab1"
        self.create_reads()

    def create_reads(self):
        # create read objects if the list is empty
        if not self: # Commonly used in python to check if lists are empty
            for file in self.files:
                path = self.directory + file
                record = SeqIO.read(path, TYPE_AB1) # With type ab1 this has all the ab1 info including the chromatogram
                self.append(Read(record=record, filename=file)) # Add the read reacord with the filename
