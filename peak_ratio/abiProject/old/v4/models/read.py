import glob
import numpy as np

from Bio import SeqIO
from Bio.Seq import Seq
from collections import UserList


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


class Read:
    def __init__(self, record, filename, template=None):
        self.filename = filename
        self.name = self.filename[:-LEN_EXTENSION_AB1]
        self.trace = dict()
        self.extract_data_from_abi(record)

        self.sequence = self.get_sequence()
        self.sequence_reverse = self.get_sequence_reverse()
        self.chromatograms = self.get_chromatograms()
        self.chromatograms_reverse = self.get_chromatograms_reverse()
        self.base_location = self.get_base_location()
        self.base_location_reverse = self.get_base_location_reverse()
        self.length = len(self.sequence)

        if template is not None:  # TODO?
            pass
        self.strand = None
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
        self.strand = aligner.get_strand(template=template.sequence, query=self.sequence)

    # def extract_trace(self):
    #     # Oblige de définir le strand avant
    #     if self.strand is None:
    #         return
    #     # Ne fais rien si une valeur est déjà présente
    #     if self.sequence is not None or self.chromatograms is not None or self.base_location is not None:
    #         return
    #     # Pour le sens "+"
    #     if self.strand == "+":
    #         self.sequence = self.get_sequence()
    #         self.chromatograms = self.get_chromatograms()
    #         self.base_location = self.get_base_location()
    #     # Pour le sens "-"
    #     else:
    #         self.sequence = self.get_sequence_reverse()
    #         self.chromatograms = self.get_chromatograms_reverse()
    #         self.base_location = self.get_base_location_reverse()

    # def set_length_sequence(self):
    #     # Oblige de définir la sequence avant
    #     if self.sequence is None:
    #         return
    #     # Ne fais rien si une valeur est déjà présente
    #     if self.length is not None:
    #         return
    #     self.length = len(self.sequence)

    def set_alignment(self, aligner, template):
        # Ne fais rien si une valeur est déjà présente
        if self.alignment is not None:
            return
        self.alignment = aligner.run(template=template.sequence, query=self.sequence)
        self.template_aligned = self.alignment.aligned[0]
        self.read_aligned = self.alignment.aligned[1]

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
        for location in sample_locations:
            for region in self.template_aligned:
                if location in range(region[0], region[1]):
                    tuple_index.append(self.template_aligned.index(region))
        return tuple_index

    def get_locations(self, sample_locations):
        tuple_index = self.find_tuple_index_in_aligned(sample_locations)
        index = tuple_index[0]
        # Vérifie que toutes les mutations sont contenues dans un tuple identique
        for value in tuple_index:
            if value != index:
                # "ERROR set_read_locations: Gap detected!")
                return None
        # converti l'index mutation_locations en index read_locations
        locations = []
        for location in sample_locations:
            locations.append(location - (self.template_aligned[index][0] - self.read_aligned[index][0]))
        return locations


class Reads(UserList):
    def __init__(self, root_dir=None, os_separator="/"):
        super().__init__()
        self.separator = os_separator
        self.root_dir = root_dir
        self.generic_files = GENERIC_FILE_AB1
        if root_dir is not None:
            self.directory = self.root_dir + self.separator
        else:
            self.directory = ""
        self.files = sorted(glob.glob(pathname=self.generic_files, root_dir=self.root_dir))
        self.create_reads()

    def create_reads(self):
        # Créer les objets Read() si la liste est vide
        if not self:
            for file in self.files:
                path = self.directory + file
                record = SeqIO.read(path, TYPE_AB1)
                self.append(Read(record=record, filename=file))
