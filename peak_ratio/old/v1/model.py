import glob
import numpy as np
import pandas as pd

from Bio import SeqIO
from Bio import Align
from Bio.Seq import Seq


FASTA_FILE = "fasta/default.fasta"
DIR_ABI = "abi"
SEPARATOR_DIR = "/"
ABI_CHANNELS = ["DATA9", "DATA10", "DATA11", "DATA12", "PLOC1", "PLOC2", "PBAS1", "PBAS2"]

SAMPLES_FILE = "samples.xlsx"
SAMPLE_PARAMETERS = ["FILE", "VARIANTS", "DAYS", "CONDITION", "REPEAT"]

MODE = 'local'
MATCH_SCORE = 1.0
MISMATCH_SCORE = -2.0
GAP_SCORE = -5


def get_abi_files():
    files = glob.glob(pathname="*.ab1", root_dir=DIR_ABI)
    return sorted(files)


def create_sample_file(file=SAMPLES_FILE):
    """Créé un fichier sample.xlsx dans le répertoire courant s'il n'est pas présent."""
    if not glob.glob(file):
        df = pd.DataFrame()
        for headers in SAMPLE_PARAMETERS:
            if headers == SAMPLE_PARAMETERS[0]:
                df[headers] = get_abi_files()
            else:
                df[headers] = None
        df.to_excel(file, index=False)
        return True
    else:
        return False


def open_fasta_file(fasta_file=FASTA_FILE):
    """
    Ouvre un fichier fasta et renvoi les séquences qu'il contient dans une liste d'objets SeqRecord.

    Args:
        fasta_file:

    Returns:

    """
    records = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        records.append(record)
    return records


def get_fasta_ids(records=None):
    """
    Obtient les identifiants des objets SeqRecord et les renvoie sous la forme d'une liste de strings.

    Args:
        records:
    """
    if records is None:
        records = open_fasta_file()
    names = []
    for record in records:
        name = str(record.id)
        names.append(name)
    return names


def get_variants(wild_type, records=None):
    """
    Créé une liste d'objets Templates à partir d'une liste d'objets SeqRecord.

    Args:
        wild_type:
        records: (list) d'objet SeqRecord. Par défaut=open_fasta_file()

    Returns:
        Une liste d'objets Template().

    """
    if records is None:
        records = open_fasta_file()
    variants = []
    for record in records:
        if record.id != wild_type.id:
            variants.append(Variant(record=record, wild_type=wild_type))
    return variants


def get_barcode_index(template, records=None):
    if records is None:
        return None
    else:
        for record in records:
            if record.id != template.id:
                if record.seq in template.sequence:
                    barcode_alignment = LocalAlignment().run(template=template.sequence, query=record.seq)
                    start_index = barcode_alignment.aligned[0][0][0]
                    stop_index = barcode_alignment.aligned[0][0][1]
                    barcode_index = list(range(start_index, stop_index))
                    return barcode_index


def get_barcodes(template, records, barcode_index):
    barcodes = []
    for record in records:
        if record.id != template.id:
            barcodes.append(Barcode(record=record, template=template, barcode_index=barcode_index))
    return barcodes


def get_reads(files, template, aligner):
    reads = []
    for file in files:
        reads.append(Read(abi_file=file, template=template, aligner=aligner))
    return reads


class Template:
    """Template"""

    def __init__(self, record):
        self.id = record.id
        self.sequence = record.seq
        self.length_sequence = len(self.sequence)


class Barcode(Template):
    def __init__(self, record, template, barcode_index):
        super().__init__(record)
        self.template = template
        self.index = barcode_index


class Gene(Template):
    def __init__(self, record):
        super().__init__(record)
        self.translation = self.sequence.translate()
        self.length_translation = len(self.translation)


class WildType(Gene):
    def __init__(self, record):
        super().__init__(record)


class Variant(Gene):
    """Template"""
    def __init__(self, record, wild_type):
        super().__init__(record)
        self.wild_type = wild_type
        self.amino_acid_index = self.get_substitution_index()
        self.amino_acid_position = self.amino_acid_index + 1
        self.amino_acid_wt = self.get_amino_acid_wild_type()
        self.amino_acid_substituted = self.get_amino_acid_substituted()
        self.substitution = self.get_substitution()
        self.codon_wt = self.get_codon_from_amino_acid_position(sequence=self.wild_type.sequence)
        self.codon_substituted = self.get_codon_from_amino_acid_position(sequence=self.sequence)
        self.mutation_index = self.get_mutation_index()
        self.mutation_position = self.mutation_index[0] + 1
        self.nucleotide_wild_type = self.get_nucleotide_wild_type()
        self.nucleotide_substituted = self.get_nucleotide_substituted()

    def get_substitution_index(self):
        """
        Find the first substitution in amino acid sequence and return the position.

        Returns:
            La position de l'acide aminé substitué dans la protéine (int).
        """
        for index in range(self.length_translation):
            if self.wild_type.translation[index] != self.translation[index]:
                return index

    def get_amino_acid_wild_type(self):
        amino_acid_wt = self.wild_type.translation[self.amino_acid_index]
        return amino_acid_wt

    def get_amino_acid_substituted(self):
        amino_acid_substituted = self.translation[self.amino_acid_index]
        return amino_acid_substituted

    def get_substitution(self):
        return str(self.amino_acid_wt) + str(self.amino_acid_position) + str(self.amino_acid_substituted)

    def get_codon_from_amino_acid_position(self, sequence):
        if self.amino_acid_position is not None:
            start = 3 * self.amino_acid_index
            stop = start + 3
            codon = sequence[start:stop]
            return codon
        else:
            return None

    def get_mutation_index(self):
        """Renvoie l'index de la première mutation sur la séquence d'ADN sous la forme d'un entier compris dans une
        liste."""
        for index in range(self.length_sequence):
            if self.wild_type.sequence[index] != self.sequence[index]:
                return [index]

    def get_nucleotide_wild_type(self):
        nucleotide_wt = self.wild_type.sequence[self.mutation_index[0]]
        return nucleotide_wt

    def get_nucleotide_substituted(self):
        nucleotide_substituted = self.sequence[self.mutation_index[0]]
        return nucleotide_substituted


class LocalAlignment(Align.PairwiseAligner):
    """PairwiseAligner pour alignement local."""
    def __init__(self, mode=MODE, match_score=MATCH_SCORE, mismatch_score=MISMATCH_SCORE, gap_score=GAP_SCORE):
        super().__init__()
        self.mode = mode
        self.match_score = match_score
        self.mismatch_score = mismatch_score
        self.gap_score = gap_score

    def run(self, template, query):
        """
        Démarre un alignement local entre deux objets Bio.SeqRecord.SeqRecord à l'aide de l'objet LocalAlignment.

        Args:
            template: (objet: Seq) qui sert de template pour l'alignement.
            query: (objet: Seq) est la séquence à aligner.

        Returns:
            Alignement de deux séquences sous la forme d'un objet Bio.Align.PairwiseAlignment.

        """
        score = self.score(template, query, "+")
        score_reverse = self.score(template, query, "-")
        if score > score_reverse:
            alignment = super().align(template, query, "+")[0]
        else:
            alignment = super().align(template, query, "-")[0]
        return alignment

    def get_strand(self, template, query):
        """Renvoie le sens du read."""
        score = self.score(template, query, "+")
        score_reverse = self.score(template, query, "-")
        if score > score_reverse:
            return "+"
        else:
            return "-"


class Read:
    def __init__(self, abi_file, template, aligner, directory=DIR_ABI, separator=SEPARATOR_DIR):
        self.file = abi_file
        self.template = template
        self.aligner = aligner
        self.path = directory + separator
        self.record = SeqIO.read(self.path + self.file, "abi")
        self.trace = self.extract_data_from_abi()
        self.strand = self.aligner.get_strand(template=self.template.sequence, query=self.get_edited_sequence())
        # Décision en fonction du sens du read
        # TODO: sortir les éléments à partir du trace en fonction du sens
        if self.strand == "+":
            self.sequence = self.get_edited_sequence()
            self.chromatograms = self.get_chromatograms()
            self.base_location = self.get_edited_base_location()
        # Pour le sens "-"
        else:
            self.sequence = self.get_edited_sequence_reverse()
            self.chromatograms = self.get_chromatograms_reverse()
            self.base_location = self.get_edited_base_location_reverse()

        self.length_sequence = len(self.sequence)
        self.alignment = aligner.run(template=self.template.sequence, query=self.sequence)
        # TODO: faire l'index => attention aux gaps!! voir remarque au controller

    def extract_data_from_abi(self):
        """
        À partir d'un fichier ab1 issu d'un séquençage, cette fonction permet d'extraire les informations explicitées
        par CHANNELS. Les choix des canaux CHANNELS se font en dehors de la fonction.

        NB : Les DATA et le PLOC1-2 sont de type tuple, mais PBAS1-2 est de type byte
        Suivre byte_to_string pour la conversion en chaîne de caractère et ainsi faciliter l'utilisation

        Sortie : Un dictionnaire qui regroupe tous les canaux sélectionnés
        """
        trace = {}
        for channel in ABI_CHANNELS:
            trace[channel] = self.record.annotations["abif_raw"][channel]
        return trace

    def get_edited_sequence(self):
        edited_sequence = self.trace["PBAS1"].decode()
        return Seq(edited_sequence)

    def get_edited_sequence_reverse(self):
        edited_sequence = self.trace["PBAS1"].decode()
        return Seq(edited_sequence).reverse_complement()

    def get_edited_base_location(self):
        edited_base_location = list(self.trace["PLOC1"])
        return edited_base_location

    def get_edited_base_location_reverse(self):
        base_location_array = np.array(self.get_edited_base_location())
        chromatogram_index_value = len(self.chromatograms["A"]) - 1
        base_location = chromatogram_index_value - base_location_array
        reversed_base_location = reversed(base_location)
        return list(reversed_base_location)

    def get_chromatograms(self):
        chromatograms = {"A": list(self.trace["DATA10"]),
                         "C": list(self.trace["DATA12"]),
                         "G": list(self.trace["DATA9"]),
                         "T": list(self.trace["DATA11"])}
        return chromatograms

    def get_chromatograms_reverse(self):
        chromatograms = {"A": list(reversed(self.trace["DATA11"])),
                         "C": list(reversed(self.trace["DATA9"])),
                         "G": list(reversed(self.trace["DATA12"])),
                         "T": list(reversed(self.trace["DATA10"]))}
        return chromatograms

    def get_tuples_of_interest(self):
        pass


class Sample:  # TODO finir la  class et faire l'appel des objets dans une fonction plus haut
    def __init__(self, name, variants, day, condition, repeat):
        self.name = name
        self.variants = variants
        self.day = day
        self.condition = condition
        self.repeat = repeat
