import glob
import os
import xlsxwriter
import pandas as pd
from Bio import SeqIO

import view
from models.read import Read
from models.fasta import Templates
from models.sample import MetadataReader
from models.aligner import LocalAlignment


ABI_DIR = "abi"
STANDARDS_DIR = "standards"
FASTA_DIR = "fasta"

GENERIC_FILE_AB1 = "*.ab1"
GENERIC_FILE_FASTA = "*.fasta"
TYPE_AB1 = "abi"
TYPE_FASTA = "fasta"

METADATA_FILENAME = "samples.xlsx"
SAMPLE_ABI_FILES = "SAMPLE_ABI_FILES"
SAMPLE_TEMPLATE_1 = "FASTA_1"
SAMPLE_TEMPLATE_2 = "FASTA_2"
SAMPLE_STANDARD_1 = "STANDARD_1"
SAMPLE_STANDARD_2 = "STANDARD_2"
SAMPLING_CHANNELS = [SAMPLE_ABI_FILES, SAMPLE_TEMPLATE_1, SAMPLE_TEMPLATE_2, SAMPLE_STANDARD_1, SAMPLE_STANDARD_2]
WIDTH_COLUMN = 30

READ_COLUMN = "READ"
TEMPLATE_1_COLUMN = "TEMPLATE_1"
TEMPLATE_2_COLUMN = "TEMPLATE_2"
STANDARD_1_COLUMN = "STANDARD_1"
STANDARD_2_COLUMN = "STANDARD_2"
SAMPLE_COLUMNS = [READ_COLUMN, TEMPLATE_1_COLUMN, TEMPLATE_2_COLUMN, STANDARD_1_COLUMN, STANDARD_2_COLUMN]


class Controller:
    def __init__(self):
        # Séparateur en fonction de l'OS
        self.separator = os.sep
        # Création des listes de fichiers de données
        self.sample_files = glob.glob(pathname=GENERIC_FILE_AB1, root_dir=ABI_DIR)
        self.standard_files = glob.glob(pathname=GENERIC_FILE_AB1, root_dir=STANDARDS_DIR)
        self.fasta_files = glob.glob(pathname=GENERIC_FILE_FASTA, root_dir=FASTA_DIR)
        # Vérification que les fichiers ont bien été trouvés
        self.verify_data_files()
        # Création d'une liste d'identifiants pour les fichiers fasta
        self.fasta_names = self.get_fasta_names_form_files()
        # Chargement des paramètres de l'expérience dans un dataframe
        self.samples = pd.DataFrame(columns=SAMPLE_COLUMNS)
        self.loading_samples_from_excel()
        # TODO

        # self.reads: list = Reads(root_dir=ABI_DIR)
        # self.standards: list = Reads(root_dir=STANDARDS_DIR)
        # self.templates: list = Templates(root_dir=FASTA_DIR)
        # self.metadata_reader = MetadataReader(reads=self.reads,
        #                                       templates=self.templates,
        #                                       standards=self.standards)
        self.aligner = LocalAlignment()

    def verify_data_files(self):
        if not self.sample_files:
            view.file_not_found(ABI_DIR)
            exit()
        if not self.standard_files:
            view.file_not_found(STANDARDS_DIR)
            exit()
        if not self.fasta_files:
            view.file_not_found(FASTA_DIR)
            exit()

    def get_fasta_names_form_files(self):
        names = []
        for file in self.fasta_files:
            path = FASTA_DIR + self.separator + file
            for record in SeqIO.parse(path, TYPE_FASTA):
                names.append(record.id)
        return names

    def loading_samples_from_excel(self):  # TODO ajouter une condition si tableau incomplet
        if not glob.glob(METADATA_FILENAME):
            self.create_setup_file()
            view.configuration_file_not_found(filename=METADATA_FILENAME)
        else:
            view.configuration_file_found(filename=METADATA_FILENAME)

        view.loading_samples()
        # Charge les objets
        metadata = pd.read_excel(io=METADATA_FILENAME)
        metadata.apply(self.set_sample, 0)
        # samples[READ_COLUMN] = metadata[SAMPLE_ABI_FILES].apply(self.set_read_from_filename, 0)
        # samples[TEMPLATE_1_COLUMN] = metadata[SAMPLE_TEMPLATE_1].apply(self.template_from_name, 0)
        # samples[TEMPLATE_2_COLUMN] = metadata[SAMPLE_TEMPLATE_2].apply(self.template_from_name, 0)
        # samples[STANDARD_1_COLUMN] = metadata[SAMPLE_STANDARD_1].apply(self.set_standard_from_filename, 0)
        # samples[STANDARD_2_COLUMN] = metadata[SAMPLE_STANDARD_2].apply(self.set_standard_from_filename, 0)
        # Charge les options
        for header in metadata.columns.tolist():
            if header not in SAMPLING_CHANNELS:
                self.samples[header] = metadata[header]

    # def set_read_from_filename(self, filename):
    #     path = ABI_DIR + self.separator + filename
    #     read = Read(record=SeqIO.read(path, TYPE_AB1), filename=filename)
    #     return read
    # TODO intégrer les FASTA mais aussi, modifier la création ou les objets reads pour une seule init !!
    def set_sample(self, metadata):
        pass
    #     return list(filter(lambda obj: obj.name == name, self.templates))[0]
    #
    # def set_standard_from_filename(self, filename):  # TODO les objets devraient être des copies !!!
    #     path = STANDARDS_DIR + self.separator + filename
    #     read = Read(record=SeqIO.read(path, TYPE_AB1), filename=filename)
    #     return read

    def create_setup_file(self):
        workbook = xlsxwriter.Workbook(METADATA_FILENAME)
        worksheet = workbook.add_worksheet()
        # Add a bold format to use to highlight cells.
        bold = workbook.add_format({'bold': True})

        for column, header in enumerate(SAMPLING_CHANNELS):
            # SAMPLING_CHANNELS = headers
            worksheet.write(0, column, header, bold)

            # SAMPLE_ABI_FILES
            if header == SAMPLE_ABI_FILES:
                for row, abi_file in enumerate(self.sample_files):
                    worksheet.write(row + 1, column, abi_file)

            # SAMPLE_TEMPLATE_1 or SAMPLE_TEMPLATE_2
            if header == SAMPLE_TEMPLATE_1 or header == SAMPLE_TEMPLATE_2:
                worksheet.data_validation(1, column, len(self.sample_files), column,
                                          {'validate': 'list', 'source': self.fasta_names})

            # SAMPLE_STANDARD_1 or SAMPLE_STANDARD_2
            if header == SAMPLE_STANDARD_1 or header == SAMPLE_STANDARD_2:
                worksheet.data_validation(1, column, len(self.sample_files), column,
                                          {'validate': 'list', 'source': self.standard_files})

        worksheet.set_column(first_col=0, last_col=len(SAMPLING_CHANNELS) - 1, width=WIDTH_COLUMN)
        workbook.close()

    # def set_strand(self, df):
    #     df[READ_COLUMN].strand = self.aligner.get_strand(template=df[TEMPLATE_1_COLUMN].sequence, query=df[READ_COLUMN].pbas)
    #     print(df[READ_COLUMN].strand)

    # def verify_reading_direction(self):
    #     if self.samples.verify_standard_directionality():
    #         view.error_standard()
    #         self.running = False
    #     if self.samples.verify_template_directionality(aligner=self.aligner):
    #         view.error_templates()
    #         self.running = False
    #     if self.samples.verify_template_length():
    #         view.error_length_templates()
    #         self.running = False

    def run(self):
        view.show_samples(self.samples)
        # samples.apply(self.set_strand, axis=1)
        # view.end()

        # self.samples.set_reads(self.aligner)
        # self.verify_reading_direction()
        # self.samples.set_locations()
        #
        # view.show_samples(df=self.samples)
