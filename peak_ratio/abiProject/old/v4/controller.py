import os

import view
import models.plot as plot

from models.read import Reads
from models.fasta import Templates
from models.aligner import LocalAlignment
from models.sample import Sampling
from models.sample import Samples
from models.sample import Excel


ABI_DIR = "abi"
STANDARDS_DIR = "standards"
FASTA_DIR = "fasta"


class Controller:
    def __init__(self):
        self.running = True  # TODO voir si encore necessaire
        self.quit = False
        # Séparateur en fonction de l'OS
        self.separator = os.sep
        # Création des listes de fichiers de données et des objets TODO read
        self.reads: Reads = Reads(root_dir=ABI_DIR, os_separator=self.separator)
        self.standards: Reads = Reads(root_dir=STANDARDS_DIR, os_separator=self.separator)
        self.templates: Templates = Templates(root_dir=FASTA_DIR, os_separator=self.separator)
        # Vérification que les fichiers ont bien été trouvés
        self.verify_data_files()
        # Chargement des paramètres de l'expérience dans un dataframe TODO
        self.create_excel_file()
        self.sampling = Sampling()
        self.aligner = LocalAlignment()
        self.samples = Samples(sampling=self.sampling,
                               reads=self.reads,
                               standards=self.standards,
                               templates=self.templates)

        self.samples.set_reads(self.aligner)
        self.verify_reading_direction()
        self.samples.set_locations()

    def create_excel_file(self):
        # Ne fais rien si le fichier existe déjà
        excel = Excel(reads=self.reads, templates=self.templates, standards=self.standards)
        if excel.excel_file_exist():
            return
        excel.create()
        view.fill_sample_file(file=excel.filename)

    def verify_data_files(self):
        for data_files in [self.reads, self.standards, self.templates]:
            if not data_files.files:
                view.file_not_found(generic_file=data_files.generic_files, directory=data_files.root_dir)
                self.quit = True
        if self.quit:
            exit()

    def verify_reading_direction(self):
        if self.samples.verify_standard_directionality():
            view.error_standard()
            self.running = False
        if self.samples.verify_template_directionality(aligner=self.aligner):
            view.error_templates()
            self.running = False
        if self.samples.verify_template_length():
            view.error_length_templates()
            self.running = False

    def run(self):
        view.show_samples(df=self.samples)
