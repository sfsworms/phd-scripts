import os
import glob

from Bio import SeqIO


TYPE_FASTA = "fasta"
GENERIC_FILE_FASTA = "*.fasta"


class Template:
    def __init__(self, record, filename):
        self.filename = filename
        self.name = record.id
        self.sequence = record.seq
        self.length = len(self.sequence)


class Templates(list):
    def __init__(self, root_dir=None):
        super().__init__()
        if root_dir is not None:
            self.directory = root_dir + os.sep
        else:
            self.directory = ""
        self.files = glob.glob(pathname=GENERIC_FILE_FASTA, root_dir=root_dir)
        self.create_templates()
        self.names = self.get_names()

    def create_templates(self):
        # Créer les objets Template() si la liste est vide
        if not self:
            for file in self.files:
                path = self.directory + file
                for record in SeqIO.parse(path, TYPE_FASTA):
                    self.append(Template(record=record, filename=file))

    def get_names(self):
        names = []
        for template in self:
            names.append(template.name)
        return names
