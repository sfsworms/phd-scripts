import glob

from Bio import SeqIO
from collections import UserList


TYPE_FASTA = "fasta"
GENERIC_FILE_FASTA = "*.fasta"


class Template:
    def __init__(self, record, filename):
        self.filename = filename
        self.name = record.id
        self.sequence = record.seq
        self.length = len(self.sequence)


class Templates(UserList):
    def __init__(self, root_dir=None, os_separator="/"):
        super().__init__()
        self.separator = os_separator
        self.root_dir = root_dir
        self.generic_files = GENERIC_FILE_FASTA
        if root_dir is not None:
            self.directory = self.root_dir + self.separator
        else:
            self.directory = ""
        self.files = glob.glob(pathname=self.generic_files, root_dir=self.root_dir)
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
