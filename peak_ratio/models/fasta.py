import os
import glob

from Bio import SeqIO
"""
Commented on Tue Aug  8 22:24:48 2023

This file is similar to the ab1. It reads all the fasta files in 

@author: worms
"""

TYPE_FASTA = "fasta"
GENERIC_FILE_FASTA = "*.fasta"


class Template:
    def __init__(self, record, filename):
        self.filename = filename
        self.name = record.id
        self.sequence = record.seq
        self.length = len(self.sequence)


class Templates(list): #Inherit from list
    def __init__(self, root_dir=None): # root_dir is assigned None per default
        super().__init__()
        if root_dir is not None:
            self.directory = root_dir + os.sep
        else:
            self.directory = ""
        self.files = glob.glob(pathname=GENERIC_FILE_FASTA, root_dir=root_dir)
        self.create_templates()
        self.names = self.get_names()

    def create_templates(self):
        #Creat  template objects if the list is empty
        if not self:
            for file in self.files:
                path = self.directory + file
                for record in SeqIO.parse(path, TYPE_FASTA): # SeqIO parse creates an iterator so we can read all the sequences in the multipart fasta
                    self.append(Template(record=record, filename=file))

    def get_names(self):
        names = []
        for template in self:
            names.append(template.name)
        return names
