import pandas as pd

import view

import models.plot as plot

from models.abi import Reads
from models.fasta import Templates
from models.aligner import LocalAlignment
from models.sample import Sampling
from models.sample import Samples
from models.sample import Excel


class Controller:
    def __init__(self):
        self.running = True
        self.reads: list = Reads(root_dir="abi")
        self.standards: list = Reads(root_dir="standards")
        self.templates: list = Templates(root_dir="fasta")
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

        for sample in self.samples["SAMPLES"]:
            plot.chromatogram_read(sample)


        read_list = []
        pic1_list = []
        pic2_list = []
        reference_list = []
        df = pd.DataFrame()
        for sample in self.samples["SAMPLES"]:
            # Obtention de la position de la référence
            NORMALIZATION_DISTANCE = 5
            if sample.read.strand == "+":
                reference_location = sample.sample_locations[0] - NORMALIZATION_DISTANCE
            else:
                reference_location = sample.sample_locations[-1] + NORMALIZATION_DISTANCE

                # Obtention des clés de chromatograms
            for sample_location in sample.sample_locations:
                key_1 = sample.template_1.sequence[sample_location]
                key_2 = sample.template_2.sequence[sample_location]
                key_reference = sample.template_1.sequence[reference_location]

                # Obtention des chromatograms
                chromatogram_1 = sample.read.chromatograms[key_1]
                chromatogram_2 = sample.read.chromatograms[key_2]
                chromatogram_reference = sample.read.chromatograms[key_reference]

                # Obtention de la hauteur des pics TODO faire une liste !!! pour mutations multiples
                pic_value_1 = chromatogram_1[sample.read.base_location[sample.read.get_locations([sample_location])[0]]]
                pic_value_2 = chromatogram_2[sample.read.base_location[sample.read.get_locations([sample_location])[0]]]
                pic_reference = chromatogram_reference[sample.read.base_location[sample.read.get_locations([reference_location])[0]]]

                read_list.append(sample.read.name)
                pic1_list.append(pic_value_1)
                pic2_list.append(pic_value_2)
                reference_list.append(pic_reference)

        df["SAMPLE"] = read_list
        df["PIC_1"] = pic1_list
        df["PIC_2"] = pic2_list
        df["REFERENCE +5"] = reference_list
        df.to_excel("PIC_SAMPLES.xlsx")

        # STD_1
        read_list = []
        pic1_list = []
        pic2_list = []
        reference_list = []
        df = pd.DataFrame()
        for sample in self.samples["SAMPLES"]:
            # Obtention de la position de la référence
            NORMALIZATION_DISTANCE = 5
            if sample.standard_1.strand == "+":
                reference_location = sample.sample_locations[0] - NORMALIZATION_DISTANCE
            else:
                reference_location = sample.sample_locations[-1] + NORMALIZATION_DISTANCE

                # Obtention des clés de chromatograms
            for sample_location in sample.sample_locations:
                key_1 = sample.template_1.sequence[sample_location]
                key_2 = sample.template_2.sequence[sample_location]
                key_reference = sample.template_1.sequence[reference_location]

                # Obtention des chromatograms
                chromatogram_1 = sample.standard_1.chromatograms[key_1]
                chromatogram_2 = sample.standard_1.chromatograms[key_2]
                chromatogram_reference = sample.standard_1.chromatograms[key_reference]

                # Obtention de la hauteur des pics TODO faire une liste !!! pour mutations multiples
                pic_value_1 = chromatogram_1[sample.standard_1.base_location[sample.standard_1.get_locations([sample_location])[0]]]
                pic_value_2 = chromatogram_2[sample.standard_1.base_location[sample.standard_1.get_locations([sample_location])[0]]]
                pic_reference = chromatogram_reference[
                    sample.standard_1.base_location[sample.standard_1.get_locations([reference_location])[0]]]

                read_list.append(sample.read.name)
                pic1_list.append(pic_value_1)
                pic2_list.append(pic_value_2)
                reference_list.append(pic_reference)

        df["SAMPLE"] = read_list
        df["PIC_1"] = pic1_list
        df["PIC_2"] = pic2_list
        df["REFERENCE +5"] = reference_list
        df.to_excel("PIC_STANDARD_1.xlsx")

        # STD_2
        read_list = []
        pic1_list = []
        pic2_list = []
        reference_list = []
        df = pd.DataFrame()
        for sample in self.samples["SAMPLES"]:
            # Obtention de la position de la référence
            NORMALIZATION_DISTANCE = 5
            if sample.standard_2.strand == "+":
                reference_location = sample.sample_locations[0] - NORMALIZATION_DISTANCE
            else:
                reference_location = sample.sample_locations[-1] + NORMALIZATION_DISTANCE

                # Obtention des clés de chromatograms
            for sample_location in sample.sample_locations:
                key_1 = sample.template_1.sequence[sample_location]
                key_2 = sample.template_2.sequence[sample_location]
                key_reference = sample.template_1.sequence[reference_location]

                # Obtention des chromatograms
                chromatogram_1 = sample.standard_2.chromatograms[key_1]
                chromatogram_2 = sample.standard_2.chromatograms[key_2]
                chromatogram_reference = sample.standard_2.chromatograms[key_reference]

                # Obtention de la hauteur des pics TODO faire une liste !!! pour mutations multiples
                pic_value_1 = chromatogram_1[sample.standard_2.base_location[sample.standard_2.get_locations([sample_location])[0]]]
                pic_value_2 = chromatogram_2[sample.standard_2.base_location[sample.standard_2.get_locations([sample_location])[0]]]
                pic_reference = chromatogram_reference[sample.standard_2.base_location[sample.standard_2.get_locations([reference_location])[0]]]

                read_list.append(sample.read.name)
                pic1_list.append(pic_value_1)
                pic2_list.append(pic_value_2)
                reference_list.append(pic_reference)

        df["SAMPLE"] = read_list
        df["PIC_1"] = pic1_list
        df["PIC_2"] = pic2_list
        df["REFERENCE +5"] = reference_list
        df.to_excel("PIC_STANDARD_2.xlsx")







        #            read_locations = sample.read.get_locations(sample.sample_locations)
#            print(sample.read.base_location)
#            print(sample.read.sequence)
#            print("A", sample.read.chromatograms["A"])
#            print("C", sample.read.chromatograms["C"])
#            print("G", sample.read.chromatograms["G"])
#            print("T", sample.read.chromatograms["T"])
