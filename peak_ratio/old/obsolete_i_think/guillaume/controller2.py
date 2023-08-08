import pandas as pd

import view

# import models.plot as plot

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
            print(sample.read.filename)
            #            print(sample.read.sequence)
            print(sample.sample_locations)
            print(sample.read.get_locations(sample.sample_locations))
        #            positions = sample.read.get_locations(sample.sample_locations)

        #        for sample in self.samples["SAMPLES"]:
        #            plot.chromatogram_read(sample)

        #  Sample mix
        NORMALIZATION_DISTANCE = 20
        read_list = []
        mix_pic1_list = []
        mix_pic2_list = []
        mix_reference_list = []
        nucleotide_positions = []
        df = pd.DataFrame()
        for sample in self.samples["SAMPLES"]:
            # Obtention de la position de la référence
            #            NORMALIZATION_DISTANCE = 5
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
                pic_value_1 = chromatogram_1[
                    sample.read.base_location[sample.read.get_locations([sample_location])[0]]]
                pic_value_2 = chromatogram_2[
                    sample.read.base_location[sample.read.get_locations([sample_location])[0]]]
                pic_reference = chromatogram_reference[
                    sample.read.base_location[sample.read.get_locations([reference_location])[0]]]
                read_list.append(sample.read.name)
                mix_pic1_list.append(pic_value_1)
                mix_pic2_list.append(pic_value_2)
                mix_reference_list.append(pic_reference)
                nucleotide_positions.append(sample_location)

        df["SAMPLE"] = read_list
        df["POSITION ND"] = nucleotide_positions
        df["PIC_1"] = mix_pic1_list
        df["PIC_2"] = mix_pic2_list
        df["REFERENCE +" + str(NORMALIZATION_DISTANCE)] = mix_reference_list
        df.to_excel("PIC_SAMPLES.xlsx")

        # STD_1
        read_list = []
        std1_pic1_list = []
        std1_pic2_list = []
        std1_reference_list = []
        df = pd.DataFrame()
        for sample in self.samples["SAMPLES"]:
            # Obtention de la position de la référence
            #            NORMALIZATION_DISTANCE = 5
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
                pic_value_1 = chromatogram_1[
                    sample.standard_1.base_location[sample.standard_1.get_locations([sample_location])[0]]]
                pic_value_2 = chromatogram_2[
                    sample.standard_1.base_location[sample.standard_1.get_locations([sample_location])[0]]]
                pic_reference = chromatogram_reference[
                    sample.standard_1.base_location[sample.standard_1.get_locations([reference_location])[0]]]
                read_list.append(sample.read.name)
                std1_pic1_list.append(pic_value_1)
                std1_pic2_list.append(pic_value_2)
                std1_reference_list.append(pic_reference)

        df["SAMPLE"] = read_list
        df["PIC_1"] = std1_pic1_list
        df["PIC_2"] = std1_pic2_list
        df["REFERENCE +" + str(NORMALIZATION_DISTANCE)] = std1_reference_list
        df.to_excel("PIC_STANDARD_1.xlsx")

        # STD_2
        read_list = []
        std2_pic1_list = []
        std2_pic2_list = []
        std2_reference_list = []
        df = pd.DataFrame()
        for sample in self.samples["SAMPLES"]:
            # Obtention de la position de la référence
            #            NORMALIZATION_DISTANCE = 5
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
                pic_value_1 = chromatogram_1[
                    sample.standard_2.base_location[sample.standard_2.get_locations([sample_location])[0]]]
                pic_value_2 = chromatogram_2[
                    sample.standard_2.base_location[sample.standard_2.get_locations([sample_location])[0]]]
                pic_reference = chromatogram_reference[
                    sample.standard_2.base_location[sample.standard_2.get_locations([reference_location])[0]]]
                read_list.append(sample.read.name)
                std2_pic1_list.append(pic_value_1)
                std2_pic2_list.append(pic_value_2)
                std2_reference_list.append(pic_reference)

        df["SAMPLE"] = read_list
        df["PIC_1"] = std2_pic1_list
        df["PIC_2"] = std2_pic2_list
        df["REFERENCE +" + str(NORMALIZATION_DISTANCE)] = std2_reference_list
        df.to_excel("PIC_STANDARD_2.xlsx")

        # From here, newly added

        # Normalisation clone 1
        NPH_list_1_final = []
        for index in range(len(read_list)):
            NPH_value_1_final = (mix_pic1_list[index] / mix_reference_list[index]) - (
                    std2_pic1_list[index] / std2_reference_list[index])
            NPH_list_1_final.append(is_in_range(NPH_value_1_final))

        NPH_list_1_expected = []
        for index in range(len(read_list)):
            NPH_value_1_expected = (std1_pic1_list[index] / std1_reference_list[index]) - (
                    std2_pic1_list[index] / std2_reference_list[index])
            NPH_list_1_expected.append(is_in_range(NPH_value_1_expected))

        ratio_list_1 = []
        for index in range(len(read_list)):
            ratio_value = NPH_list_1_final[index] / NPH_list_1_expected[index]
            ratio_list_1.append(ratio_value)

        # Normalisation clone 2
        NPH_list_2_final = []
        for index in range(len(read_list)):  # mix_read_list est le nb d'échantillon
            NPH_value_2_final = (mix_pic2_list[index] / mix_reference_list[index]) - (
                    std1_pic2_list[index] / std1_reference_list[index])
            NPH_list_2_final.append(is_in_range(NPH_value_2_final))

        NPH_list_2_expected = []
        for index in range(len(read_list)):
            NPH_value_2_expected = (std2_pic2_list[index] / std2_reference_list[index]) - (
                    std1_pic2_list[index] / std1_reference_list[index])
            NPH_list_2_expected.append(is_in_range(NPH_value_2_expected))

        ratio_list_2 = []
        for index in range(len(read_list)):
            ratio_value = NPH_list_2_final[index] / NPH_list_2_expected[index]
            ratio_list_2.append(ratio_value)

        # Calcul du Copy Number Proportion
        CNP_list_1 = []
        for index in range(len(read_list)):
            CNP_value_1 = ratio_list_1[index] / (ratio_list_1[index] + ratio_list_2[index])
            CNP_list_1.append(CNP_value_1)

        # Création fichier excel
        df = pd.DataFrame()

        df["SAMPLE"] = read_list
        df["POSITION ND"] = nucleotide_positions
        df["PIC_1"] = mix_pic1_list
        df["PIC_2"] = mix_pic2_list
        df["REFERENCE +" + str(NORMALIZATION_DISTANCE)] = mix_reference_list

        df["NPH_1_final"] = NPH_list_1_final
        df["NPH_1_expected"] = NPH_list_1_expected
        df["ratio_1"] = ratio_list_1

        df["NPH_2_final"] = NPH_list_2_final
        df["NPH_2_expected"] = NPH_list_2_expected
        df["ratio_2"] = ratio_list_2

        df["CNP_1"] = CNP_list_1
        df.to_excel("PIC NORMALIZED.xlsx")


def is_in_range(self):
    if self > 1:
        self = 1
    if self < 0:
        self = 0
    return self
