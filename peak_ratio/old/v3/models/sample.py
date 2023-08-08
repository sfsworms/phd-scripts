import pandas as pd
import xlsxwriter
import glob


SETUP_FILENAME = "samples.xlsx"
WIDTH_COLUMN = 30

SAMPLE_ABI_FILES = "SAMPLE_ABI_FILES"
SAMPLE_TEMPLATE_1 = "FASTA_1"
SAMPLE_TEMPLATE_2 = "FASTA_2"
SAMPLE_STANDARD_1 = "STANDARD_1"
SAMPLE_STANDARD_2 = "STANDARD_2"

SAMPLING_CHANNELS = [SAMPLE_ABI_FILES, SAMPLE_TEMPLATE_1, SAMPLE_TEMPLATE_2, SAMPLE_STANDARD_1, SAMPLE_STANDARD_2]

# HEADERS de l'objet Samples
READ_COLUMN = "READ"
TEMPLATE_1_COLUMN = "TEMPLATE_1"
TEMPLATE_2_COLUMN = "TEMPLATE_2"
STANDARD_1_COLUMN = "STANDARD_1"
STANDARD_2_COLUMN = "STANDARD_2"
SAMPLE_COLUMNS = [READ_COLUMN, TEMPLATE_1_COLUMN, TEMPLATE_2_COLUMN, STANDARD_1_COLUMN, STANDARD_2_COLUMN]


class Sample:  # TODO à repenser => mettre dans la classe samples sous forme d'objet?
    def __init__(self, read, template_1, template_2, standard_1, standard_2):
        self.read = read
        self.template_1 = template_1
        self.template_2 = template_2
        self.standard_1 = standard_1
        self.standard_2 = standard_2
        self.sample_locations = []

    def set_strands(self, aligner):
        self.read.set_strand(aligner=aligner, template=self.template_1)
        self.standard_1.set_strand(aligner=aligner, template=self.template_1)
        self.standard_2.set_strand(aligner=aligner, template=self.template_2)

    def extract_traces(self):
        self.read.extract_trace()
        self.standard_1.extract_trace()
        self.standard_2.extract_trace()

    def set_length_sequences(self):
        self.read.set_length_sequence()
        self.standard_1.set_length_sequence()
        self.standard_2.set_length_sequence()

    def set_alignments(self, aligner):
        self.read.set_alignment(aligner=aligner, template=self.template_1)
        self.standard_1.set_alignment(aligner=aligner, template=self.template_1)
        self.standard_2.set_alignment(aligner=aligner, template=self.template_2)

    def set_sample_locations(self):
        # Ne fais rien si des valeurs sont présentes
        if len(self.sample_locations) != 0:
            return
        # self.mutation_locations = mutations entre les templates
        for index in range(self.template_1.length):
            if self.template_1.sequence[index] != self.template_2.sequence[index]:
                self.sample_locations.append(index)


class Samples(pd.DataFrame):
    def __init__(self):
        super().__init__()
        # self[SAMPLES_COLUMN] = sampling.get_samples_from_sampling(reads=reads,
        #                                                           standards=standards,
        #                                                           templates=templates)
        # self.add_variants(sampling=sampling)
        # self.add_conditions(sampling=sampling)

    # def add_variants(self, sampling):
    #     self[VARIANT_1_COLUMN] = sampling.df[SAMPLE_TEMPLATE_1]
    #     self[VARIANT_2_COLUMN] = sampling.df[SAMPLE_TEMPLATE_2]
    #
    # def add_conditions(self, sampling):
    #     for condition in sampling.get_condition_names():
    #         self[condition] = sampling.df[condition]

    def set_reads(self, aligner):
        for sample in self[SAMPLES_COLUMN]:
            sample.set_strands(aligner=aligner)
            sample.extract_traces()
            sample.set_length_sequences()
            sample.set_alignments(aligner=aligner)

    def set_locations(self):
        for sample in self[SAMPLES_COLUMN]:
            sample.set_sample_locations()

    def verify_standard_directionality(self):
        """Vérifie que la direction des reads des standards sont dans le même sens que le read de séquençage."""
        for sample in self[SAMPLES_COLUMN]:
            if sample.standard_1.strand != sample.read.strand or sample.standard_2.strand != sample.read.strand:
                return True
        return False

    def verify_template_directionality(self, aligner):
        """Vérifie que les templates sont dans le même sens."""
        for sample in self[SAMPLES_COLUMN]:
            if aligner.get_strand(sample.template_1.sequence, sample.template_2.sequence) == "-":
                return True
        return False

    def verify_template_length(self):
        """Vérifie que les templates ont la même longueur."""
        for sample in self[SAMPLES_COLUMN]:
            if sample.template_1.length != sample.template_2.length:
                return True
        return False


class MetadataReader:
    def __init__(self, reads, templates, standards):
        self.excel_filename = SETUP_FILENAME
        self.reads = reads
        self.templates = templates
        self.standards = standards

    def excel_file_exist(self):
        if glob.glob(self.excel_filename):
            return True

    def create_setup_file(self):  # TODO fonctionnel mais à formater
        workbook = xlsxwriter.Workbook(self.excel_filename)
        worksheet = workbook.add_worksheet()
        # Add a bold format to use to highlight cells.
        bold = workbook.add_format({'bold': True})

        for column, header in enumerate(SAMPLING_CHANNELS):
            # SAMPLING_CHANNELS = headers
            worksheet.write(0, column, header, bold)

            # SAMPLE_ABI_FILES
            if header == SAMPLE_ABI_FILES:
                for row, abi_file in enumerate(self.reads.files):
                    worksheet.write(row+1, column, abi_file)

            # SAMPLE_TEMPLATE_1 or SAMPLE_TEMPLATE_2
            if header == SAMPLE_TEMPLATE_1 or header == SAMPLE_TEMPLATE_2:
                worksheet.data_validation(1, column, len(self.reads.files), column,
                                          {'validate': 'list', 'source': self.templates.names})

            # SAMPLE_STANDARD_1 or SAMPLE_STANDARD_2
            if header == SAMPLE_STANDARD_1 or header == SAMPLE_STANDARD_2:
                worksheet.data_validation(1, column, len(self.reads.files), column,
                                          {'validate': 'list', 'source': self.standards.files})

        worksheet.set_column(first_col=0, last_col=len(SAMPLING_CHANNELS)-1, width=WIDTH_COLUMN)
        workbook.close()

    def loading_samples_from_excel(self):
        samples = pd.DataFrame()
        metadata = pd.read_excel(io=self.excel_filename)
        # Charge les objets
        samples[READ_COLUMN] = metadata[SAMPLE_ABI_FILES].apply(self.read_from_filename, 0)
        samples[TEMPLATE_1_COLUMN] = metadata[SAMPLE_TEMPLATE_1].apply(self.template_from_name, 0)
        samples[TEMPLATE_2_COLUMN] = metadata[SAMPLE_TEMPLATE_2].apply(self.template_from_name, 0)
        samples[STANDARD_1_COLUMN] = metadata[SAMPLE_STANDARD_1].apply(self.standard_from_filename, 0)
        samples[STANDARD_2_COLUMN] = metadata[SAMPLE_STANDARD_2].apply(self.standard_from_filename, 0)
        # Charge les options
        for header in metadata.columns.tolist():
            if header not in SAMPLING_CHANNELS:
                samples[header] = metadata[header]
        return samples

    def read_from_filename(self, filename):
        return list(filter(lambda obj: obj.filename == filename, self.reads))[0]

    def template_from_name(self, name):
        return list(filter(lambda obj: obj.name == name, self.templates))[0]

    def standard_from_filename(self, filename):
        return list(filter(lambda obj: obj.filename == filename, self.standards))[0]
