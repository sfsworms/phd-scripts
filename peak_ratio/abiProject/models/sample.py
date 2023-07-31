import pandas as pd
import xlsxwriter
import glob


EXCEL_FILENAME = "samples.xlsx"
WIDTH_COLUMN = 30

SAMPLE_ABI_FILES = "SAMPLE_ABI_FILES"
SAMPLE_TEMPLATE_1 = "FASTA_1"
SAMPLE_TEMPLATE_2 = "FASTA_2"
SAMPLE_STANDARD_1 = "STANDARD_1"
SAMPLE_STANDARD_2 = "STANDARD_2"

SAMPLING_CHANNELS = [SAMPLE_ABI_FILES, SAMPLE_TEMPLATE_1, SAMPLE_TEMPLATE_2, SAMPLE_STANDARD_1, SAMPLE_STANDARD_2]

SAMPLES_COLUMN = "SAMPLES"
VARIANT_1_COLUMN = "VARIANT 1"
VARIANT_2_COLUMN = "VARIANT 2"


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


class Sampling:
    def __init__(self):
        self.filename = EXCEL_FILENAME
        self.df = pd.read_excel(io=self.filename)
        self.samples_files = self.get_sample_files()

    def get_sample_files(self):
        return self.df[SAMPLE_ABI_FILES].tolist()

    def get_value(self, sample, channel):
        return self.df[self.df[SAMPLE_ABI_FILES] == sample][channel].values[0]

    def get_read(self, sample, reads):
        filename = self.get_value(sample=sample, channel=SAMPLE_ABI_FILES)
        return list(filter(lambda obj: obj.filename == filename, reads))[0]

    def get_template_1(self, sample, templates):
        name = self.get_value(sample=sample, channel=SAMPLE_TEMPLATE_1)
        return list(filter(lambda obj: obj.name == name, templates))[0]

    def get_template_2(self, sample, templates):
        name = self.get_value(sample=sample, channel=SAMPLE_TEMPLATE_2)
        return list(filter(lambda obj: obj.name == name, templates))[0]

    def get_standard_1(self, sample, standards):
        filename = self.get_value(sample=sample, channel=SAMPLE_STANDARD_1)
        return list(filter(lambda obj: obj.filename == filename, standards))[0]

    def get_standard_2(self, sample, standards):
        filename = self.get_value(sample=sample, channel=SAMPLE_STANDARD_2)
        return list(filter(lambda obj: obj.filename == filename, standards))[0]

    def get_samples_from_sampling(self, reads, standards, templates):
        samples = []
        for sample in self.samples_files:
            read = self.get_read(sample=sample, reads=reads)
            template_1 = self.get_template_1(sample=sample, templates=templates)
            template_2 = self.get_template_2(sample=sample, templates=templates)
            standard_1 = self.get_standard_1(sample=sample, standards=standards)
            standard_2 = self.get_standard_2(sample=sample, standards=standards)
            samples.append(Sample(read=read,
                                  template_1=template_1,
                                  template_2=template_2,
                                  standard_1=standard_1,
                                  standard_2=standard_2))
        return samples

    def get_condition_names(self):
        headers = self.df.columns.tolist()
        for channel in SAMPLING_CHANNELS:
            headers.pop(headers.index(channel))
        return headers


class Samples(pd.DataFrame):
    def __init__(self, sampling: Sampling, reads, standards, templates):
        super().__init__()
        self[SAMPLES_COLUMN] = sampling.get_samples_from_sampling(reads=reads,
                                                                  standards=standards,
                                                                  templates=templates)
        self.add_variants(sampling=sampling)
        self.add_conditions(sampling=sampling)

    def add_variants(self, sampling):
        self[VARIANT_1_COLUMN] = sampling.df[SAMPLE_TEMPLATE_1]
        self[VARIANT_2_COLUMN] = sampling.df[SAMPLE_TEMPLATE_2]

    def add_conditions(self, sampling):
        for condition in sampling.get_condition_names():
            self[condition] = sampling.df[condition]

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


class Excel:
    def __init__(self, reads, templates, standards):
        self.filename = EXCEL_FILENAME
        self.abi_files = reads.files
        self.template_names = templates.names
        self.standards_files = standards.files

    def excel_file_exist(self):
        if not glob.glob(self.filename):
            return None
        else:
            return True

    def create(self):  # TODO fonctionnel mais à formater
        workbook = xlsxwriter.Workbook(self.filename)
        worksheet = workbook.add_worksheet()
        # Add a bold format to use to highlight cells.
        bold = workbook.add_format({'bold': True})

        for column, header in enumerate(SAMPLING_CHANNELS):
            # SAMPLING_CHANNELS = headers
            worksheet.write(0, column, header, bold)

            # SAMPLE_ABI_FILES
            if header == SAMPLE_ABI_FILES:
                for row, abi_file in enumerate(self.abi_files):
                    worksheet.write(row+1, column, abi_file)

            # SAMPLE_TEMPLATE_1 or SAMPLE_TEMPLATE_2
            if header == SAMPLE_TEMPLATE_1 or header == SAMPLE_TEMPLATE_2:
                worksheet.data_validation(1, column, len(self.abi_files), column,
                                          {'validate': 'list', 'source': self.template_names})

            # SAMPLE_STANDARD_1 or SAMPLE_STANDARD_2
            if header == SAMPLE_STANDARD_1 or header == SAMPLE_STANDARD_2:
                worksheet.data_validation(1, column, len(self.abi_files), column,
                                          {'validate': 'list', 'source': self.standards_files})

        worksheet.set_column(first_col=0, last_col=len(SAMPLING_CHANNELS)-1, width=WIDTH_COLUMN)
        workbook.close()
