"""
Commented on 09/08/2023

Sample is a storage class to store the details associated with each samples ( template, standards...)
Sampling gets the values for the lines of the ab1 files names in an excel file
Samples store a seried of values in a pd data frame
Excel gets the data needed for outputting in an excel file
    create creates the excel files using the data in Excel

@author: worms
"""


import pandas as pd
import xlsxwriter
import glob

# Loads of constant used to name the files with the sample information

EXCEL_FILENAME = "samples.xlsx"
WIDTH_COLUMN = 30

# Constants defining column headers in the Excel file
SAMPLE_ABI_FILES = "SAMPLE_ABI_FILES"
SAMPLE_TEMPLATE_1 = "FASTA_1"
SAMPLE_TEMPLATE_2 = "FASTA_2"
SAMPLE_STANDARD_1 = "STANDARD_1"
SAMPLE_STANDARD_2 = "STANDARD_2"

SAMPLING_CHANNELS = [SAMPLE_ABI_FILES, SAMPLE_TEMPLATE_1, SAMPLE_TEMPLATE_2, SAMPLE_STANDARD_1, SAMPLE_STANDARD_2]

SAMPLES_COLUMN = "SAMPLES"
VARIANT_1_COLUMN = "VARIANT 1"
VARIANT_2_COLUMN = "VARIANT 2"

# Sample class encapsulates the details of each sample
class Sample:  # TODO à repenser => mettre dans la classe samples sous forme d'objet?
    # Storing the provided attributes for each sample
    def __init__(self, read, template_1, template_2, standard_1, standard_2):
        self.read = read
        self.template_1 = template_1
        self.template_2 = template_2
        self.standard_1 = standard_1
        self.standard_2 = standard_2
        # This intialize as an empty list, but it's supposed to contain the positions of the mutations relevant to the analysis of the sample
        self.sample_locations = []

    # Sets the strand direction for reads and standards based on alignment
    # Aligning all of those on the templates let us pick positions at the right place
    def set_strands(self, aligner):
        self.read.set_strand(aligner=aligner, template=self.template_1)
        self.standard_1.set_strand(aligner=aligner, template=self.template_1)
        self.standard_2.set_strand(aligner=aligner, template=self.template_2)

    # Extracts the traces for reads and standards
    def extract_traces(self):
        self.read.extract_trace()
        self.standard_1.extract_trace()
        self.standard_2.extract_trace()
    
    # Sets the length of sequences for reads and standards
    def set_length_sequences(self):
        self.read.set_length_sequence()
        self.standard_1.set_length_sequence()
        self.standard_2.set_length_sequence()

    # Aligns the reads and standards with their respective templates
    def set_alignments(self, aligner):
        self.read.set_alignment(aligner=aligner, template=self.template_1)
        self.standard_1.set_alignment(aligner=aligner, template=self.template_1)
        self.standard_2.set_alignment(aligner=aligner, template=self.template_2)

    # Identifies and stores the locations of mutations between templates
    def set_sample_locations(self):
        # Don't do anything if values already exist
        if len(self.sample_locations) != 0:
            print("Error : sample_locations didn't intialize as empty.")
            return
        # Identify mutation locations by comparing templates this gives a list of
        # mutation position
        for index in range(self.template_1.length): # For each base of the template
            if self.template_1.sequence[index] != self.template_2.sequence[index]: 
                #This fails!
                #If the base at that location are different, print the base
                
                self.sample_locations.append(index)
        # Check that some mutations were identified
        if len(self.sample_locations) == 0:
            print("Error: no difference were found between the templates")
            exit()

                
# Sampling class handles reading and extracting sample information from an Excel file
class Sampling:
    def __init__(self):
        self.filename = EXCEL_FILENAME
        # Read the Excel file into a pandas DataFrame
        self.df = pd.read_excel(io=self.filename)
        # Extract list of sample files
        self.samples_files = self.get_sample_files()

    # Extracts a list of sample files from the DataFrame from the excel column specified
    def get_sample_files(self):
        return self.df[SAMPLE_ABI_FILES].tolist()

    # Gets the values for each lane
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

    # Constructs a list of Sample objects based on the data read from the Excel file
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

    # Extracts condition names, which are headers not included in the predefined SAMPLING_CHANNELS list
    def get_condition_names(self):
        headers = self.df.columns.tolist()
        for channel in SAMPLING_CHANNELS:
            headers.pop(headers.index(channel))
        return headers


# Extends pandas DataFrame to create a specialized DataFrame for storing and processing samples
class Samples(pd.DataFrame): # Samples inherit from data frame class
    #sampling:Sampling means that sampling should be of Sampling class, but doesn't actually enforce it
    def __init__(self, sampling: Sampling, reads, standards, templates):
        super().__init__()
        
        # Populate the DataFrame with the output of the sampling.get_samples_from_sampling method and assign the results to the SAMPLES_COLUMN column
        self[SAMPLES_COLUMN] = sampling.get_samples_from_sampling(reads=reads,
                                                                  standards=standards,
                                                                  templates=templates)
        self.add_variants(sampling=sampling)
        self.add_conditions(sampling=sampling)
        
    # This just seems to create a variant columns that copies from "SAMPLE_TEMPLATE3
    def add_variants(self, sampling):
        self[VARIANT_1_COLUMN] = sampling.df[SAMPLE_TEMPLATE_1]
        self[VARIANT_2_COLUMN] = sampling.df[SAMPLE_TEMPLATE_2]

    #Add conditions to the DF if they are presets.
    def add_conditions(self, sampling):
        for condition in sampling.get_condition_names():
            self[condition] = sampling.df[condition]

    #For each sample in the list, gets the alignment, traces, sequence length and alignment
    def set_reads(self, aligner):
        for sample in self[SAMPLES_COLUMN]:
            sample.set_strands(aligner=aligner)
            sample.extract_traces()
            sample.set_length_sequences()
            sample.set_alignments(aligner=aligner)

    #Identify the locations
    def set_locations(self):
        for sample in self[SAMPLES_COLUMN]:
            sample.set_sample_locations()
            
    def verify_only_n(self):
        """Check that the reads doesn't only have N (aka sequencing failed)"""
        for sample in self[SAMPLES_COLUMN]:
            s = sample.read.get_sequence()
            # If all the letters in the sequence are N, return True
            if all(char == 'N' for char in s):
                return sample.read.filename

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

# Excel class handles creation and management of Excel files
class Excel:
    def __init__(self, reads, templates, standards):
        self.filename = EXCEL_FILENAME
        # Storing file and name lists for reads, templates, and standards
        self.abi_files = reads.files
        self.template_names = templates.names
        self.standards_files = standards.files

    # Checks if an Excel file already exists
    def excel_file_exist(self):
        # Returns None if no file found, True otherwise
        if not glob.glob(self.filename):
            return None
        else:
            return True
    # Creates a new Excel file with pre-defined structure
    def create(self):  # TODO fonctionnel mais à formater
        workbook = xlsxwriter.Workbook(self.filename) # Gets the name needed
        worksheet = workbook.add_worksheet()
        # Add a bold format to use to highlight cells.
        bold = workbook.add_format({'bold': True})

        for column, header in enumerate(SAMPLING_CHANNELS):
            # SAMPLING_CHANNELS = headers
            worksheet.write(0, column, header, bold)

            # SAMPLE_ABI_FILES Make a row for eahc abi files
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
