import pandas as pd
import view
import models.plot as plot
from models.abi import Reads
from models.fasta import Templates
from models.aligner import LocalAlignment
from models.sample import Sampling, Samples, Excel #Doesn't import create?

class Controller:
    def __init__(self):
        # This flag indicate whether the controller is actively running or not.
        self.running = True 
        
        # Initializing lists for reads and standards from ABI files and templates from FASTA files.
        self.reads: list = Reads(root_dir="abi") #Get the ab1 files name of the reads in "abi"
        self.standards: list = Reads(root_dir="standards") # Get the ab1 files name of the standards in "standards". They're needed because the height of a peak vary based on the local environment
        self.templates: list = Templates(root_dir="fasta") # Get the templates sequences from the fasta. I they're used to fix the position on the ab1 files

        # Check for an Excel file that has sequence info, or prompt the user to create on (see def block)
        # It gives an excel file with samples (ab1), fasta1, fasta2, standard1 and standard2
        self.create_excel_file()
        
        # Create objects for sampling, alignment, and samples.
        self.sampling = Sampling() #Creates a Sampling() class object with the values for all the ab1 file. As I understand it, this wouldn't call the functions outside init yet
        self.aligner = LocalAlignment() # Create an aligner for the controller to use
        self.samples = Samples(sampling=self.sampling,
                               reads=self.reads,
                               standards=self.standards,
                               templates=self.templates) #Store all the previous values in a Samples df.
        
        self.verify_reads_quality()

        # Set reads in samples and verify their direction.
        self.samples.set_reads(self.aligner) # This use the set reads method of samples to align reads and return the alignment on the templates, trace file, length of sequence and other, for each file in the SAMPLES colums
        self.verify_reading_direction() # Runs a series of check in "verify reading direction"
        self.samples.set_locations() # For each sample, grab the two corresponding templates and identify all the mutations positions

    def create_excel_file(self):
        # Create an Excel file if it doesn't already exist. That excel file should contains the templates and standards for all reads and are then used in the script
        excel = Excel(reads=self.reads, templates=self.templates, standards=self.standards)
        if excel.excel_file_exist():
            return
        excel.create() # This print out the excel file it it exists. At that point it usually would be empty.
        view.fill_sample_file(file=excel.filename)
    
    #Check the quality of the reads. If they are only N calls an error
    # Somehow doesn't stop the execution
    def verify_reads_quality(self):
        problematic_file = self.samples.verify_only_n()
        if problematic_file:
            view.error_read_quality(problematic_file)
            self.running = False
            exit()

    def verify_reading_direction(self):
        # Verify the directionality of standard reads, template reads, and template lengths.
        # If there's an issue, an error is shown via the view module and the running flag is set to False.
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
        view.show_samples(df=self.samples) # Print the samples to the console

        for sample in self.samples["SAMPLES"]:
            plot.chromatogram_read(sample) # Plot the chromatograms around the selected bit for all samples. That's in model.plot. All the outputs should be png. 

        # Extract and save peak values for the samples to an Excel file.
        # Repeated logic for sample reads, standard_1, and standard_2.
        # These blocks capture peak values based on the sample's location and reference.
        
        read_list = []
        pic1_list = []
        pic2_list = []
        reference_list = []
        df = pd.DataFrame()
        for sample in self.samples["SAMPLES"]: #Look at the various samples. Just take the "samples" column.
           
            # Define a normalization distance constant. This is how many peaks before the barcoding peak the reference peak is.
            NORMALIZATION_DISTANCE = 5
            
            # Determine the reference location based on the strand of the read.
            # If the strand is positive, use the first location minus the normalization distance.
            # If the strand is negative, use the last location plus the normalization distance.
    
            if sample.read.strand == "+":
                reference_location = sample.sample_locations[0] - NORMALIZATION_DISTANCE
            else:
                reference_location = sample.sample_locations[-1] + NORMALIZATION_DISTANCE

            # For each sample location, determine the chromatogram keys for both templates and the reference.
            for sample_location in sample.sample_locations:
                key_1 = sample.template_1.sequence[sample_location]
                key_2 = sample.template_2.sequence[sample_location]
                key_reference = sample.template_1.sequence[reference_location]

                # Fetch the chromatograms for the determined keys.
                chromatogram_1 = sample.read.chromatograms[key_1]
                chromatogram_2 = sample.read.chromatograms[key_2]
                chromatogram_reference = sample.read.chromatograms[key_reference]
        
                # Determine the peak values based on the sample location.
                # TODO: Consider using a list for handling multiple mutations.
                pic_value_1 = chromatogram_1[sample.read.base_location[sample.read.get_locations([sample_location])[0]]]
                pic_value_2 = chromatogram_2[sample.read.base_location[sample.read.get_locations([sample_location])[0]]]
                pic_reference = chromatogram_reference[sample.read.base_location[sample.read.get_locations([reference_location])[0]]]

                read_list.append(sample.read.name)
                pic1_list.append(pic_value_1)
                pic2_list.append(pic_value_2)
                reference_list.append(pic_reference)
                    
        # Add all the values to the df before printing it out.
        df["SAMPLE"] = read_list # Names of the samples
        df["PIC_1"] = pic1_list # Values of the peak fitting to the first reference
        df["PIC_2"] = pic2_list # Values of the peak fitting to the first reference
        df["REFERENCE +5"] = reference_list # Value of the reference?
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

                # Obtention des chromatograms des standards 1
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