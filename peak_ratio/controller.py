import pandas as pd
import view
import models.plot as plot
from models.abi import Reads
from models.fasta import templates_list
from models.aligner import LocalAlignment
from models.sample import Sampling, Samples, Excel #Doesn't import create?


# Define a normalization distance constant. This is how many peaks before the barcoding peak the reference peak is.
NORMALIZATION_DISTANCE = 5 

class Controller:
    def __init__(self):
        # This flag indicate whether the controller is actively running or not.
        self.running = True 
        
        # Initializing lists for reads and standards from ABI files and templates from FASTA files.
        self.reads: list = Reads(root_dir="abi") #Get the ab1 files name of the reads in "abi"
        self.standards: list = Reads(root_dir="standards") # Get the ab1 files name of the standards in "standards". They're needed because the height of a peak vary based on the local environment
        self.templates_list: list = templates_list(root_dir="fasta") # Get the templates sequences from the fasta. They're used to fix the position on the ab1 files

        # Check for an Excel file that has sequence info, or prompt the user to create on (see def block)
        # It gives an excel file with samples (ab1), fasta1, fasta2, standard1 and standard2
        self.create_excel_file()
        
        # Create objects for sampling, alignment, and samples.
        self.sampling = Sampling() #Creates a Sampling() class object with the values for all the ab1 file. As I understand it, this wouldn't call the functions outside init yet
        self.aligner = LocalAlignment() # Create an aligner for the controller to use
        self.samples = Samples(sampling=self.sampling,
                               reads=self.reads,
                               standards=self.standards,
                               templates=self.templates_list) #Store all the previous values in a Samples df, which has a column of sample objects, and then two columns called VARIANTS which contains the template names
        
        self.verify_reads_quality()

        # Set reads in samples and verify their direction.
        self.samples.set_reads(self.aligner) # This use the set reads method of samples to align reads and return the alignment on the templates, trace file, length of sequence and other, for each file in the SAMPLES colums
        self.verify_reading_direction() # Runs a series of check in "verify reading direction"
        self.samples.set_locations() # For each sample, grab the two corresponding templates and identify all the mutations positions

    def create_excel_file(self):
        # Create an Excel file if it doesn't already exist. That excel file should contains the templates and standards for all reads and are then used in the script
        excel = Excel(reads=self.reads, templates=self.templates_list, standards=self.standards)
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
        
    # This methods should create a single DF that include the values for all peaks for the read 
    #TODO: and standards in the same DF
    # Equally, this should add the key of the value of each standards
    #TODO: Get a way to have more than 2 standards
    
    def get_peaks_read(self):
        # Define a function to get the peak value for a given nucleotide at a specific location
        def get_peak_value(chromatogram, read, location):
            return chromatogram[read.base_location[read.get_locations_on_read([location])[0]]]
        
        all_rows = []  # List to hold all rows before converting to DataFrame
    
        # Iterate over each sample
        for sample in self.samples["SAMPLES"]:
            # Determine the reference location and column based on the strand of the read
            if sample.standard_2.strand == "+":
                reference_location = sample.sample_locations[0] - NORMALIZATION_DISTANCE
                REFERENCE_COLUMN = f"Reference (+{NORMALIZATION_DISTANCE})"
            else:
                reference_location = sample.sample_locations[-1] + NORMALIZATION_DISTANCE
                REFERENCE_COLUMN = f"Reference (-{NORMALIZATION_DISTANCE})"
            
            # Get base and chromatogram of the reference
            base_reference = sample.template_1.sequence[reference_location]
            chromatogram_reference = sample.read.chromatograms[base_reference]
            
            # Fetch chromatograms for each nucleotide
            chromatograms = {
                "A": sample.read.chromatograms["A"],
                "C": sample.read.chromatograms["C"],
                "G": sample.read.chromatograms["G"],
                "T": sample.read.chromatograms["T"]
            }
            
            # For each mutation location in the sample
            for sample_location in sample.sample_locations:
                row_data = {
                    "SAMPLE": sample.read.name,
                    "POSITION (FASTA)": sample_location,
                    "BASE TEMPLATE 1": sample.template_1.sequence[sample_location],
                    "BASE TEMPLATE 2": sample.template_2.sequence[sample_location],
                    REFERENCE_COLUMN: get_peak_value(chromatogram_reference, sample.read, reference_location)
                }
                
                # Get peak values for each nucleotide
                for nucleotide in ["A", "C", "G", "T"]:
                    row_data[f"PEAK_{nucleotide}"] = get_peak_value(chromatograms[nucleotide], sample.read, sample_location)
                
                # Append the row data to the all_rows list
                all_rows.append(row_data)
        
        # Convert the all_rows list to a DataFrame using concat
        df = pd.concat([pd.DataFrame([row]) for row in all_rows], ignore_index=True)
        
        return df

# NOTE: This refactored method assumes that 'self.samples' is a DataFrame or a structure that supports dictionary-like indexing.
# It's based on the provided code, but might need adjustments based on the complete program context.

        view.show_samples(df=self.samples) # Print the samples to the console

        for sample in self.samples["SAMPLES"]:
            plot.chromatogram_read(sample) # Plot the chromatograms around the selected bit for all samples. That's in model.plot. All the outputs should be png. 

        # Extract and save peak values for the samples to an Excel file.
        # Repeated logic for sample reads, standard_1, and standard_2.
        # These blocks capture peak values based on the sample's location and reference.
        
        #Create dfs with all the calls for the read, template 1 and template 2
        
        df_reads = self.get_peaks_read()

        #TODO: Edit the templates bit, export it to the same excel file
        # Merge the df
        #     #Make a list of the dfs to merge
        # dfs = [df_reads, df_template1, df_template2]
        # results_df = reduce(lambda left,right: pd.merge(left,right,on=["SAMPLE", "POSITION_FASTA"]), dfs)

        # Send to excel
        df_reads.to_excel("PIC_SAMPLES_concat.xlsx")

        
        
        
        # This is the same logic, using the traces of the standards ab1 files.

        # STD_1
        read_list = []
        pic1_list = []
        pic2_list = []
        reference_list = []
        df = pd.DataFrame()
        for sample in self.samples["SAMPLES"]: # For each individual sample in the self.samples (here 1 i guess.)

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
                pic_value_1 = chromatogram_1[sample.standard_1.base_location[sample.standard_1.get_locations_on_read([sample_location])[0]]]
                pic_value_2 = chromatogram_2[sample.standard_1.base_location[sample.standard_1.get_locations_on_read([sample_location])[0]]]
                pic_reference = chromatogram_reference[
                    sample.standard_1.base_location[sample.standard_1.get_locations_on_read([reference_location])[0]]]

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
                pic_value_1 = chromatogram_1[sample.standard_2.base_location[sample.standard_2.get_locations_on_read([sample_location])[0]]]
                pic_value_2 = chromatogram_2[sample.standard_2.base_location[sample.standard_2.get_locations_on_read([sample_location])[0]]]
                pic_reference = chromatogram_reference[sample.standard_2.base_location[sample.standard_2.get_locations_on_read([reference_location])[0]]]

                read_list.append(sample.read.name)
                pic1_list.append(pic_value_1)
                pic2_list.append(pic_value_2)
                reference_list.append(pic_reference)

        df["SAMPLE"] = read_list
        df["PIC_1"] = pic1_list
        df["PIC_2"] = pic2_list
        df["REFERENCE +5"] = reference_list
        df.to_excel("PIC_STANDARD_2.xlsx")
        
        
        
        # This is a refactore version of get peaks. HAvent checked it yet.
        
 
        