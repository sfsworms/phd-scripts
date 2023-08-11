def error_templates():
    print("An unexpected error has been detected. Please check the reading direction of your templates in the fasta "
          "file.")


def error_length_templates():
    print("An unexpected error has been detected: Templates must be exactly the same length!")


def error_standard():
    print("An unexpected error has been detected: Please make sure the .ab1 file of the read is in the same direction as the reads of the standards ab1.")


def show_settings_read(name, strand):
    print(f'{name} is done ({strand}).')


def show_samples(df):
    print("List of samples being treated:")
    print(df.to_string())


def fill_sample_file(file):
    input(f'Please fill the {file} file and press enter to continue.')
    
    
def error_read_quality(file):
    print(f"{file} is of low quality (Only 'N'), please check the ab1 file and remove if sequencing failed.")
