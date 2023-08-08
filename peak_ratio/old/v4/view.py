def file_not_found(generic_file, directory):
    print(f'No {generic_file} file found in the {directory} directory.')


def error_templates():
    print("An unexpected error has been detected. Please check the reading direction of your templates in the fasta "
          "file.")


def error_length_templates():
    print("An unexpected error has been detected: Templates must be exactly the same length!")


def error_standard():
    print("An unexpected error has been detected: Please check the reading direction of your abi files.")


def show_settings_read(name, strand):
    print(f'{name} is done ({strand}).')


def show_samples(df):
    print(df.to_string())


def fill_sample_file(file):
    input(f'Please fill the {file} file and press enter to continue.')
