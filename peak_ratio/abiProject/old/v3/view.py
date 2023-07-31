def file_not_found(directory):
    print(f'No files in "{directory}" directory.')


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


def configuration_file_not_found(filename):
    input(f'Please fill the {filename} file and press enter to continue.')


def configuration_file_found(filename):
    print(f'The configuration file "{filename}" has been found.')


def loading_samples():
    print("Loading samples from excel file ...")


def end():
    print("Done")
