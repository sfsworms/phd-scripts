def fill_sample_file(file):
    input(f'Please fill the {file} file and press enter to continue.')


def get_template_index(id_records):
    sequences_length = len(id_records)

    # Affichage des séquences à sélectionner.
    for index, id_record in enumerate(id_records):
        print(f'{index}: {id_record}')

    # Boucle while qui s'arrête lorsque l'utilisateur donne une entrée valide.
    while True:
        back = input(f'Select the wild-type or template sequence (0 - {sequences_length - 1}): ')

        # Gestion de l'erreur au cas où l'entrée utilisateur n'est pas un nombre entier.
        try:
            wild_type_index = int(back)
        except ValueError:
            # Set une valeur impossible à atteindre.
            wild_type_index = -1

        # Arrête la boucle si l'entrée est valide et la renvoie, sinon print un message d'erreur.
        if wild_type_index in range(sequences_length):
            return wild_type_index
        else:
            print("That was no valid number.  Try again...")


def show_variants(variants):
    for variant in variants:
        print(f'{variant.id}: {variant.substitution}, {variant.codon_wt} -> {variant.codon_substituted}')
        # print(variant.mutation_index, variant.mutation_position, variant.nucleotide_wild_type, variant.nucleotide_substituted)

