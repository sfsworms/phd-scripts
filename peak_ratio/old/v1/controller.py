import model
import view


def get_barcodes(records, id_records):
    # Choix du template pour les alignements
    template_index = view.get_template_index(id_records=id_records)
    # Création des objets Template()
    template = model.Template(record=records[template_index])
    barcode_index = model.get_barcode_index(template=template, records=records)
    barcodes = model.get_barcodes(template=template, records=records, barcode_index=barcode_index)
    return template, barcodes


def get_variants(records, id_records):
    # Choix du wild-type
    wild_type_index = view.get_template_index(id_records=id_records)
    # Création des objets Template()
    wild_type = model.WildType(record=records[wild_type_index])
    variants = model.get_variants(wild_type=wild_type, records=records)
    view.show_variants(variants=variants)
    return wild_type, variants


def initialization():
    # Création du fichier sample.xlsx et demande de le remplir si celui-ci n'existait pas auparavant.
    if model.create_sample_file():
        view.fill_sample_file(model.SAMPLES_FILE)


def run():
    initialization()
    # Créer un objet d'alignement
    aligner = model.LocalAlignment()
    # Overture du fichier fasta
    records = model.open_fasta_file()
    # Liste des séquences fasta
    id_records = model.get_fasta_ids(records=records)

    template, variants = get_variants(records=records, id_records=id_records)
#    template, barcodes = get_barcodes(records=records, id_records=id_records)

    # Création des objets Read() au départ des fichiers "*.ab1"
    abi_files = model.get_abi_files()
    # TODO -> faire le lien entre les fichier abi et les conditions
    reads = model.get_reads(files=abi_files, template=template, aligner=aligner)

    for read in reads:  # TODO à finir/supprimer
        print(read.alignment)
        print(read.alignment.aligned)
        print(read.sequence)

        first_base_template = read.alignment.aligned[0][0][0]  # TODO: attention aux gaps!!!
        first_base_read = read.alignment.aligned[1][0][0]  # TODO: attention aux gaps!!!

        print(first_base_template)
        print(first_base_read)
        print(read.base_location)
        print(len(read.chromatograms["A"]))
        print(read.length)
        print(len(read.base_location))
#        print(read.sequence[first_base_read])


# TODO générale faire gaffe au type de retour de fonction (1 seul type ou None)
