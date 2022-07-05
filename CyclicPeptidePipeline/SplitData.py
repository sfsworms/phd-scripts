from Bio import SeqIO

FILE_TO_SPLIT = "Cytoplasmic-NNK-Gen-1-LB_R1_001_peptide3.fasta" # Nom du fichier qu'on veut split, qui doit être dans le même répertoire
FILE_TYPE = "fasta"
BATCH_SIZE = 1000000  # nb de peptides par fichier


def batch_iterator(iterator, batch_size):
    """Returns lists of length batch_size.

    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.AlignIO.parse(...), or simply
    lines from a file handle.

    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.
    """
    entry = True  # Make sure we loop once
    it = iter(iterator)
    while entry:
        batch = []
        while len(batch) < batch_size:
            try:
                entry = next(it)
            except StopIteration:
                entry = None
            if entry is None:
                # End of file
                break
            batch.append(entry)
        if batch:
            yield batch


def split_data():
    record_iter = SeqIO.parse(open(FILE_TO_SPLIT), FILE_TYPE)
    for i, batch in enumerate(batch_iterator(record_iter, BATCH_SIZE)):
        filename = "group_%i.fasta" % (i + 1)
        with open(filename, "w") as handle:
            count = SeqIO.write(batch, handle, "fasta")
        print("Wrote %i records to %s" % (count, filename))

split_data()