import Bio
import sys

FILE_TYPE = "fasta"


# Check the length of sys.argv to see if arguments have been provided
if len(sys.argv) == 3:
    # Use values provided by the user
    file_to_split = sys.argv[1]
    batch_size = sys.argv[2]

if len(sys.argv) == 2:
    # Assign default values for missing arguments
    file_to_split = sys.argv[1] # Nom du fichier qu'on veut split, qui doit être dans le même répertoire
    batch_size = 10^6 # nb de peptides par fichier
    print("Size of batches unspecified, running with 1000 000 sequence per file.")

else:
    print("Please provide the name of the file to split as first argument and the size of batches as second argument, for exemple:")
    print("python fasta_splitter.py big_batch.fa 100000")
    sys.exit()


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
    record_iter = Bio.SeqIO.parse(open(file_to_split), FILE_TYPE)
    for i, batch in enumerate(batch_iterator(record_iter, batch_size)):
        filename = "group_%i.fasta" % (i + 1)
        with open(filename, "w") as handle:
            count = Bio.SeqIO.write(batch, handle, "fasta")
        print("Wrote %i records to %s" % (count, filename))

split_data()