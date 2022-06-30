# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 16:27:46 2022

@author: worms
"""


import sys, math
from Bio import SeqIO

def batch_iterator(iterator, batch_size):
    """
    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.
    (derived from https://biopython.org/wiki/Split_large_file)
    """
    entry = True  # Make sure we loop once
    while entry:
        batch = []
        while len(batch) < batch_size:
            try:
                entry = iterator.__next__()
            except StopIteration:
                entry = None
            if entry is None:
                break # EOF = end of file
            batch.append(entry)
        if batch:
            yield batch

if(len(sys.argv) != 3):
        sys.exit("usage: split_fasta.py MULTI_FASTA_FILE N_CHUNKS")

ffile=sys.argv[1]  # fasta file
chunks=sys.argv[2] # number of chunks
nseq = len([1 for line in open(ffile) if line.startswith(">")])
chunksize=math.ceil(nseq/int(chunks))
print("Splitting multi-fasta file of", nseq, "sequences into chunks of size", chunksize)

records = SeqIO.parse(open(ffile), "fasta")
for i, batch in enumerate(batch_iterator(records, chunksize)):
        filename = "chunk_%i.fasta" % (i + 1)
        with open(filename, "w") as handle:
                count = SeqIO.write(batch, handle, "fasta")
        print("Wrote %i sequences to %s" % (count, filename))

sys.exit("Done.")