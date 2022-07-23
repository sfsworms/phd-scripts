# -*- coding: utf-8 -*-
"""
Created on Thu Jul  7 15:35:37 2022

@author: worms
"""
from Bio import SeqIO
record = SeqIO.read("NC_000913.gbk", "genbank")
output_handle = open("NC_000913_cds.fasta", "w")
count = 0
for feature in record.features:
    if feature.type == "CDS":
        count = count + 1
        feature_name = feature.qualifiers["locus_tag"][0]
        feature_start = feature.location.start
        feature_end = feature.location.end
        feature_seq = feature.extract(record.seq)
        # Simple FASTA output without line wrapping:
        output_handle.write(">gi|556503834|ref|NC_000913.3|:"+ str(feature_start) + "-" + str(feature_end)+" Escherichia coli str. K-12 substr. MG1655, complete genome" + "\n" + str(feature_seq) + "\n")
output_handle.close()
print(str(count) + " CDS sequences extracted")

test = record.features[345]

print(test.location.start)

print(record.id)
