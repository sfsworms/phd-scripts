# -*- coding: utf-8 -*-
"""
Created on Thu Jul  7 15:35:37 2022

@author: worms
"""
filename = "Testwriting.fa"


with open(filename, 'w') as f:
    print(">JP1 Partial Genome", file = f)
    f.write("ATG")
    f.write("CACACACACACADDADA")
    