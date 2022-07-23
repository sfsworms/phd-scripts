# -*- coding: utf-8 -*-
"""
Created on Thu Jul  7 15:35:37 2022

@author: worms
"""

import sys
 
if len(sys.argv) == 1:
    print('Please provide a command line argument as a file name!')
    sys.exit()
else:
    myfile = sys.argv[1]
 
sequence = ''
 
try:
    with open(myfile, 'r') as f:
        for line in f:
            if '>' not in line:
                sequence = sequence + line.strip()
except OSError as oserr:
    print('OS error:', oserr)
else:
    print(sequence)
 