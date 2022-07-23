# -*- coding: utf-8 -*-
"""
Created on Mon Jul 11 11:25:33 2022

@author: worms
"""
import sys
import os

folder = sys.argv[1]

test = os.listdir(folder)

print(test)