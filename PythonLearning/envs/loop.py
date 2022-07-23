# -*- coding: utf-8 -*-
"""
Created on Thu Jul  7 14:09:14 2022

@author: worms
"""

### We're going to play with loops. Indents are important.

x = 3


if(x == 2):
    print("If was TRUE")
elif(x == 3):
    print("x is true")
else:
    print("I do not know x")
    
## While loop

i = 0
while(i < 10):
    i = i+1
    print(i)
    
## For loop

bases = ['T' , 'C', 'G', 'A']

for base in bases:
    print(base)
    
## range()

population = [425]
r = 0.0197
numbersOfYears=30
years = range(0,numbersOfYears+1, 1)
print(list(years))

for year in years:
    population.append(population[year]*(r+1))
    
for year in years:
    print("At years %2d, the population is %7d" %(year, population[year]))
    
## Dictionary count

dnaSeq = 'AGAGCTGGTGCACAGTCGTAGCCAGCA'

codonNumbers = int(len(dnaSeq)/3)

codonList = []

for i in range(codonNumbers):
    codonList[i] = 