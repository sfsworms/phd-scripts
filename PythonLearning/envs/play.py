# -*- coding: utf-8 -*-
"""
Created on Wed Jul  6 16:31:33 2022

@author: worms
"""
import math
import operator

testExponent = 54e6
print(testExponent)

string1 = 'This'
string2 = "something"
print(string1 + 'is' + string2)

string3 = string1 + 'is' + string2
string3 = string3.upper()
print(string3)


#List just have multiple variables
thisIsList = ["XYZ", "XAA"]
print(thisIsList)

#You can modify the
thisIsList[1] = "PLOP"
print(thisIsList)

#You can add elements
thisIsList.append("couocu")
print(thisIsList)

testString = ''.join(thisIsList)
print(testString)

listAgain = list(testString)
print(listAgain)
partOfList = listAgain[4:6]
print(partOfList)

partOfList.insert(1, "Prout")
print(partOfList)

partOfList.pop(1)
print(partOfList)

### Tuples
#Tuples are basically immutable lists
#They can be declared without parenthesis, but () are good practice
Histidine = 'H','CAT','CAC'
Lysine = ('K', 'AAA', 'AAG')

print(Histidine)

#They can be assigned to variables

item1, item2, item3 = Lysine
print(item1)
print(item2)
print(item3)

### Dictionary
#They are kinds of lists but where each entry is associated to a 'key'
#They are decalred with {}

myDict = {'Axe':'Cutting',
          'Sword':'Cutting',
          'Spear':'Piercing'}

print(myDict)

keys2 = myDict.keys()
keys3 = list(keys2)

values = list(myDict.values())
print(values)

myKey = 'Axe'

check = myKey in myDict
print(check)

#You can look up valure from keys

myDict['Axe']


###Exercice
x = 0.1 * 10.0
y = 0.1 + 0.1 + 0.1 + 0.1 + 0.1 + 0.1 + 0.1 + 0.1 + 0.1 + 0.1 
y
x

string = "gggtgcgacg attcattgtt ttcggacaag tggataggca accactaccg gtggattgtc"
print(string.replace(" ",""))
stringNoSpace = string.replace(" ","")
stringNoSpace.upper()

## Exercie


mydictionary = {"c" : "2", 
                "b" : "3", 
                "a" : "4", 
                "d" : "1"
                }
 
mydictionary_sorted_by_keys = sorted(mydictionary.items(), key=operator.itemgetter(0)) #Sorted works on all kinds of lists
mydictionary_sorted_by_values = sorted(mydictionary.items(), key=operator.itemgetter(1))
sorted(string)
print(mydictionary_sorted_by_keys)

len(mydictionary)

## Find operator

mouse = "mpltncrmvp varplslllt fflcacaetp prftrtpvdq tgvsggvasf icqatgdprp"
mouse = mouse.replace(" ","") 

pattern = "f"

pattern in mouse

mouse.count(pattern)

mouse.find(pattern)

secondF = mouse.find(pattern, 20+ len(pattern))
print(secondF)

thirdF = mouse.find(pattern,secondF+len(pattern))
print(thirdF)

## Slicing facility

codon = 'gac'

indexCodon = string.find(codon)
print(indexCodon)
nextCodon = string[indexCodon+3:indexCodon+6]

print(nextCodon)

## TRanslate lets you change file

dnaSeq = stringNoSpace.upper()

complement = dnaSeq.maketrans('acgtACGT', 'tgcaTGCA')

compSeq = dnaSeq.translate(complement) #Dictionary has to be made with a call to the string
print(complement)
print(compSeq)

dnaSeq2 = "GATAGGTGTGTAGAG"
compSeq2 = dnaSeq2.translate(complement) #But can be used on ither strings
compSeq2

bond_line = "|"*len(dnaSeq2)

print(dnaSeq2)
print(bond_line)
print(compSeq2)

compSeq2[::-1] #This reverse the string
