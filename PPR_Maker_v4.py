#!/usr/bin/env python
# coding: utf-8

# In[110]:


# WELCOME

print('')
print('Welcome to the PPR Maker')
print('')
print('Loading...')
print('')

# LIBRARIES

import pandas as pd
import numpy as np
import datetime as dt

# FUNCTIONS

# Integer validation
def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        pass
    try:
        import unicodedata
        unicodedata.numeric(s)
        return True
    except (TypeError, ValueError):
        pass
    return False

# Find the coordinates of an item in a 2D list
def coordinates(myList, v):
    boolean = False
    for i, x in enumerate(myList):
        if v in x:
            boolean = True
            return (i, x.index(v))
    if boolean == False:
        print('Element not in list')

# INITIALIZATIONS

# RNA bases and PPR motifs
RNAbases = ['U','C','A','G']
PPRmotifs = ['ND','NS','SN','TD']

# All plasmids
allPlasmids = []
for i in range(0,6):
    for j in range(0,48):
        allPlasmids.append('R' + str(i+1) + '-' + str(j+1))
allPlasmidsX = [allPlasmids[i:i + 12] for i in range(0, len(allPlasmids), 12)]

# Block 1 (structuraly identical to blocks 3 and 5)
Block1 = [[],[],[],[]]
for i in range(0,48):
    Block1[0].append(i+1)
for i in range(0,12):
    Block1[1].append('SN')
for i in range(12,24):
    Block1[1].append('NS')
for i in range(24,36):
    Block1[1].append('TD')   
for i in range(36,48):
    Block1[1].append('ND')
for i in range(0,3):
    Block1[2].append('SN')
for i in range(3,6):
    Block1[2].append('NS')
for i in range(6,9):
    Block1[2].append('TD')
for i in range(9,12):
    Block1[2].append('ND')
for i in range(12,24):
    Block1[2].append(Block1[2][i-12])
for i in range(24,48):
    Block1[2].append(Block1[2][i-24])
Block1[3].append('S')
Block1[3].append('N')
Block1[3].append('T')
for i in range(3,48):
    Block1[3].append(Block1[3][i-3])
Block1X = [[],[]]
for i in range(0,48):
    Block1X[0].append(Block1[0][i])
    Block1X[1].append(Block1[1][i] + Block1[2][i] + Block1[3][i])
    
# Block 2 (structurally identical to blocks 4 and 6)
Block2 = [[],[],[],[]]
for i in range(0,48):
    Block2[0].append(i+1)
for i in range(0,16):
    Block2[1].append('N')
for i in range(16,32):
    Block2[1].append('S')
for i in range(32,48):
    Block2[1].append('D')
for i in range(0,4):
    Block2[2].append('SN')
for i in range(4,8):
    Block2[2].append('NS')
for i in range(8,12):
    Block2[2].append('TD')
for i in range(12,16):
    Block2[2].append('ND')
for i in range(16,48):
    Block2[2].append(Block2[2][i-16])
Block2[3].append('SN')
Block2[3].append('NS')
Block2[3].append('TD')
Block2[3].append('ND')
for i in range(4,48):
    Block2[3].append(Block2[3][i-4])
Block2X = [[],[]]
for i in range(0,48):
    Block2X[0].append(Block2[0][i])
    Block2X[1].append(Block2[1][i] + Block2[2][i] + Block2[3][i])

print('Done')
print('')
    
howMany = 0
while howMany == 0:
    howMany = input('How many PPR proteins do you want to make? ')
    if is_number(howMany) == True:
        if float(howMany) <= 0:
            howMany = 0
            print('Must be a positive integer')
        else:
            if (float(howMany)).is_integer() == False:
                howMany = 0
                print('Must be a positive integer')
            else:
                howMany = int(howMany)
    else:
        howMany = 0
        print('Must be a positive integer')
print('')
    
finalPlasmids, plates, wells, PPRID = [], [], [], []
for k in range(0, howMany):
    
    # Binding sequence?
    bindingSequence = 0
    while bindingSequence == 0:
        if howMany == 1:
            bindingSequence = input('Please enter the binding sequence: ')
        else:
            bindingSequence = input('Please enter binding sequence number ' + str(k+1) + ': ')
        print('')
        if bindingSequence == '':
            bindingSequence = 0
            print('Invalid entry')
            print('')
        elif len(bindingSequence) != 10 and len(bindingSequence) != 15:
            bindingSequence = 0
            print('The length of the binding sequence must be either 10 or 15')
            print('')
        else:
            bindingSequenceX = bindingSequence
            bindingSequence = bindingSequence.upper()
            i = 0
            while i in range(0, len(bindingSequenceX)):
                if bindingSequence[i] not in {'T','U','C','A','G'}:
                    i = len(bindingSequence)
                    bindingSequence = 0
                    print('Invalid sequence')
                    print('')
                else:
                    i = i + 1
    if 'T' in bindingSequence:
        bindingSequence = bindingSequence.replace('T','U')
            
    # Split the binding sequence
    bindingSequenceX = ''
    for i in range(0, len(bindingSequence)):
        if i == len(bindingSequence) - 1:
            bindingSequenceX = bindingSequenceX + bindingSequence[i]
        else:
            bindingSequenceX = bindingSequenceX + bindingSequence[i] + '.'
    listBindingSequence = bindingSequenceX.split('.')
        
    # Build the PPR
    i, j, PPR = 0, 0, ''
    while i <= len(bindingSequence) - 1:
        if bindingSequence[i] == RNAbases[j]:
            if i == 0:
                PPR = PPRmotifs[j]
                i = i + 1
                j = 0
            else:
                PPR = PPR + ' ' + PPRmotifs[j]
                i = i + 1
                j = 0
        else:
            j = j + 1
    listPPR = PPR.split()
        
    # Divide listPPR in 2.5mers
    PPRStr = ''
    counter = 0
    for i in range(0,len(listPPR)):
        for j in range(0,2):
            if counter < 5:
                PPRStr = PPRStr + listPPR[i][j]
                counter = counter + 1
            else:
                PPRStr = PPRStr + '.' + listPPR[i][j]
                counter = 1
    PPRStr = PPRStr.split('.')
        
    # Plasmids
    plasmids = []
    if len(bindingSequence) == 5:
        plasmids.append(Block1X[0][Block1X[1].index(PPRStr[0])])
        plasmids.append(Block2X[0][Block2X[1].index(PPRStr[1])])
    elif len(bindingSequence) == 10:
        plasmids.append(Block1X[0][Block1X[1].index(PPRStr[0])])
        plasmids.append(Block2X[0][Block2X[1].index(PPRStr[1])])
        plasmids.append(Block1X[0][Block1X[1].index(PPRStr[2])])
        plasmids.append(Block2X[0][Block2X[1].index(PPRStr[3])])
    elif len(bindingSequence) == 15:
        plasmids.append(Block1X[0][Block1X[1].index(PPRStr[0])])
        plasmids.append(Block2X[0][Block2X[1].index(PPRStr[1])])
        plasmids.append(Block1X[0][Block1X[1].index(PPRStr[2])])
        plasmids.append(Block2X[0][Block2X[1].index(PPRStr[3])])
        plasmids.append(Block1X[0][Block1X[1].index(PPRStr[4])])
        plasmids.append(Block2X[0][Block2X[1].index(PPRStr[5])])
    plasmidsLen = len(plasmids)
    for i in range(0,plasmidsLen):
        finalPlasmids.append('R' + str(i+1) + '-' + str(plasmids[i]))
            
    # Plates
    plates += 2 * [1]
    if plasmidsLen == 4:
        plates += 2 * [2]
    if plasmidsLen == 6:
        plates += 2 * [2]
        plates += 2 * [3]
        
    # Wells
    for i in range(0,plasmidsLen):
        x, y = coordinates(allPlasmidsX,finalPlasmids[i+len(PPRID)])[0], coordinates(allPlasmidsX,finalPlasmids[i+len(PPRID)])[1]
        if x == 0 or x == 8 or x == 16:
            wells.append('A' + str(y+1))
        if x == 1 or x == 9 or x == 17:
            wells.append('B' + str(y+1))
        if x == 2 or x == 10 or x == 18:
            wells.append('C' + str(y+1))
        if x == 3 or x == 11 or x == 19:
            wells.append('D' + str(y+1))
        if x == 4 or x == 12 or x == 20:
            wells.append('E' + str(y+1))
        if x == 5 or x == 13 or x == 21:
            wells.append('F' + str(y+1))
        if x == 6 or x == 14 or x == 22:
            wells.append('G' + str(y+1))
        if x == 7 or x == 15 or x == 23:
            wells.append('H' + str(y+1))
        
    for i in range(0, plasmidsLen):
        PPRID.append(k+1)
        
# Display results
df = pd.DataFrame({
    'PPR':PPRID,
    'Plasmids':finalPlasmids,
    'Plate':plates,
    'Well':wells})
print(df.to_string(index = False))
print('')

# -----------------------
# OT-2 protocol generator
# -----------------------

# Extract experiment variables from dataframe
numPPR = df['PPR'].iloc[-1]

# Map PPR output wells
outwells = []
for i in ['E','F','G','H']:
    for j in range(12):
        outwells.append(i + str(j+1))

# Read head and commands files
with open('Head.py') as txt:
    head = txt.read()
with open('Commands.py') as txt:
    commands = txt.read()

# Concatenate head and command files to generate protocol
for i in range(numPPR):
    # Select subset of PPR info from df
    dfPPR = df[df['PPR'] == i+1]
    
    # Change flag10PPR to False if assembling a 15nt-targeting PPR
    if len(dfPPR.index) == 6:
        commands = commands.replace('True', 'False')
    
    # Add commands to head for each PPR    
    head += '\n'
    head += commands
            
    # Replace tags in commands with info on dfPPR
    for j in range(len(dfPPR.index)) :
        head = head.replace('WELL' + str(j+1), dfPPR['Well'].iloc[j])
        head = head.replace('OUT', outwells[i])
    
    # Generate output file with timestamp
    with open(dt.datetime.now().strftime('%Y-%m-%d_%H.%M.%S') + 
              '_OTprotocol.py', 'w') as txt:
        txt.write(head)
   
print('Protocol file has been created')
print('')
input('Press ENTER to exit')


# In[ ]:




