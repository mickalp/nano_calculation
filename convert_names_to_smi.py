#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  5 10:09:40 2022

@author: admin
"""

import pandas as pd
import requests
import time
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import PyMol
from rdkit.Chem.Draw.MolDrawing import MolDrawing, DrawingOptions




df = pd.read_excel('BAZA DANYCH.xlsx', sheet_name='WHIM')


ids = df['IONIC LIQIUD']
identifires =  ids.to_numpy()


def CIRconvert(comp_name):
    try:
        url = 'https://cactus.nci.nih.gov/chemical/structure/' + comp_name + '/smiles'
        smiles = requests.get(url)
        smiles = smiles.text
        return smiles
    except:
        return 'Did not work'
    
smiles_list = []

start_time = time.time()

print('start time:\t', start_time)


for i in identifires:
    print(i)
    
    while True:
        smiles_list.append(CIRconvert(i))
        current_time = time.time()
        elapsed_time = current_time-start_time
        
        if elapsed_time > 6:
            print('Next iteration, time too long, not find SMILES for \t', i)
            del elapsed_time
            break
        
        

from collections import Counter

wdk = Counter(smiles_list)

sum('.' in s for s in smiles_list)

add_sm = pd.Series(smiles_list, name='SMI')

df2 = pd.concat([add_sm, df], ignore_index=True, axis=1)

df2 = df2[df2[0].str.contains("Page") == False]

df2.to_excel('smi_and_whim.xlsx')

#%%
