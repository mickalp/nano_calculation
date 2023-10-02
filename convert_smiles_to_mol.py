#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  2 15:36:34 2023

@author: michal
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

# =============================================================================
# generate 3d structures
# =============================================================================


#read excel file

sh = 'smiles'

df = pd.read_excel('smi_and_whim.xlsx', sheet_name = sh)


#number of compound from reduced matrix of compounds
number_of_compound = df[df.columns[0]]

#list of reduced list of smiles of ionic liquids
list_of_smi = df[df.columns[1]]



# mol = Chem.MolFromSmiles(list_of_smi[0])

list_of_mol = []
list_of_failed = []
w = 0


for i in list_of_smi:
    try:
        mol = Chem.MolFromSmiles(i)
        mol_with_H = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol_with_H)
        AllChem.MMFFOptimizeMolecule(mol_with_H)
# =============================================================================
#     be careful here to not save it at previous directory after optimization
# =============================================================================
        # Draw.MolToFile(mol,f'{i}.png')
        # fout = Chem.SDWriter(f'./{i}.mol')
        # fout.write(mol_with_H)
        # fout.close()
        
        
    except Exception:
        wdk = (f'BADBADBAD    {i}')
        list_of_failed.append(wdk)
        continue



mol = Chem.MolFromSmiles(list_of_smi[115])
# mol_with_H = Chem.AddHs(mol)
mol_with_H = mol
AllChem.EmbedMolecule(mol_with_H)
AllChem.UFFOptimizeMolecule(mol_with_H, numThreads=0)