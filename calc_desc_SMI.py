#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 11 21:02:35 2023

@author: michal
"""


import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors


class RDKit_2D:
    def __init__(self, smiles):
        self.mols = [Chem.MolFromSmiles(i) for i in smiles]
        self.smiles = smiles
        
    def compute_2Drdkit(self, name):
        rdkit_2d_desc = []
        calc = MoleculeDescriptors.MolecularDescriptorCalculator([x[0] for x in Descriptors._descList])
        header = calc.GetDescriptorNames()
        
        for i in range(len(self.mols)):
            ds = calc.CalcDescriptors(self.mols[i])
            rdkit_2d_desc.append(ds)
            
        df = pd.DataFrame(rdkit_2d_desc,columns=header)
        df.insert(loc=0, column='smiles', value=self.smiles)
        df.dropna(axis=1, inplace=True)
        df.to_excel(name[:-4]+'/Users/michal/OneDrive - University of Gdansk/OneDrive - University of Gdansk (for Students)/Jakub Rudzinski/dane_do_modelowania/_RDKit_2D_DM.xlsx', index=False)

        
new = RDKit_2D(without_ext_str)

new.compute_2Drdkit('')


#%%

from mordred import descriptors, Calculator

calc = Calculator(descriptors, ignore_3D=False)
mol = Chem.MolFromSmiles('c1ccccc1')

to_c = calc(mol).fill_missing().drop_missing().name

dict_col = list(to_c._name_to_value.keys())



mols = [Chem.MolFromSmiles(smi) for smi in without_ext_str]
emp = []


for i in mols: 
    cl = calc(i).fill_missing()[0:]
        
    emp.append((cl))


def replace_nan_with_mean(df):
    for col in df.columns:
        num_nan = df[col].isna().sum()
        if num_nan == 1 or num_nan == 2:
            mean_value = df[col].mean()
            df[col].fillna(value=mean_value, inplace=True)
        elif num_nan > 2:
            df.drop(col, axis=1, inplace=True)
    return df




df = pd.DataFrame(emp, columns=dict_col)
df.dropna(axis=1, inplace=True)
df.to_excel('/Users/michal/OneDrive - University of Gdansk/OneDrive - University of Gdansk (for Students)/Jakub Rudzinski/dane_do_modelowania/mordred_DM.xlsx')


    






