#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 20 10:55:01 2023

@author: admin
"""

import os
import pandas as pd
import numpy as np



#folder path to check token
dir_path = '/Users/admin/OneDrive - University of Gdansk/OneDrive - University of Gdansk (for Students)/Jakub Rudzinski/out_v1/out_files/' #big computer
# dir_path = '/Users/michal/OneDrive - University of Gdansk/OneDrive - University of Gdansk (for Students)/Jakub Rudzinski/out_v1/out_files/' #small computer

#list to store names of the files
res = []
os.chdir(dir_path)

for path in os.listdir(dir_path): #files in my path
    if os.path.isfile(os.path.join(dir_path, path)):
        res.append(path)

without_ext_str = []



res.sort()

    

homo_energy_list = []
lumo_energy_list = []


start_phrase = "                        Standard orientation:"

cuyt_phr = '---------------------------------------------------------------------'
str_end_ind = ' ---------------------------------------------------------------------'

band_gap_list = []
i=0
hartree_to_ev = 27.2114

not_work = []
ind_to_del = []




for _ in res:
    with open(_, "r") as input_file:
        file_contents = input_file.read()
        
        print(i)
        # print(res[i])
        i+=1
        
        
    
    # Split the file contents into sections separated by the end phrase
    sections = file_contents.split('Standard orientation:')
    
    # Iterate through the sections and find the last one that starts with the start phrase
    for section in reversed(sections):
        if section.strip().startswith(cuyt_phr):
            # Cut the text from the start phrase to the end of the section
            cut_text = section.strip()[len(start_phrase):]
            break
        
        
    cut_txt = cut_text.split('\n')
    
    index_HOMO = [index for (index, item) in enumerate(cut_txt) if 'Alpha  occ. eigenvalues' in item]
    index_LUMO  = [index for (index, item) in enumerate(cut_txt) if 'Alpha virt. eigenvalues' in item]
    
    try:
        HOMO_en = float(cut_txt[index_HOMO[-1]][-8:])*hartree_to_ev
        LUMO_en = float(cut_txt[index_LUMO[0]][30:40])*hartree_to_ev
        homo_energy_list.append(HOMO_en)
        lumo_energy_list.append(LUMO_en)

        
        band_gap = abs(HOMO_en-LUMO_en)
        band_gap_list.append(band_gap)

    except Exception:
        not_work.append(res[i-1])
        print('NOT FINAL COORDS!!!!!!!!!!!!!:', res[i-1])
        ind_to_del.append(i-1)
        

        continue


   
for item in not_work: res.remove(item)

for i in res:
    a = i[:-8]
    without_ext_str.append(a)
    
    
#%%

# =============================================================================
# DIPOLE MOMENT CUTTING    
# =============================================================================

import os
import pandas as pd
import numpy as np



#folder path
# dir_path = '/Users/admin/OneDrive - University of Gdansk/OneDrive - University of Gdansk (for Students)/Jakub Rudzinski/out_v1/out_files/' #big computer
dir_path = '/Users/michal/OneDrive - University of Gdansk/OneDrive - University of Gdansk (for Students)/Jakub Rudzinski/out_v1/out_files/' #small computer

#list to store names of the files
res = []
os.chdir(dir_path)

for path in os.listdir(dir_path): #files in my path
    if os.path.isfile(os.path.join(dir_path, path)):
        res.append(path)

without_ext_str = []



res.sort()

    

homo_energy_list = []
lumo_energy_list = []


start_phrase = " Charge=              0.0000 electrons\n"
cuyt_phr_2 = " Charge=              -0.0000 electrons\n"

cuyt_phr = ' Dipole moment (field-independent basis, Debye):'

str_end_ind = ' Dipole moment (field-independent basis, Debye):'
cheking_ph = 'Stationary point'

band_gap_list = []
i=0


not_work = []
ind_to_del = []



list_md = []

for _ in res:
    with open(_, "r") as input_file:
        file_contents = input_file.read()
        
        print(i)
        # print(res[i])
        i+=1
        
    if cheking_ph in file_contents:   
    
        # Split the file contents into sections separated by the end phrase
        sections = file_contents.split(' Charge=')
        
    

    
        # Iterate through the sections and find the last one that starts with the start phrase
        for section in reversed(sections):
            if cuyt_phr in section:
                
                # Cut the text from the start phrase to the end of the section


                dipole_mom = section[176:183]
                print(dipole_mom)
                try:
                    mm = float(dipole_mom)
                    list_md.append(mm)
                    
                except Exception:
                    print('pa')
                    continue
                break
    else:
        not_work.append(res[i])




for item in not_work: res.remove(item)

for i in res:
    a = i[:-8]
    without_ext_str.append(a)
    
    
#%%
    
print('IN NEXT CELL WILL BE SAVING OF FILES')    
    
#%%

to_sv = [without_ext_str, list_md]

data_save = {'smiles':to_sv[0], 'dipole moment [Debye]':to_sv[1]}
to_sv_df = pd.DataFrame(data_save)

to_sv_df.to_excel('/Users/michal/OneDrive - University of Gdansk/OneDrive - University of Gdansk (for Students)/Jakub Rudzinski/dane_do_modelowania/DM.xlsx')
# to_sv_df.to_csv('/Users/michal/OneDrive - University of Gdansk/OneDrive - University of Gdansk (for Students)/Jakub Rudzinski/dane_do_modelowania/homo.csv')




#%%

import pandas as pd
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
        # df.to_excel(name[:-4]+'/Users/michal/Library/CloudStorage/OneDrive-UniversityofGdansk/OneDrive - University of Gdansk (for Students)/Jakub Rudzinski/out_v1/not_finshed/_RDKit_2D.xlsx', index=False)
        # df.to_csv(name[:-4]+'/Users/michal/Library/CloudStorage/OneDrive-UniversityofGdansk/OneDrive - University of Gdansk (for Students)/Jakub Rudzinski/out_v1/not_finshed/_RDKit_2D.csv', index=False)
        df.to_excel(name[:-4]+'/Users/michal/OneDrive - University of Gdansk/OneDrive - University of Gdansk (for Students)/Jakub Rudzinski/dane_do_modelowania/_RDKit_2D_DM.xlsx', index=False)
        # df.to_csv(name[:-4]+'/Users/admin/OneDrive - University of Gdansk/OneDrive - University of Gdansk (for Students)/Jakub Rudzinski/dane_do_modelowania/_RDKit_2D.csv', index=False)
        # return df
        
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
#%%

for i in mols: 
    cl = calc(i).fill_missing()[0:]
        
    emp.append((cl))

#%%


import pandas as pd
import numpy as np

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
# df.to_csv('/Users/admin/OneDrive - University of Gdansk/OneDrive - University of Gdansk (for Students)/Jakub Rudzinski/dane_do_modelowania/mordred.csv')



#%%

replace_nan_with_mean(df)
# df.dropna(axis=1, inplace=True)

    






