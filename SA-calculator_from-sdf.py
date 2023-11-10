#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  8 09:45:50 2023

@author: Ehsan
"""
import os
from rdkit import Chem
from rdkit.Chem import PandasTools
import pandas as pd
import argparse
# from openbabel.pybel import *
import time
import glob
from rdkit.Chem import RDConfig
from rdkit.Chem import Draw
import sys
sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))
import sascorer
import matplotlib.pyplot as plt
import pandas as pd
from rdkit.Chem import AllChem


directory = ' '
path = '/media/arma/DATA/S-Ehsan/SARCOPENIA/3D-MCTS/record/'
extension = 'sdf'
os.chdir(path)
result = glob.glob('*.{}'.format(extension))

obabel='/usr/local/bin/obabel'
parser = argparse.ArgumentParser(description='SA calculation by converting sdf to smi')
parser.add_argument('--obabel', action="store", type=str, help='the path for obabel',
                    default='/usr/local/bin/obabel')

args = parser.parse_args()

os.system('mkdir SA-pass SmileFiles')

smi = ''
smilist=[]
smilescore=[]
filename = []
threshold = 5

# Function to convert SDF to SMILES
def convert_sdf_to_smiles(result):

    for i in result:
        s = i.replace('.sdf', '')
        os.system(rf'{obabel} {i} -O {s}.smi')
    os.system(rf'mv *.smi SmileFiles/')
    os.chdir(path+'SmileFiles')
    sresult = glob.glob('*.{}'.format('smi'))
    for z in sresult:
        read=pd.read_csv(z,header=None)
        smilist.append(read[0][0])
        m = Chem.MolFromSmiles(read[0][0])
        name=z.replace('.smi', '')
        filename.append(name)
        sa_score= "%.2f"%sascorer.calculateScore(m)
        smilescore.append(sa_score)
    smi= pd.DataFrame(filename, columns=['name']) 
    smi['Smiles']= smilist
    smi['SA-score']= smilescore
    smi['SA-score'] = pd.to_numeric(smi['SA-score'], errors='coerce', downcast='float')
    # Draw.MolsToGridImage([m], legends=["%.2f"%sascorer.calculateScore(m)])
    sa_pass = pd.DataFrame(smi[smi['SA-score'] <= threshold])
    sa_pass.to_csv('sa_pass.csv', index=False)
    return sa_pass

def get_pass_files(sa_pass):
    for l in sa_pass.iloc[:,0]:
        os.chdir(path)
        file=l+'.sdf'
        os.system(rf'mv {file} SA-pass/')
    return print(f'Your small molecules with SA-score below {threshold} move to SA-pass directory')
    

print("Conversion complete. SMILES data saved to SmileFiles Directory...\n")
print('for all Smiles SA score calculated ...\n')


if __name__ == "__main__":
    sa_pass=convert_sdf_to_smiles(result)
    get_pass_files(sa_pass)
