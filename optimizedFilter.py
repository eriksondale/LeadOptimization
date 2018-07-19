from __future__ import print_function
# import rdkit components
from rdkit import rdBase
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Descriptors
from rdkit.Chem import rdchem
import sys
from sys import argv as arg
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit import DataStructs

leadFile = open(arg[1],"r")
for line in leadFile:
    tempMol = Chem.MolFromSmiles(line.strip("\n"))
    if(Descriptors.MolWt(tempMol, 0) <= 500):
        print(Chem.MolToSmiles(tempMol))
