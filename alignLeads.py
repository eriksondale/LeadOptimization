import rdkit
import sys
from sys import argv as arg
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw,PyMol,MCS
from rdkit import rdBase

print('NOTE: This will align substructures but remove any difference conformers may have: will fix in future!')
if(len(arg)!=3):
    print("Script requires 1 .sdf file and 1 .smi file: one aligned and one to be aligned")
else:
    leadSuppl = Chem.SDMolSupplier(arg[1])
    for mol in leadSuppl:
        lead = mol
        Chem.SanitizeMol(lead)
    if(arg[2].endswith('.smi')):
        optSuppl = Chem.SmilesMolSupplier(arg[2])
    elif(arg[2].endswith('sdf')):
        optSuppl = Chem.SDMolSupplier(arg[2])
    else:
        print('File type not supported')
        sys.exit()

    w = Chem.SDWriter('alignLeadsOut.sdf')

    mols = [lead, optSuppl[0]]
    res=MCS.FindMCS(mols,threshold=0.9,completeRingsOnly=True) # Calculates most common substructure and outputs SMARTS pattern
    p = Chem.MolFromSmarts(res.smarts) # Creates mol object of most common substructure
    core = AllChem.DeleteSubstructs(AllChem.ReplaceSidechains(Chem.RemoveHs(lead),p),Chem.MolFromSmiles('*'))
    core.UpdatePropertyCache()
    for mol in optSuppl:
        try:
            AllChem.ConstrainedEmbed(mol, core)
            w.write(mol)
        except:
            pass
