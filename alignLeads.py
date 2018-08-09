import rdkit
from sys import argv as arg
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw,PyMol,MCS
#from rdkit.Chem.Draw import IPythonConsole
from rdkit import rdBase

if(len(arg)!=3):
    print("Script requires 1 .sdf file and 1 .smi file: one aligned and one to be aligned")
else:
    leadSuppl = Chem.SDMolSupplier(arg[1])
    for mol in leadSuppl:
        lead = mol
        Chem.SanitizeMol(lead)
    optSuppl = Chem.SmilesMolSupplier(arg[2])
    w = Chem.SDWriter('opleadswcoords.sdf')

    mols = [lead, optSuppl[0]]
    res=MCS.FindMCS(mols,threshold=0.9,completeRingsOnly=True)
    p = Chem.MolFromSmarts(res.smarts)
    core = AllChem.DeleteSubstructs(AllChem.ReplaceSidechains(Chem.RemoveHs(lead),p),Chem.MolFromSmiles('*'))
    core.UpdatePropertyCache()
    for mol in optSuppl:
        try:
            AllChem.ConstrainedEmbed(mol, core)
            w.write(mol)
        except:
            pass
