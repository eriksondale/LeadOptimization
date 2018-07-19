

#Note this is just a simple linear approach to get preliminary results
# it can easily be expanded upon to get more quality leads

# Things to improve:
#   - Detecting functional groups
#   - Settling equal situations when determining scaffold and additions


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

# Reaction checker
def validRxn(reactant, reaction, revRxn):
	#print(reactant)
	#print(reaction)
	try:
		product = None
		product = reaction.RunReactants(reactant)
		print("Reaction was tested...")
		if(len(product) == 0):
			print('product is none')
			return None
		else:
			#reverseReactionPlan = reactionPlan[reactionPlan.find('>>')+2:] + '>>' + reactionPlan[0:reactionPlan.find('>>')]
			#reverseReaction = AllChem.ReactionFromSmarts(reverseReactionPlan)
			for prod in product:
				try:
					reactant = revRxn.RunReactants(prod)
					for pairs in reactant:
						for molecule in pairs:
							moleculeBit = FingerprintMols.FingerprintMol(molecule)
							compoundBit= FingerprintMols.FingerprintMol()
							similarity = DataStructs.FingerprintSimilarity(moleculeBit, compoundBit)
							if(similarity == 1): # Rxn is valid b/c product and reverse product is found to be same
								return prod
				except:
					pass
		  	print("reactant not reformed")
			return None
	except Exception as e:
		print("Overall error: "+ str(e))
		return None

# Lead
with open(arg[1],"r") as leadFile:
    lead = Chem.MolFromSmiles(leadFile.read())
    #AllChem.SanitizeMol(lead)
    leadFile.close()

# Rxn
with open(arg[2],"r") as rxnFile:
    rxnText = rxnFile.read()
    rxn = Chem.AllChem.ReactionFromSmarts(rxnText)
    rxnReversed = Chem.AllChem.ReactionFromSmarts(rxnText[rxnText.find('>>')+2:] + '>>' + rxnText[0:rxnText.find('>>')])
    rxnFile.close()

reactants = []
reactants.append(lead)
products = validRxn(reactants, rxn)

if products is None:
	print("Reaction not applicable to lead...")
	sys.exit()

# Scaffold
scaffold = None
for prods in products:
        AllChem.SanitizeMol(prods)
        if(scaffold is None):
            scaffold = prods
        else:
            if(Descriptors.MolWt(prods, 0) > Descriptors.MolWt(scaffold, 0)):
                scaffold = prods

print("Scaffold Established: " + Chem.MolToSmiles(scaffold))
# Scaffold established
fragList = []

with open(arg[3],"r") as smallMolFile:
    for line in smallMolFile:
        tempMol = Chem.MolFromSmiles(line)
	reactants = []
	reactants.append(tempMol)
        tempProducts = validRxn(reactants, rxn)
        for products in tempProducts:
                fragList.append(products)

# Saving Fragments for Later
with open("fragments.smi","a") as fragFile:
    for frag in fragList:
        try:
            Chem.SanitizeMol(frag)
            fragFile.write(Chem.MolToSmiles(frag) + "\n")
        except:
            pass
    fragFile.close()

print("Number of Fragments: " +  str(len(fragList)))

# Rebuilding Molecule
print("Optimizing Lead....")

with open("optimizedLeads.smi","a") as leadFile:
    for frag in fragList:
        try:
            newLead = rxnReversed.RunReactants([scaffold, frag],rxnReversed)
            for leads in newLead:
                    try:
                        Chem.SanitizeMol(leads)
                        leadFile.write(Chem.MolToSmiles(leads) + "\n")
                    except:
                        pass
        except:
            pass


# Outputing Optimized Leads
print("New Leads Ouput to Optimized Lead File")
