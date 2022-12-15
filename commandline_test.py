#!/usr/bin/env python
import os
import sys
import csv
import pandas as pd
import numpy as np
import sqlite3 as sq
from tqdm import tqdm
from docopt import docopt
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import AllChem, Descriptors, FilterCatalog



cmd_str = """Usage: 
final.py --classification <MOLECULE> --in <INFILE_NAME> --out <OUTFILE_NAME> --save_as <TYPE>

--classification type small_molecule | peptide | big_molecule
--in INFILE_NAME  input SMILES(default),Amino Acid sequence or Nuclic Acid Sequence in single column CSV file
--save_as csv(default) | excel | sqlite3 | 
--out OUTFILE_NAME  output csv 

Example : python commandline_test.py --classification peptide --in ../data/test_peptide.csv --out output --save_as sqlite3
python commandline_test.py --classification peptide --in ../data/test_peptide.csv --out output --save_as excel
python commandline_test.py --classification peptide --in ../data/test_peptide.csv --out output --save_as csv
python commandline_test.py --classification small_molecule --in ../data/temp.smi --out output --save_as sqlite3
python commandline_test.py --classification small_molecule --in ../data/temp.smi --out output --save_as excel
python commandline_test.py --classification small_molecule --in ../data/temp.smi --out output --save_as csv

"""

def lipinski_drug_like_ness(H_bond_acceptors,H_bond_doner,Molecular_Weight,LogP):
	if H_bond_acceptors > 10 or H_bond_doner > 5 or Molecular_Weight > 500 or LogP > 5 :
		return 1
	else:
		return 0

def ghose_drug_like_ness_prefer(Molecular_Weight,LogP,molecular_refractivity,number_of_atoms):
	if (1.3 > LogP or LogP > 4.1) and (230 < Molecular_Weight or  Molecular_Weight > 390) and (70 > molecular_refractivity or molecular_refractivity > 110) and (30 > number_of_atoms or number_of_atoms > 55):
		return 1
	else:
		return 0

def veber_drug_like_ness(tPSA,rotatable_bonds):
	if tPSA > 140 or rotatable_bonds > 10 :
		return 1
	else :
		return 0

def muegge_drug_like_ness(Molecular_Weight, H_bond_acceptors, H_bond_doner, tPSA, LogP, number_of_rings,number_of_hetro_atoms,number_of_carbon,rotatable_bonds):
	if 200 > Molecular_Weight < 600  or H_bond_acceptors > 10 or H_bond_doner > 5 or tPSA > 150 and -2 > LogP > 5 or number_of_rings > 7 or number_of_hetro_atoms < 2 or number_of_carbon < 5 or rotatable_bonds > 15:
		return 1
	else :
		return 0

def egan_drug_like_ness(tPSA,LogP):
	if tPSA > 131.6 or LogP > 5.88 :
		return 1
	else:
		return 0

params1 = FilterCatalog.FilterCatalogParams()
params1.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
catalog1 = FilterCatalog.FilterCatalog(params1)
def pains(mol2):
	if catalog1.HasMatch(mol2):
		return 1
	else:
		return 0

params = FilterCatalog.FilterCatalogParams()
params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.BRENK)
catalog = FilterCatalog.FilterCatalog(params)
def brenk(mol1):
	if catalog.HasMatch(mol1) :
		return 1
	else:
		return 0

def read_data(infile_name,classification="small_molecule"):
	df = pd.DataFrame()
	if not os.path.exists(infile_name):
		print(f"Error: {infile_name} not found")
		sys.exit(0)
	if classification == "peptide":
		all_molecule = np.array([])
		with open(infile_name,'r') as csvfile:
			csvreader = csv.reader(csvfile)
			for each_row in csvreader:
				mol = Chem.MolFromSequence(each_row[0])
				all_molecule = np.append(all_molecule,[(Chem.MolToSmiles(mol,isomericSmiles=False))])

		df['SMILES']=pd.Series(all_molecule)
	else :
		df["SMILES"] = pd.read_csv(infile_name,header=None)
	print(f"Read {df.shape[0]} records from {infile_name}", file=sys.stderr)
	return df

def physicochemical_property(data_file):
	df = pd.DataFrame()
	df["SMILES"] = data_file
	mols = [Chem.MolFromSmiles(x) for x in data_file]
	df["Formula"] = [rdMolDescriptors.CalcMolFormula(x) for x in mols]
	df["Molecular_Weight"] = [Descriptors.ExactMolWt(x) for x in mols]
	df["Number_of_Atoms"] = [x.GetNumAtoms() for x in mols]
	df["tPSA"] = [Descriptors.TPSA(x) for x in mols]
	df["LogP"] = [Descriptors.MolLogP(x) for x in mols]
	df["FCSP3"] = [rdMolDescriptors.CalcFractionCSP3(x) for x in mols]
	df["H_bond_doner"] = [Chem.Lipinski.NumHDonors(x) for x in mols]
	df["H_bond_acceptors"] = [Chem.Lipinski.NumHAcceptors(x) for x in mols]
	df["Aromatic_Atoms"] = [len(list(x.GetAromaticAtoms())) for x in mols]
	df["Molecular_Refractivity"] = [Chem.Crippen.MolMR(x) for x in mols]
	df["Number_of_Rings"] = [rdMolDescriptors.CalcNumRings(x) for x in mols]
	df["Number_of_hetro_atoms"] = [rdMolDescriptors.CalcNumHeteroatoms(x) for x in mols]
	df["Number_of_RotatableBonds"] = [Chem.Lipinski.NumRotatableBonds(x) for x in mols]
	return df

def medicinal_chemistry(data_file):
	df = pd.DataFrame(columns=["Lipinski_Violations","Ghose_Violations","Veber_Violations","Muegge_Violations","Egan_Violations","PAINS_Alerts","BRENK_Alerts"])
	i = 1
	for row in data_file.itertuples(index=True, name='Pandas'):
		Lipinski_Violations=lipinski_drug_like_ness(row.H_bond_acceptors,row.H_bond_doner,row.Molecular_Weight,row.LogP)
		Ghose_Violations = ghose_drug_like_ness_prefer(row.Molecular_Weight,row.LogP,row.Molecular_Refractivity,row.Number_of_Atoms)
		Veber_Violations = veber_drug_like_ness(row.tPSA,row.Number_of_RotatableBonds)
		carbon = Chem.rdqueries.AtomNumEqualsQueryAtom(6)
		Muegge_Violations = muegge_drug_like_ness(row.Molecular_Weight,row.H_bond_acceptors,row.H_bond_doner, row.tPSA, row.LogP, row.Number_of_Rings, row.Number_of_hetro_atoms, len(Chem.MolFromSmiles(row.SMILES).GetAtomsMatchingQuery(carbon)), row.Number_of_RotatableBonds )
		Egan_Violations = egan_drug_like_ness(row.tPSA,row.LogP)
		PAINS_Alerts = pains(Chem.MolFromSmiles(row.SMILES))
		BRENK_Alerts = brenk(Chem.MolFromSmiles(row.SMILES))
		df.loc[i] = [Lipinski_Violations,Ghose_Violations,Veber_Violations,Muegge_Violations,Egan_Violations,PAINS_Alerts,BRENK_Alerts ]
		i = i +1
	return df

def output_save(outfile_name, save_as,physicochemical_property_data_frame,medicinal_chemistry_data_frame):
	if save_as == "excel":
		writer = pd.ExcelWriter(outfile_name+'.xlsx', engine='xlsxwriter')
		physicochemical_property_data_frame.to_excel(writer, sheet_name='Sheet1')
		medicinal_chemistry_data_frame.to_excel(writer, sheet_name='Sheet2')
		writer.close()
	elif save_as == "sqlite3" :
		conn = sq.connect(outfile_name+'.sqlite'.format("physicochemical_property"))
		physicochemical_property_data_frame.to_sql("physicochemical_property", conn, if_exists='replace')
		conn = sq.connect(outfile_name+'.sqlite'.format("medicinal_property"))
		medicinal_chemistry_data_frame.to_sql("medicinal_property", conn, if_exists='replace')
		conn.close()
	else :
		physicochemical_property_data_frame.to_csv(outfile_name+"_physicochemical_property.csv")
		medicinal_chemistry_data_frame.to_csv(outfile_name+"_medicinal_property.csv")



def main(cmd_input_str):
	cmd_input = docopt(cmd_input_str)
	classification = cmd_input.get("--classification")
	infile_name = cmd_input.get("--in")
	outfile_name = cmd_input.get("--out")
	save_as = cmd_input.get("--save_as")
	datafile1 = read_data(infile_name,classification)
	physicochemical_property_data_frame = physicochemical_property(datafile1["SMILES"])
	medicinal_chemistry_data_frame = medicinal_chemistry(physicochemical_property_data_frame)
	output_save(outfile_name,save_as,physicochemical_property_data_frame,medicinal_chemistry_data_frame)

if __name__ == "__main__":
	main(cmd_str)