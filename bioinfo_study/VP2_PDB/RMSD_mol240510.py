import os
from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd

def calculate_rmsd(file1, file2):
    mol1 = Chem.MolFromMolFile(file1)
    mol2 = Chem.MolFromMolFile(file2)
    
    if mol1 is None or mol2 is None:
        print(f"Failed to load molecules from {file1} or {file2}")
        return None
    
    # Align the molecules using RDKit's built-in function
    try:
        mol1 = Chem.AddHs(mol1)
        mol2 = Chem.AddHs(mol2)
        AllChem.EmbedMolecule(mol1, randomSeed=0xf00d)
        AllChem.EmbedMolecule(mol2, randomSeed=0xf00d)
        rmsd = AllChem.GetBestRMS(mol1, mol2)
        return rmsd
    except Exception as e:
        print(f"Error in RMSD calculation for {file1} and {file2}: {e}")
        return None

def compare_all_mols(folder_path):
    mol_files = [os.path.join(folder_path, f) for f in os.listdir(folder_path) if f.endswith('.mol')]
    results = []
    for i in range(len(mol_files)):
        for j in range(i + 1, len(mol_files)):
            rmsd_value = calculate_rmsd(mol_files[i], mol_files[j])
            if rmsd_value is not None:
                results.append({
                    'File1': os.path.basename(mol_files[i]),
                    'File2': os.path.basename(mol_files[j]),
                    'RMSD': rmsd_value
                })
                print(f"RMSD between {mol_files[i]} and {mol_files[j]}: {rmsd_value}")
    return results

# Example usage
folder_path = '/home/kang/bioinfo_study/VP2_PDB'
results = compare_all_mols(folder_path)
results_df = pd.DataFrame(results)
print(results_df)

