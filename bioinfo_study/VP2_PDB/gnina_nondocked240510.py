import subprocess
import csv
import os
from Bio import PDB

def extract_ligand(pdb_file, output_dir):
    """Extract ligands from PDB files and save as separate files."""
    parser = PDB.PDBParser()
    structure = parser.get_structure('PDB_structure', pdb_file)
    io = PDB.PDBIO()
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0][0] == 'H':  # HETATM records
                    io.set_structure(structure)
                    ligand_path = os.path.join(output_dir, f"ligand_{os.path.basename(pdb_file).replace('.pdb', '.mol')}")
                    io.save(ligand_path, select=LigandSelect())
                    return ligand_path
    return None

class LigandSelect(PDB.Select):
    def accept_residue(self, residue):
        return residue.id[0][0] == 'H'

def calculate_cnn_affinity(receptor, ligand):
    """Calculate CNN affinity using gnina."""
    cmd = [
        '/home/kang/bioinfo_study/gnina',
        '-r', receptor,
        '-l', ligand,
        '--score_only',
        '--cnn_model', 'general_default2018'
    ]
    result = subprocess.run(cmd, stdout=subprocess.PIPE, text=True)
    return result.stdout

def extract_affinity(output):
    """ Extract the CNN affinity score from gnina output. """
    lines = output.split('\n')
    for line in lines:
        if "CNNscore" in line:
            parts = line.split()
            score = parts[-1]
            return float(score)
    return None

def process_pdbs(folder_path, output_csv):
    """Process all PDB files in a folder to calculate CNN affinity."""
    receptors = [os.path.join(folder_path, f) for f in os.listdir(folder_path) if f.endswith('.pdb')]
    results = []

    with open(output_csv, 'w', newline='') as csvfile:
        fieldnames = ['Receptor', 'Ligand', 'CNN_Affinity']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for receptor in receptors:
            ligand = extract_ligand(receptor, folder_path)
            if ligand:
                output = calculate_cnn_affinity(receptor, ligand)
                affinity = extract_affinity(output)
                if affinity is not None:
                    writer.writerow({'Receptor': receptor, 'Ligand': ligand, 'CNN_Affinity': affinity})
                    results.append((receptor, ligand, affinity))
                    print(f"Processed {receptor}: CNN Affinity = {affinity}")
                else:
                    print(f"Failed to extract affinity for {receptor}")
            else:
                print(f"No ligand found in {receptor}")

    return results

# Example usage
folder_path = '/home/kang/bioinfo_study/VP2_PDB'
output_csv = '/home/kang/bioinfo_study/VP2_PDB/affinity_results.csv'
process_pdbs(folder_path, output_csv)

