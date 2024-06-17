from Bio import PDB
import os
import subprocess
import csv

def extract_and_save_parts(pdb_file, output_dir):
    parser = PDB.PDBParser()
    structure = parser.get_structure("PDB_structure", pdb_file)
    io = PDB.PDBIO()
    
    receptor_file = os.path.join(output_dir, f"receptor_{os.path.basename(pdb_file)}")
    ligand_file = os.path.join(output_dir, f"ligand_{os.path.basename(pdb_file)}")
    
    # Extract and save receptor (ATOM records)
    io.set_structure(structure)
    io.save(receptor_file, select=ReceptorSelect())
    
    # Extract and save ligand (HETATM records)
    io.set_structure(structure)
    io.save(ligand_file, select=LigandSelect())
    
    return receptor_file, ligand_file

class ReceptorSelect(PDB.Select):
    def accept_residue(self, residue):
        return residue.id[0] == " "  # Select only standard residues (ATOM records)

class LigandSelect(PDB.Select):
    def accept_residue(self, residue):
        return residue.id[0].startswith("H_")  # Select only HETATM records, typically ligands

def calculate_cnn_affinity(receptor, ligand):
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
    lines = output.split('\n')
    for line in lines:
        if "CNNscore" in line:
            parts = line.split()
            score = parts[-1]
            return float(score)
    return None

def process_pdbs(folder_path, output_csv):
    with open(output_csv, 'w', newline='') as csvfile:
        fieldnames = ['Receptor', 'Ligand', 'CNN_Affinity']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        
        for filename in os.listdir(folder_path):
            if filename.endswith(".pdb"):
                pdb_path = os.path.join(folder_path, filename)
                receptor_path, ligand_path = extract_and_save_parts(pdb_path, folder_path)
                output = calculate_cnn_affinity(receptor_path, ligand_path)
                affinity = extract_affinity(output)
                if affinity is not None:
                    writer.writerow({'Receptor': receptor_path, 'Ligand': ligand_path, 'CNN_Affinity': affinity})
                    print(f"Processed {filename}: CNN Affinity = {affinity}")
                else:
                    print(f"Failed to extract affinity for {filename}")

# Example usage
folder_path = '/home/kang/bioinfo_study/VP2_PDB/'
output_csv = '/home/kang/bioinfo_study/VP2_PDB//affinity_results.csv'
process_pdbs(folder_path, output_csv)

