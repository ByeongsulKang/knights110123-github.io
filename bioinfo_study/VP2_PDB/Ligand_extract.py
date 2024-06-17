from Bio import PDB
import os

def extract_and_save_parts(pdb_file, output_dir):
    """Extract and save receptor and ligand parts from a PDB file."""
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

class ReceptorSelect(PDB.Select):
    def accept_residue(self, residue):
        return residue.id[0] == " "  # Select only standard residues (ATOM records)

class LigandSelect(PDB.Select):
    def accept_residue(self, residue):
        return residue.id[0].startswith("H_")  # Select only HETATM records, typically ligands

def process_pdbs(folder_path, output_folder):
    for filename in os.listdir(folder_path):
        if filename.endswith(".pdb"):
            pdb_path = os.path.join(folder_path, filename)
            extract_and_save_parts(pdb_path, output_folder)
            print(f"Processed {filename}")

# Example usage
folder_path = '/home/kang/bioinfo_study/VP2_PDB'
output_folder = '/home/kang/bioinfo_study/VP2_PDB/Receptor_Ligand.pdb'
process_pdbs(folder_path, output_folder)

