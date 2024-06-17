import os
from Bio import PDB
from openbabel import openbabel

def convert_pdb_to_mol(pdb_path, mol_path):
    """Convert PDB files to MOL format using Open Babel."""
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("pdb", "mol")
    mol = openbabel.OBMol()
    
    if obConversion.ReadFile(mol, pdb_path):
        obConversion.WriteFile(mol, mol_path)
    else:
        print(f"Failed to read {pdb_path}")

def extract_atom_records(pdb_file, output_file):
    """Extract ATOM and HETATM records from a PDB file."""
    parser = PDB.PDBParser()
    structure = parser.get_structure("Structure", pdb_file)
    
    with open(output_file, 'w') as out_f:
        io = PDB.PDBIO()
        io.set_structure(structure)
        io.save(out_f, select=AtomSelect())

class AtomSelect(PDB.Select):
    """Custom selection class to include ATOM and HETATM records."""
    def accept_residue(self, residue):
        return residue.id[0] != ' '  # Accept residues that are ATOMs or HETATMs

def process_pdb_files(folder_path):
    for filename in os.listdir(folder_path):
        if filename.endswith(".pdb"):
            pdb_path = os.path.join(folder_path, filename)
            temp_pdb_path = os.path.join(folder_path, "temp.pdb")  # Temporary file for filtered ATOM/HETATM
            mol_path = os.path.join(folder_path, filename.replace(".pdb", ".mol"))
            
            # Extract ATOM and HETATM records and save to a temporary PDB
            extract_atom_records(pdb_path, temp_pdb_path)
            
            # Convert the temporary PDB to MOL format
            convert_pdb_to_mol(temp_pdb_path, mol_path)
            
            # Remove the temporary file
            os.remove(temp_pdb_path)
            print(f"Converted {pdb_path} to {mol_path}")

# Example usage
folder_path = '/home/kang/bioinfo_study/VP2_PDB'
process_pdb_files(folder_path)

