import os

def extract_records(pdb_file_path, output_file_path, record_type):
    with open(pdb_file_path, 'r') as infile, open(output_file_path, 'w') as outfile:
        for line in infile:
            if line.startswith(record_type):
                outfile.write(line)

# TACE.pdb 파일에서 ATOM 열만 추출하여 receptor.pdb로 저장
pdb_file_path = "/home/kang/TNF-a/TACE.pdb"
receptor_output_file_path = "/home/kang/TNF-a/TACE_receptor.pdb"
extract_records(pdb_file_path, receptor_output_file_path, "ATOM")

# ligands_3D 디렉토리 내의 모든 PDB 파일에서 HETATM 열만 추출하여 각각 저장
ligands_dir = "/home/kang/TNF-a/ligands_3D"
ligand_files = [f for f in os.listdir(ligands_dir) if f.endswith('.pdb')]

for ligand_file in ligand_files:
    ligand_file_path = os.path.join(ligands_dir, ligand_file)
    ligand_output_file_path = os.path.join(ligands_dir, f"HETATM_{ligand_file}")
    extract_records(ligand_file_path, ligand_output_file_path, "HETATM")

print("Extraction complete.")

