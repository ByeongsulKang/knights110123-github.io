import os
import subprocess
from rdkit import Chem
from rdkit.Chem import AllChem


# 입력 및 출력 폴더 설정
input_folder = "/home/kang/TNF-a/Newtrail/ligands_pdb"
reconstructed_folder = "/home/kang/TNF-a/Newtrail/reconstructed_pdb"
output_folder = "/home/kang/TNF-a/Newtrail/ligands_pdbqt"
os.makedirs(reconstructed_folder, exist_ok=True)
os.makedirs(output_folder, exist_ok=True)

# PDB 파일의 CONECT 레코드를 제거하는 함수
def remove_conect_records(input_pdb, output_pdb):
    with open(input_pdb, 'r') as infile, open(output_pdb, 'w') as outfile:
        for line in infile:
            if not line.startswith('CONECT'):
                outfile.write(line)

# PDB 파일을 재구성하여 저장하는 함수
def reconstruct_pdb(input_pdb, output_pdb):
    print(f"Reading PDB file: {input_pdb}")
    mol = Chem.MolFromPDBFile(input_pdb, removeHs=False)
    if mol is None:
        raise RuntimeError(f"Error: Could not read molecule from {input_pdb}")
    print(f"Successfully read PDB file: {input_pdb}")

    with open(output_pdb, 'w') as f:
        f.write(Chem.MolToPDBBlock(mol))
    print(f"Reconstructed PDB file saved to {output_pdb}")

# OpenBabel을 사용하여 PDB 파일을 PDBQT 파일로 변환하는 함수
def convert_pdb_to_pdbqt(input_pdb, output_pdbqt):
    command = ['obabel', '-ipdb', input_pdb, '-opdbqt', '-O', output_pdbqt]
    print(f"Running command: {' '.join(command)}")
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error converting {input_pdb} to PDBQT: {result.stderr}")
        raise RuntimeError(f"Error converting {input_pdb} to PDBQT: {result.stderr}")
    if not os.path.isfile(output_pdbqt):
        print(f"Command executed: {' '.join(command)}")
        print(f"Standard output: {result.stdout}")
        print(f"Standard error: {result.stderr}")
        raise RuntimeError(f"Error: {output_pdbqt} was not created.")
    print(f"Successfully converted {input_pdb} to {output_pdbqt}")

# 입력 폴더에서 PDB 파일을 찾아 CONECT 레코드를 제거하고 재구성한 후 PDBQT 파일로 변환
for file_name in os.listdir(input_folder):
    if file_name.endswith('.pdb'):
        input_pdb_path = os.path.join(input_folder, file_name)
        temp_pdb_path = os.path.join(reconstructed_folder, f"no_conect_{file_name}")
        reconstructed_pdb_path = os.path.join(reconstructed_folder, file_name)
        output_pdbqt_path = os.path.join(output_folder, f"{os.path.splitext(file_name)[0]}.pdbqt")
        
        try:
            # CONECT 레코드 제거
            remove_conect_records(input_pdb_path, temp_pdb_path)
            # PDB 파일 재구성
            reconstruct_pdb(temp_pdb_path, reconstructed_pdb_path)
            # 재구성된 PDB 파일을 PDBQT로 변환
            convert_pdb_to_pdbqt(reconstructed_pdb_path, output_pdbqt_path)
        except Exception as e:
            print(e)

