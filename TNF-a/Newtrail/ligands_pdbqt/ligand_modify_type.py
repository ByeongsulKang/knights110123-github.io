import os
from rdkit import Chem
from rdkit.Chem import AllChem

# 입력 폴더 및 출력 폴더 설정
input_folder = "/home/kang/TNF-a/Newtrail/ligands_pdb"
reconstructed_folder = "/home/kang/TNF-a/Newtrail/reconstructed_pdb"
os.makedirs(reconstructed_folder, exist_ok=True)

# PDB 파일을 재구성하여 저장하는 함수
def reconstruct_pdb(input_pdb, output_pdb):
    mol = Chem.MolFromPDBFile(input_pdb, removeHs=False)
    if mol is None:
        raise RuntimeError(f"Error: Could not read molecule from {input_pdb}")
    with open(output_pdb, 'w') as f:
        f.write(Chem.MolToPDBBlock(mol))
    print(f"Reconstructed PDB file saved to {output_pdb}")

# 입력 폴더에서 PDB 파일을 찾아 재구성하고 출력 폴더에 저장
for file_name in os.listdir(input_folder):
    if file_name.endswith('.pdb'):
        input_pdb_path = os.path.join(input_folder, file_name)
        output_pdb_path = os.path.join(reconstructed_folder, file_name)
        
        try:
            # PDB 파일 재구성
            reconstruct_pdb(input_pdb_path, output_pdb_path)
        except Exception as e:
            print(e)

