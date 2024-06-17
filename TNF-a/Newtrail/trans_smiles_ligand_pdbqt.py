import os
from rdkit import Chem
from rdkit.Chem import AllChem
import subprocess

# 리간드 SMILES 목록
ligands = {
    "GW3333": "CCC(C)C(C(=O)NC1=CC=CC=N1)NC(=O)C(CC(C)C)C(C(C)C)N(C=O)O",
    "TMI-1": "CC#CCOC1=CC=C(C=C1)S(=O)(=O)N2CCSC([C@@H]2C(=O)NO)(C)C",
    "BMS-561392": "CC1=NC2=CC=CC=C2C(=C1)COC3=CC=C(C=C3)[C@@]4(CCN(C4=O)[C@H](CC(C)C)C(=O)NO)N",
    "TMI-005": "CC1([C@@H](N(CCS1)S(=O)(=O)C2=CC=C(C=C2)OCC#CCO)C(=O)NO)C"
}

# 출력 디렉토리
output_dir = "ligands_pdbqt"
os.makedirs(output_dir, exist_ok=True)

for ligand_name, smiles in ligands.items():
    # RDKit을 사용하여 초기 PDB 파일 생성
    ligand = Chem.MolFromSmiles(smiles)
    if ligand is None:
        print(f"Error: Could not create ligand from SMILES: {smiles}")
        continue
    
    # 리간드에 수소 추가
    ligand = Chem.AddHs(ligand)
    
    # 초기 PDB 파일로 저장
    temp_pdb_path = os.path.join(output_dir, f"{ligand_name}_temp.pdb")
    with open(temp_pdb_path, 'w') as f:
        f.write(Chem.MolToPDBBlock(ligand))
    
    # OpenBabel을 사용하여 PDBQT로 변환
    ligand_pdbqt_path = os.path.join(output_dir, f"{ligand_name}.pdbqt")
    subprocess.run(['obabel', temp_pdb_path, '-O', ligand_pdbqt_path])

    print(f"Successfully saved {ligand_name} to {ligand_pdbqt_path}")

