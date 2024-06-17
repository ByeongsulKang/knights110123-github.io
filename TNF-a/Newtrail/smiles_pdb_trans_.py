# 코드를 실행합니다.
import os

from rdkit import Chem
from rdkit.Chem import AllChem

# 리간드 SMILES 목록
ligands = {
    "GW3333": "CCC(C)C(C(=O)NC1=CC=CC=N1)NC(=O)C(CC(C)C)C(C(C)C)N(C=O)O",
    "TMI-1": "CC#CCOC1=CC=C(C=C1)S(=O)(=O)N2CCSC([C@@H]2C(=O)NO)(C)C",
    "BMS-561392": "CC1=NC2=CC=CC=C2C(=C1)COC3=CC=C(C=C3)[C@@]4(CCN(C4=O)[C@H](CC(C)C)C(=O)NO)N",
    "TMI-005": "CC1([C@@H](N(CCS1)S(=O)(=O)C2=CC=C(C=C2)OCC#CCO)C(=O)NO)C"
}

# 출력 디렉토리
output_dir = "ligands_pdb"
os.makedirs(output_dir, exist_ok=True)

for ligand_name, smiles in ligands.items():
    ligand = Chem.MolFromSmiles(smiles)
    if ligand is None:
        print(f"Error: Could not create ligand from SMILES: {smiles}")
        continue
    
    # 리간드에 수소 추가
    ligand = Chem.AddHs(ligand)
    
    try:
        # 리간드 3D 구조 생성 및 최적화
        AllChem.EmbedMolecule(ligand, AllChem.ETKDG())
        result = AllChem.UFFOptimizeMolecule(ligand)
        if result != 0:
            print(f"Error: {ligand_name} could not be optimized with UFF. Trying MMFF94...")
            result = AllChem.MMFFOptimizeMolecule(ligand)
            if result != 0:
                print(f"Error: {ligand_name} could not be optimized with MMFF94 either. Result code: {result}")
                continue
        
        # PDB 파일로 저장
        pdb_file_path = os.path.join(output_dir, f"{ligand_name}.pdb")
        with open(pdb_file_path, 'w') as f:
            f.write(Chem.MolToPDBBlock(ligand))
        
        print(f"Successfully saved {ligand_name} to {pdb_file_path}")
    except Exception as e:
        print(f"An error occurred while processing {ligand_name}: {e}")

