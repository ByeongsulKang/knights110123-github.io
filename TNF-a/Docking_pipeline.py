from rdkit import Chem
from rdkit.Chem import AllChem
import os
import pandas as pd

# 수정된 리셉터 파일 로드
receptor_file = "/home/kang/TNF-a/TACE_receptor.pdb"
receptor = Chem.MolFromPDBFile(receptor_file, removeHs=False)
if receptor is None:
    raise ValueError(f"Error: Could not load receptor from {receptor_file}")

# 리간드 파일들이 저장된 디렉토리
ligands_dir = "/home/kang/TNF-a/ligands_3D"

# 결과 저장
results = []

# HETATM 레코드가 있는 리간드 파일들
ligand_files = [f for f in os.listdir(ligands_dir) if f.startswith('HETATM_') and f.endswith('.pdb')]

for ligand_file in ligand_files:
    ligand_file_path = os.path.join(ligands_dir, ligand_file)
    ligand = Chem.MolFromPDBFile(ligand_file_path, removeHs=False, sanitize=False)
    if ligand is None:
        print(f"Error: Could not load ligand from {ligand_file_path}")
        continue
    
    # 리간드 최적화
    AllChem.EmbedMolecule(ligand, AllChem.ETKDG())
    result = AllChem.UFFOptimizeMolecule(ligand)
    if result != 0:
        print(f"Error: {ligand_file} could not be optimized.")
        continue
    
    # 리셉터와 리간드를 결합
    complex_mol = Chem.CombineMols(receptor, ligand)
    if complex_mol is None:
        print(f"Error: {ligand_file} complex could not be created.")
        continue
    
    # 복합체 에너지 계산
    AllChem.EmbedMolecule(complex_mol)
    energy = AllChem.UFFGetMoleculeForceField(complex_mol).CalcEnergy()
    
    ligand_name = ligand_file.replace('HETATM_', '').replace('.pdb', '')
    results.append({"ligand": ligand_name, "energy": energy})

# 결과를 데이터프레임으로 변환
df_results = pd.DataFrame(results)
df_results.to_csv("/home/kang/TNF-a/docking_results.csv", index=False)

# 결과 출력
print(df_results)

