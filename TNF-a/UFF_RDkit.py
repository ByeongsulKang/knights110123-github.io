from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolAlign
import pandas as pd

# 수정된 리셉터 파일 로드
receptor_file = "/home/kang/TNF-a/TACE_modified.pdb"
receptor = Chem.MolFromPDBFile(receptor_file, removeHs=False)
if receptor is None:
    raise ValueError(f"Error: Could not load receptor from {receptor_file}")

# 리간드 SMILES 목록
ligands = {
    "GW3333": "CC(C)C1CC(C(=O)NC1=O)NC(=O)C(CC2=CC=CC=C2)N",
    "TMI-1": "CC(C)C1CC(C(=O)NC1=O)NC(=O)C(CC2=CC=CC=C2)NCC#C",
    "BMS-561392 (DPC-333)": "CC(C)C1CC(C(=O)NC1=O)NC(=O)C(CC2=CC=CC=C2)NC3=NC=CC(=N3)C4=CC=CC=C4",
    "TMI-2": "CC(C)C1CC(C(=O)NC1=O)NC(=O)C(CC2=CC=CC=C2)NCC(C3=CC=CC=C3)O",
    "BMS-566394": "CC(C)C1CC(C(=O)NC1=O)NC(=O)C(CC2=CC=CC=C2)NC3=NC=CC(=N3)C(C4=CC=CC=C4)OC",
    "TMI-005 (apratastat)": "CC(C)C1CC(C(=O)NC1=O)NC(=O)C(CC2=CC=CC=C2)NCC#CC",
    "GW4459": "CC(C)C1CC(C(=O)NC1=O)NC(=O)C(CC2=CC=CC=C2)NCC#N",
    "Compound #51 in Duan et al. (2007)": "CC(C)C1CC(C(=O)NC1=O)NC(=O)C(CC2=CC=CC=C2)NCC#CO",
    "W-3646": "CC(C)C1CC(C(=O)NC1=O)NC(=O)C(CC2=CC=CC=C2)NCC#CN"
}

# 결과 저장
results = []

for name, smiles in ligands.items():
    # 리간드 3D 구조 생성 및 최적화
    ligand = Chem.MolFromSmiles(smiles)
    if ligand is None:
        print(f"Error: {name} could not be converted to a molecule.")
        continue
    ligand = Chem.AddHs(ligand)
    AllChem.EmbedMolecule(ligand, AllChem.ETKDG())
    result = AllChem.UFFOptimizeMolecule(ligand)
    if result != 0:
        print(f"Error: {name} could not be optimized.")
        continue
    
    # 리셉터와 리간드를 결합
    complex_mol = Chem.CombineMols(receptor, ligand)
    if complex_mol is None:
        print(f"Error: {name} complex could not be created.")
        continue
    
    # 복합체 에너지 계산
    AllChem.EmbedMolecule(complex_mol)
    energy = AllChem.UFFGetMoleculeForceField(complex_mol).CalcEnergy()
    
    results.append({"ligand": name, "energy": energy})

# 결과를 데이터프레임으로 변환
df_results = pd.DataFrame(results)
df_results.to_csv("docking_results.csv", index=False)

# 결과 출력
print(df_results)

