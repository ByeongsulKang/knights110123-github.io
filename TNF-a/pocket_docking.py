from rdkit import Chem
from rdkit.Chem import AllChem
import os
import pandas as pd

# 포켓 파일 디렉토리
pockets_dir = "TACE_out/pockets"
pocket_files = [f for f in os.listdir(pockets_dir) if f.endswith('_atm.pdb')]

# 리간드 SMILES 목록
ligands = {
    "GW3333": "CCC(C)C(C(=O)NC1=CC=CC=N1)NC(=O)C(CC(C)C)C(C(C)C)N(C=O)O",
    "TMI-1": "CC#CCOC1=CC=C(C=C1)S(=O)(=O)N2CCSC([C@@H]2C(=O)NO)(C)C",
    "BMS-561392": "CC1=NC2=CC=CC=C2C(=C1)COC3=CC=C(C=C3)[C@@]4(CCN(C4=O)[C@H](CC(C)C)C(=O)NO)N",
    "TMI-005": "CC1([C@@H](N(CCS1)S(=O)(=O)C2=CC=C(C=C2)OCC#CCO)C(=O)NO)C"
}

# 결과 저장
results = []

for pocket_file in pocket_files:
    pocket_file_path = os.path.join(pockets_dir, pocket_file)
    pocket = Chem.MolFromPDBFile(pocket_file_path, removeHs=False)
    if pocket is None:
        print(f"Error: Could not load pocket from {pocket_file_path}")
        continue
    print(f"Loaded pocket: {pocket_file_path}")
    
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
                    # 최적화 실패 시 구조를 저장하여 문제 파악
                    Chem.MolToMolFile(ligand, f"{ligand_name}_failed.mol")
                    continue
            
            # 포켓과 리간드를 결합
            complex_mol = Chem.CombineMols(pocket, ligand)
            if complex_mol is None:
                print(f"Error: {ligand_name} complex with {pocket_file} could not be created.")
                continue
            
            # 복합체 에너지 계산
            AllChem.EmbedMolecule(complex_mol)
            energy = AllChem.UFFGetMoleculeForceField(complex_mol).CalcEnergy()
            
            pocket_name = pocket_file.replace('_atm.pdb', '')
            results.append({"pocket": pocket_name, "ligand": ligand_name, "energy": energy})
        except Exception as e:
            print(f"An error occurred while processing {ligand_name} with {pocket_file}: {e}")

# 결과를 데이터프레임으로 변환
df_results = pd.DataFrame(results)
df_results.to_csv("/home/kang/TNF-a/docking_results.csv", index=False)

# 결과 출력
print(df_results)

