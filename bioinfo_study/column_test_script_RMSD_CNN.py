import pandas as pd

# 데이터 로드
data_rmsd = pd.read_csv('/home/kang/bioinfo_study/ligand_docked_rmsd_results0425.csv')
data_affinity = pd.read_csv('/home/kang/bioinfo_study/extracted_data.csv')

# 컬럼 이름 출력
print("RMSD Data Columns:", data_rmsd.columns)
print("Affinity Data Columns:", data_affinity.columns)

