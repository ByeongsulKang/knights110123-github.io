import pandas as pd

# 파일 경로
path_csv1 = '/home/kang/bioinfo_study/ligand_docked_rmsd_results0425.csv'
path_csv2 = '/home/kang/bioinfo_study/extracted_data.csv'

# 데이터 로드
data1 = pd.read_csv(path_csv1)
data2 = pd.read_csv(path_csv2)

# 데이터 병합
merged_data = pd.concat([data1, data2], axis=1)  # axis=1은 열 방향으로 병합을 의미함

# 병합된 데이터 확인
print(merged_data.head())

# 결과 데이터 CSV 파일로 저장
merged_data.to_csv('/home/kang/bioinfo_study/scatter_base.csv', index=False)

