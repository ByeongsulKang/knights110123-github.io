import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# 결과 파일 로드
file_path = '/home/kang/bioinfo_study/240502-PDB/RMSD_result.txt'
df = pd.read_csv(file_path, delimiter=':', names=['pair', 'rmsd'])

# 데이터 정리
df['rmsd'] = df['rmsd'].str.replace(' Å', '').astype(float)
df[['base', 'target']] = df['pair'].str.split(' and ', expand=True)
df.drop(columns=['pair'], inplace=True)

# 데이터 확인
print(df.head())

# 고유 파일 이름 추출
files = np.unique(df[['base', 'target']].values)

# 히트맵 데이터 행렬 초기화
heat_matrix = pd.DataFrame(index=files, columns=files, data=np.nan)

# 히트맵 데이터 채우기
for index, row in df.iterrows():
    heat_matrix.at[row['base'], row['target']] = row['rmsd']
    heat_matrix.at[row['target'], row['base']] = row['rmsd']

# 히트맵 그리기
plt.figure(figsize=(10, 8))
sns.heatmap(heat_matrix, annot=True, fmt=".2f", cmap='viridis')
plt.title('RMSD Heatmap')
plt.show()

