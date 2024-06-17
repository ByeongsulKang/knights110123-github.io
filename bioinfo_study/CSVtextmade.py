import pandas as pd

# 파일 경로
file_path = '/home/kang/bioinfo_study/240425 gnina_cnnscore.txt'

# 데이터 로드, 공백으로 구분된 데이터를 읽어 들이기
# 데이터가 1부터 99번까지만 포함되어야 함
data = pd.read_csv(file_path, delim_whitespace=True, header=None, skiprows=4, nrows=99)

# 컬럼 선택 (세 번째와 네 번째 컬럼)
selected_data = data.iloc[:, [2, 3]]

# 컬럼 이름 설정
selected_data.columns = ['CNN pose score', 'CNN affinity']

# 결과 데이터 CSV 파일로 저장
output_path = '/home/kang/bioinfo_study/extracted_data.csv'
selected_data.to_csv(output_path, index=False)

print(f"Data extracted and saved to {output_path}")

