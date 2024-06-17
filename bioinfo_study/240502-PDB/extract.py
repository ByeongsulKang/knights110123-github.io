import pandas as pd
import umap.umap_ as umap
import matplotlib.pyplot as plt

def parse_rmsd_results(file_path):
    # 결과를 저장할 리스트 초기화
    data = []

    # 파일 열기 및 줄별로 처리
    with open(file_path, 'r') as file:
        for line in file:
            # 필요한 정보 추출
            if line.startswith("RMSD between"):
                parts = line.split()
                file1 = parts[2]
                file2 = parts[4]
                rmsd_value = float(parts[-2])  # 마지막에서 두 번째 요소는 RMSD 값

                # 데이터 리스트에 사전 형태로 추가
                data.append({
                    'File1': file1,
                    'File2': file2,
                    'RMSD (Å)': rmsd_value
                })

    # 데이터 프레임 생성
    df = pd.DataFrame(data)

    # 데이터 프레임을 피벗 테이블로 변환
    pivot_df = df.pivot(index='File1', columns='File2', values='RMSD (Å)')
   
    pivot_df = pivot_df.clip(upper=100)
   
   # NaN 값을 100으로 채우기
    pivot_df.fillna(100, inplace=True)
    
    # RMSD 값을 100 - RMSD로 변환
    pivot_df = 100 - pivot_df

    # 피벗 테이블 출력
    print(pivot_df)

    # 결과를 CSV 파일로 저장
    pivot_df.to_csv('rmsd_pivot_results_transformed.csv')


# UMAP 임베딩
    reducer = umap.UMAP(random_state=42)
    embedding = reducer.fit_transform(pivot_df)

    # 임베딩 결과 시각화
    plt.figure(figsize=(12, 10))
    plt.scatter(embedding[:, 0], embedding[:, 1], c='blue', cmap='Spectral', s=50)
    plt.title('UMAP Projection of RMSD Data')
    plt.xlabel('UMAP Component 1')
    plt.ylabel('UMAP Component 2')
    plt.colorbar()
    plt.show()

# 파일 경로 설정
file_path = '/home/kang/bioinfo_study/240502-PDB/RMSD_result.txt'  # 파일 경로를 실제 파일 위치로 변경하세요.

# 함수 호출
parse_rmsd_results(file_path)

