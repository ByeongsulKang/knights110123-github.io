import pandas as pd
import numpy as np
import math

def calculate_rmsd(file1, file2):
    with open(file1, 'r') as f1, open(file2, 'r') as f2:
        lines1 = [line.strip() for line in f1.readlines() if line.strip() and not line.startswith(('OpenBabel'))]
        lines2 = [line.strip() for line in f2.readlines() if line.strip() and not line.startswith(('OpenBabel'))]

    distances_squared = []
    for line1, line2 in zip(lines1[2:], lines2[2:]):  # 첫 두 줄은 헤더 정보를 건너뜁니다.
        parts1 = line1.split()
        parts2 = line2.split()
        if len(parts1) > 3 and len(parts2) > 3:  # 좌표 데이터가 있는지 확인
            x1, y1, z1 = map(float, parts1[:3])
            x2, y2, z2 = map(float, parts2[:3])
            distances_squared.append((x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2)
    
    if distances_squared:
        rmsd = math.sqrt(sum(distances_squared) / len(distances_squared))
    else:
        rmsd = -1  # 유효한 좌표 데이터가 없는 경우
    return rmsd

# 파일 목록 생성
file_list = [f"FMM2_{i}.mol" for i in range(1, 101)]

# RMSD 결과를 저장할 빈 리스트
rmsd_results = []

# 모든 파일 쌍에 대해 RMSD 계산
for i in range(len(file_list)):
    for j in range(i+1, len(file_list)):
        rmsd = calculate_rmsd(file_list[i], file_list[j])
        rmsd_results.append({'File1': file_list[i], 'File2': file_list[j], 'RMSD': rmsd})

# 결과를 pandas DataFrame으로 변환
df = pd.DataFrame(rmsd_results)

# DataFrame을 CSV 파일로 저장
df.to_csv('rmsd_results.csv', index=False)

print("RMSD calculation completed and results saved to rmsd_results.csv.")

