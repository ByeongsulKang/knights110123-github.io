import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import linregress

# 파일 경로
file_path = '/home/kang/bioinfo_study/scatter_base.csv'

# 데이터 로드
data = pd.read_csv(file_path)

# 데이터 확인 (컬럼 이름 포함)
print("Data Columns:", data.columns)
print(data.head())  # 데이터의 상위 몇 행을 출력하여 실제 데이터 확인

# 산점도 및 회귀선 표시
plt.figure(figsize=(10, 6))
try:
    plot = sns.scatterplot(x='CNN affinity', y='RMSD (nm)', data=data)
    plot.set_title('CNN Affinity vs RMSD')
    plot.set_xlabel('CNN Affinity')
    plot.set_ylabel('RMSD (nm)')

    # 선형 회귀 계산
    slope, intercept, r_value, p_value, std_err = linregress(data['CNN affinity'], data['RMSD (nm)'])
    plt.plot(data['CNN affinity'], intercept + slope * data['CNN affinity'], 'r', label=f'Line: y={slope:.2f}x+{intercept:.2f}')
    plt.legend()

    # 결과 출력
    print(f"Regression line: y = {slope:.2f}x + {intercept:.2f}")
    print(f"Pearson Correlation Coefficient: {r_value:.2f}")
    print(f"P-value: {p_value:.3f}")

    # 플롯 표시
    plt.show()
except Exception as e:
    print(f"An error occurred: {e}")

