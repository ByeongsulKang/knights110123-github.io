# 파일의 첫 5줄을 출력하여 컬럼 확인
with open('240425 gnina_cnnscore.txt', 'r') as file:
    for _ in range(5):
        print(file.readline())

