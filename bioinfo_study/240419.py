def extract_structure(input_file, output_file, keywords):
    """특정 키워드를 포함하는 줄만 추출하여 새 파일에 저장"""
    with open(input_file, 'r') as file:
        lines = file.readlines()
    with open(output_file, 'w') as file:
        for line in lines:
            # 키워드 리스트 중 하나라도 해당 라인에 포함되어 있는지 확인
            if any(keyword in line for keyword in keywords):
                file.write(line)

# 입력 파일 이름을 '1kxx.pdb'로 설정
input_pdb_file = '1KXX.pdb'

# 리셉터 구조 추출
# 'ATOM'으로 시작하는 라인만 추출하여 'receptor.pdb'로 저장
extract_structure(input_pdb_file, 'receptor.pdb', ['ATOM'])

# 리간드 구조 추출
# 'HETATM'으로 시작하며 'LIG' (리간드 이름, 실제 이름으로 변경 필요)를 포함하는 라인만 추출하여 'ligand.pdb'로 저장
extract_structure(input_pdb_file, 'ligand.pdb', ['HETATM', 'LIG'])  # 'LIG'는 실제 리간드 이름으로 변경

