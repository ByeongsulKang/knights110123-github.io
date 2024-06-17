import os
import subprocess

# 입력 폴더 및 출력 폴더 설정
input_folder = "/home/kang/TNF-a/Newtrail/ligands_pdb"
output_folder = "/home/kang/TNF-a/Newtrail/ligands_pdb"
os.makedirs(output_folder, exist_ok=True)

# OpenBabel을 사용하여 PDB 파일을 PDBQT 파일로 변환하는 함수
def convert_pdb_to_pdbqt(input_pdb, output_pdbqt):
    command = ['obabel', input_pdb, '-O', output_pdbqt, '-xr']
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"Error converting {input_pdb} to PDBQT: {result.stderr}")
    print(f"Successfully converted {input_pdb} to {output_pdbqt}")

# 입력 폴더에서 PDB 파일을 찾아 PDBQT 파일로 변환
for file_name in os.listdir(input_folder):
    if file_name.endswith('.pdb'):
        input_pdb_path = os.path.join(input_folder, file_name)
        output_pdbqt_path = os.path.join(output_folder, f"{os.path.splitext(file_name)[0]}.pdbqt")
        convert_pdb_to_pdbqt(input_pdb_path, output_pdbqt_path)

