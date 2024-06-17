import os

# 리간드 파일이 있는 디렉토리
ligand_dir = "ligands_3D"

# 리간드 파일 목록 가져오기
ligand_files = [f for f in os.listdir(ligand_dir) if f.endswith(".pdb")]

# 파일 이름 출력
for pdb_file in ligand_files:
    print(f"Found file: {pdb_file}")
import os

# 리간드 파일이 있는 디렉토리
ligand_dir = "ligands_3D"

# 리간드 파일 목록 가져오기 및 경로 출력
ligand_files = [f for f in os.listdir(ligand_dir) if f.endswith(".pdb")]

for pdb_file in ligand_files:
    pdb_file_path = os.path.join(ligand_dir, pdb_file)
    print(f"File path: {pdb_file_path}")

