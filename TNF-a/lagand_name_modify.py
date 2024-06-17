import os
import re

# 리간드 파일 목록
ligand_dir = "/home/kang/TNF-a/ligands_3D"
ligand_files = [os.path.join(ligand_dir, f) for f in os.listdir(ligand_dir) if f.endswith(".pdb")]

# 파일 이름 정리 함수
def clean_filename(filename):
    return re.sub(r'[^a-zA-Z0-9_]', '_', filename)

# 파일 이름 변경
for pdb_file in ligand_files:
    base_name = os.path.basename(pdb_file)
    clean_name = clean_filename(base_name)
    new_path = os.path.join(ligand_dir, clean_name)
    os.rename(pdb_file, new_path.replace('.pdb', '_clean.pdb'))

