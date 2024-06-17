import pymol
from pymol import cmd, stored
import os

def calculate_all_rmsd(directory, output_file):
    """
    Calculates RMSD for all combinations of PDB files in the given directory.
    
    Args:
    - directory (str): Directory containing PDB files.
    - output_file (str): File path fot the output text file.
    """
    # 폴더에서 모든 .pdb 파일 찾기
    pdb_files = [f for f in os.listdir(directory) if f.endswith('.pdb')]
    pdb_files = [os.path.join(directory, f) for f in pdb_files]

    # 초기화
    cmd.reinitialize()

    # 파일 로드 및 RMSD 계산
    rmsd_results = []
    for i, base in enumerate(pdb_files):
        cmd.load(base, "base")
        for target in pdb_files[i+1:]:
            cmd.load(target, "target")
            # 정렬 및 RMSD 계산
            cmd.align("target", "base")
            rmsd = cmd.rms_cur("target and name CA", "base and name CA")
            rmsd_results.append((base, target, rmsd))
            cmd.delete("target")
        cmd.delete("base")
    
    # 결과 출력
    for base, target, rmsd in rmsd_results:
        print(f"RMSD between {os.path.basename(base)} and {os.path.basename(target)}: {rmsd:.2f} Å")
    
    with open(output_file, "w") as file:
        for base, target, rmsd in rmsd_results:
            file.write(f"RMSD between {os.path.basename(base)} and {os.path.basename(target)}: {rmsd:.2f} Å\n")

# 폴더 경로 설정
pdb_directory = "/home/kang/bioinfo_study/240502-PDB"
output_path = "/home/kang/bioinfo_study/240502-PDB/RMSD_result.txt"
# 함수 호출
calculate_all_rmsd(pdb_directory, output_path)

