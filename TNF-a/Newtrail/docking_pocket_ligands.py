import os
import subprocess
import pandas as pd

# 포켓 PDBQT 파일 디렉토리
pockets_dir = "/home/kang/TNF-a/TACE_out/pockets"
pocket_files = [f for f in os.listdir(pockets_dir) if f.endswith('.pdbqt')]

# 리간드 PDBQT 파일 디렉토리
ligands_dir = "ligands_pdbqt"
ligand_files = [f for f in os.listdir(ligands_dir) if f.endswith('.pdbqt')]

# Gnina 실행 함수
def run_gnina(receptor, ligand, out_file, log_file):
    command = [
        'gnina', 
        '--receptor', receptor, 
        '--ligand', ligand, 
        '--out', out_file, 
        '--log', log_file, 
        '--autobox_ligand', ligand,  # 자동으로 도킹 박스를 설정합니다.
        '--exhaustiveness', '8'
    ]
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error running Gnina for {ligand} with {receptor}: {result.stderr}")
    else:
        print(f"Gnina ran successfully for {ligand} with {receptor}")
    return result

# 결과 저장
results = []

# 포켓 PDBQT 파일 리스트
pocket_pdbqt_files = [os.path.join(pockets_dir, f) for f in os.listdir(pockets_dir) if f.endswith('.pdbqt')]

# 리간드 PDBQT 파일 리스트
ligand_pdbqt_files = [os.path.join(ligands_dir, f) for f in os.listdir(ligands_dir) if f.endswith('.pdbqt')]

for pocket_pdbqt_file in pocket_pdbqt_files:
    for ligand_pdbqt_file in ligand_pdbqt_files:
        # 도킹 결과 파일명 설정
        out_file = f"docked_{os.path.splitext(os.path.basename(ligand_pdbqt_file))[0]}_to_{os.path.splitext(os.path.basename(pocket_pdbqt_file))[0]}.sdf"
        log_file = f"log_{os.path.splitext(os.path.basename(ligand_pdbqt_file))[0]}_to_{os.path.splitext(os.path.basename(pocket_pdbqt_file))[0]}.txt"
        
        # Gnina 실행
        result = run_gnina(pocket_pdbqt_file, ligand_pdbqt_file, out_file, log_file)
        
        if result.returncode == 0:
            with open(log_file, 'r') as log:
                lines = log.readlines()
                for line in lines:
                    print(line.strip())  # 로그 파일 내용을 모두 출력하여 친화도 점수가 있는 부분 확인
                    if "CNNscore" in line:
                        parts = line.split()
                        affinity = float(parts[1])
                        results.append({
                            "pocket": os.path.splitext(os.path.basename(pocket_pdbqt_file))[0],
                            "ligand": os.path.splitext(os.path.basename(ligand_pdbqt_file))[0],
                            "affinity": affinity
                        })
                        break
                else:
                    print(f"Warning: No CNNscore found in log file {log_file}")
        else:
            print(f"Error: Gnina failed for {ligand_pdbqt_file} with {pocket_pdbqt_file}. Output: {result.stderr}")

# 결과를 데이터프레임으로 변환
df_results = pd.DataFrame(results)
df_results.to_csv("docking_results_gnina.csv", index=False)

# 결과 출력
print(df_results)

