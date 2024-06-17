import subprocess

# 도킹 실행
def run_gnina_docking():
    docking_cmd = [
        '/home/kang/bioinfo_study/gnina',
        '-r', 'receptor.pdb',
        '-l', 'ligand.pdb',
        '--autobox_ligand', 'ligand.pdb',
        '--num_modes', '100',
        '--cnn_scoring',
        '--out', 'ligand_docked.pdb'
    ]
    subprocess.run(docking_cmd)

# 점수 저장
def save_gnina_scores():
    score_cmd = 'gnina -r receptor.pdb -l ligand.pdb --autobox_ligand ligand.pdb --score_only > gnina_scoring.txt'
    subprocess.run(score_cmd, shell=True)

# 함수 실행
run_gnina_docking()
save_gnina_scores()

