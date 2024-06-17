import subprocess

def setup_paths():
    """경로 설정 함수"""
    paths = {
        'gnina_executable': '/home/kang/bioinfo_study/gnina',
        'receptor': '/home/kang/bioinfo_study/receptor.pdb',
        'ligand': '/home/kang/bioinfo_study/ligand.pdb',
        'output_docked': '/home/kang/bioinfo_study/ligand_docked.pdb',
        'output_scores': '/home/kang/bioinfo_study/gnina_scoring.txt'
    }
    return paths

def run_gnina_docking(paths):
    """도킹 실행 함수"""
    docking_cmd = [
        paths['gnina_executable'],
        '-r', paths['receptor'],
        '-l', paths['ligand'],
        '--autobox_ligand', paths['ligand'],
        '--num_modes', '100',
        '--cnn_scoring',
        '--out', paths['output_docked']
    ]
    subprocess.run(docking_cmd)
    print("도킹이 완료되었습니다.")

def save_gnina_scores(paths):
    """점수 저장 함수"""
    score_cmd = [
        paths['gnina_executable'],
        '-r', paths['receptor'],
        '-l', paths['ligand'],
        '--autobox_ligand', paths['ligand'],
        '--score_only'
    ]
    with open(paths['output_scores'], 'w') as file:
        subprocess.run(score_cmd, stdout=file)
    print("점수가 저장되었습니다.")

def main():
    """메인 실행 함수"""
    paths = setup_paths()
    run_gnina_docking(paths)
    save_gnina_scores(paths)

if __name__ == "__main__":
    main()

