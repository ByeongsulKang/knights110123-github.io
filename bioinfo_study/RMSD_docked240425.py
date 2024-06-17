import mdtraj as md
import pandas as pd

def calculate_rmsd(file_path, reference_index=0):
    # 트레젝토리 로드
    traj = md.load(file_path)

    # 참조 포즈 설정
    reference_pose = traj.slice(reference_index)

    # 모든 포즈에 대한 RMSD 계산
    rmsd_values = md.rmsd(traj, reference_pose, frame=reference_index)

    # RMSD 데이터프레임 생성
    data = {
        'Pose': range(1, len(rmsd_values) + 1),
        'RMSD (nm)': rmsd_values
    }
    df = pd.DataFrame(data)

    return df

# 파일 경로와 참조 포즈 인덱스 설정
file_path = '/home/kang/bioinfo_study/ligand_docked.pdb'
reference_index = 0  # 첫 번째 포즈를 참조 포즈로 사용

# RMSD 계산
rmsd_df = calculate_rmsd(file_path, reference_index)

# 결과를 CSV 파일로 저장
csv_file_path = 'ligand_docked_rmsd_results0425.csv'
rmsd_df.to_csv(csv_file_path, index=False)

print(f"RMSD data saved to {csv_file_path}")

