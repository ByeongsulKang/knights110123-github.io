import os
import csv
from Bio import PDB

def calculate_rmsd(pdb_file1, pdb_file2):
    parser = PDB.PDBParser()
    structure1 = parser.get_structure('PDB1', pdb_file1)
    structure2 = parser.get_structure('PDB2', pdb_file2)

    # 첫 번째 모델과 체인 추출
    ref_model = structure1[0]
    sample_model = structure2[0]

    ref_chain = next(ref_model.get_chains())
    sample_chain = next(sample_model.get_chains())

    # CA 원자만 선택
    ref_atoms = [atom for atom in ref_chain.get_atoms() if atom.name == 'CA']
    sample_atoms = [atom for atom in sample_chain.get_atoms() if atom.name == 'CA']

    # 일치하는 원자만을 사용
    ref_atoms = sorted(ref_atoms, key=lambda x: x.get_id())
    sample_atoms = sorted(sample_atoms, key=lambda x: x.get_id())
    
    min_length = min(len(ref_atoms), len(sample_atoms))
    ref_atoms = ref_atoms[:min_length]
    sample_atoms = sample_atoms[:min_length]

    # RMSD 계산
    if ref_atoms and sample_atoms and len(ref_atoms) == len(sample_atoms):
        super_imposer = PDB.Superimposer()
        super_imposer.set_atoms(ref_atoms, sample_atoms)
        super_imposer.apply(sample_model.get_atoms())
        return super_imposer.rms
    else:
        return 100  # 적절한 원자가 매칭되지 않는 경우

def compare_all_pdbs(folder_path, output_csv):
    pdb_files = [os.path.join(folder_path, f) for f in os.listdir(folder_path) if f.endswith('.pdb')]
    with open(output_csv, 'w', newline='') as csvfile:
        fieldnames = ['File1', 'File2', 'RMSD']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for i in range(len(pdb_files)):
            for j in range(i + 1, len(pdb_files)):
                rmsd = calculate_rmsd(pdb_files[i], pdb_files[j])
                if rmsd is not None:
                    writer.writerow({'File1': os.path.basename(pdb_files[i]), 'File2': os.path.basename(pdb_files[j]), 'RMSD': rmsd})
                    print(f'RMSD between {os.path.basename(pdb_files[i])} and {os.path.basename(pdb_files[j])}: {rmsd}')
                else:
                    print(f'No valid atom match between {os.path.basename(pdb_files[i])} and {os.path.basename(pdb_files[j])}')

# 폴더 경로와 출력 CSV 파일 경로 설정
folder_path = '/home/kang/bioinfo_study/VP2_PDB'  # 실제 PDB 파일이 저장된 경로로 변경
output_csv = '/home/kang/bioinfo_study/VP2_PDB/RMSD_results.csv'  # 결과 CSV 파일의 경로

# 함수 호출
compare_all_pdbs(folder_path, output_csv)

