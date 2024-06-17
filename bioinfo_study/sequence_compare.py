import os
import csv
from Bio import PDB
from Bio.pairwise2 import align as pairwise2_align
from Bio.SeqUtils import seq1

def extract_sequence(pdb_file):
    parser = PDB.PDBParser()
    structure = parser.get_structure('PDB', pdb_file)
    for model in structure:
        for chain in model:
            sequence = []
            for residue in chain.get_residues():
                if PDB.is_aa(residue, standard=True):
                    sequence.append(seq1(residue.resname))
            return ''.join(sequence)

def compare_sequences(folder_path, output_csv):
    sequences = {}
    pdb_files = [f for f in os.listdir(folder_path) if f.endswith('.pdb')]

    # 서열 추출
    for pdb_file in pdb_files:
        path = os.path.join(folder_path, pdb_file)
        seq = extract_sequence(path)
        sequences[pdb_file] = seq
    
    # 결과를 CSV 파일에 저장
    with open(output_csv, 'w', newline='') as csvfile:
        fieldnames = ['File1', 'File2', 'Score']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        # 서열 비교
        for file1, seq1 in sequences.items():
            for file2, seq2 in sequences.items():
                if file1 < file2:  # 중복 비교 방지
                    alignments = pairwise2_align.globalxx(seq1, seq2)
                    if alignments:  # 비어 있지 않은지 확인
                        score = alignments[0].score
                    else:
                        score = 0  # 유사도가 전혀 없는 경우, 점수를 0으로 처리
                    writer.writerow({'File1': file1, 'File2': file2, 'Score': score})
                    print(f"Comparison between {file1} and {file2}: Score = {score}")

# 폴더 경로와 출력 CSV 파일 경로 설정
folder_path = '/home/kang/bioinfo_study/VP2_PDB'  # 실제 PDB 파일이 저장된 경로로 변경
output_csv = '/home/kang/bioinfo_study/VP2_PDB/sequence_compare.csv'  # 결과 CSV 파일의 경로

# 함수 호출
compare_sequences(folder_path, output_csv)

