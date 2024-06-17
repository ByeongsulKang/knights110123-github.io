import subprocess

# 리셉터 PDB 파일을 PDBQT 파일로 변환
subprocess.run(["prepare_receptor4.py", "-r", "TACE.pdb", "-o", "TACE.pdbqt"])

# 리셉터 PDBQT 파일 수정: 아연 원자의 Gasteiger 매개변수를 추가
def modify_zinc_parameters(pdbqt_file):
    lines = []
    with open(pdbqt_file, 'r') as file:
        for line in file:
            if "ZN" in line:
                line = line.replace("0.000  0.000  0.000", "1.00  0.00")
            lines.append(line)
    with open(pdbqt_file, 'w') as file:
        file.writelines(lines)

modify_zinc_parameters('TACE.pdbqt')
