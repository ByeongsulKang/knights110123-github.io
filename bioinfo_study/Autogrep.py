import sys

def process_pdb(pdb_filename):
    # 입력 파일 이름에서 출력 파일 이름 생성
    output_filename = f"{pdb_filename.split('.')[0]}_ADP.pdb"

    # 파일 처리
    with open(pdb_filename, 'r') as file:
        with open(output_filename, 'w') as outfile:
            for line in file:
                if line.startswith("ATOM") and " A " in line:
                    outfile.write(line)
            # 'TER' 추가
            outfile.write("TER\n")
            # 파일 포인터를 다시 시작으로 옮김
            file.seek(0)
            for line in file:
                if line.startswith("HETATM") and "ADP" in line and " A " in line:
                    outfile.write(line)
            # 'END' 추가
            outfile.write("END\n")

    print(f"Processed file saved as {output_filename}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <pdb_filename>")
        sys.exit(1)
    process_pdb(sys.argv[1])
