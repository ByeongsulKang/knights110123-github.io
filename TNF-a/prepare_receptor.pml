from pymol import cmd

# 리셉터 파일 처리
def prepare_receptor():
    # Load the receptor PDB file
    cmd.load("TACE.pdb")

    # Remove zinc atoms if present
    cmd.select("zinc", "resn ZN")
    cmd.remove("zinc")

    # Add hydrogens
    cmd.h_add()

    # Add Gasteiger charges
    cmd.set("pqr_charges")

    # Save as PDBQT file
    cmd.save("TACE_G.pdbqt")

# 리셉터 파일 준비 함수 실행
prepare_receptor()
cmd.quit()

