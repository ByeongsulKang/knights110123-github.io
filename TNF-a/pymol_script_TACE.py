import subprocess

# PyMOL을 사용하여 PDB 파일을 확인 및 수정
pymol_script = """
home/kang/TNF-a/TACE.pdb
remove resn HOH
remove solvent
h_add
remove elem O
remove elem C and not (name CA or name CB or name CG or name CD or name CE or name NZ or name NE or name ND or name NE2 or name ND2 or name CZ or name CH or name CH2 or name CD1 or name CD2)
save /home/kang/TNF-a/TACE_modified.pdb
quit
"""

with open("/home/kang/TNF-a/modify_pdb.pml", "w") as file:
    file.write(pymol_script)

subprocess.run(["pymol", "-cq", "/home/kang/TNF-a/modify_pdb.pml"])

