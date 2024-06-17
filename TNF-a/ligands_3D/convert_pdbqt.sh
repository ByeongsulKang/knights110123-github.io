# 리간드 PDB 파일을 PDBQT 파일로 변환하는 스크립트
#!/bin/bash

for pdb_file in ligands_3D/*.pdb; do
    pdbqt_file="${pdb_file%.pdb}.pdbqt"
    prepare_ligand4.py -l "$pdb_file" -o "$pdbqt_file"
done

