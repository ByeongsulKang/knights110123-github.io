import subprocess

def run_gnina_and_save_output():
    cmd = [
        '/home/kang/bioinfo_study/gnina',
        '-r','/home/kang/bioinfo_study/receptor.pdb',
        '-l', '/home/kang/bioinfo_study/ligand.pdb',
        '--cnn_scoring', 'all',
        '--out', '/home/kang/bioinfo_study/ligand_docked100-1.pdb'
    ]
    with open('results.txt', 'w') as file:
        subprocess.run(cmd, stdout=file, text=True)

run_gnina_and_save_output()

