import subprocess

def run_gnina_and_save_output():
    cmd = ['/home/kang/bioinfo_study/gnina',
        '-r', '/home/kang/bioinfo_study/receptor.pdb',
        '-l', '/home/kang/bioinfo_study/ligand.pdb',
        '--autobox_ligand', '/home/kang/bioinfo_study/ligand.pdb',
        '--num_modes', '100',
        '--cnn_scoring', 'all',
        '--out', '/home/kang/bioinfo_study/ligand_docked.pdb'
    ]
    with open('/home/kang/bioinfo_study/gnina_scoring.txt', 'w') as file:
        subprocess.run(cmd, stdout=file)

run_gnina_and_save_output()

