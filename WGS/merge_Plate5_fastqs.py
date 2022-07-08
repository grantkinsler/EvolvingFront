from glob import glob

import subprocess

data_dir = 'AllRawData/'
prefix_file = 'EvolvingFront_WGS_Plate5_prefixes.inp'
merged_dir = 'Plate5_merged/'

with open(prefix_file,'r') as f:
    for prefix in f:
        for read in [1,2]:
            prefix = prefix.replace('\n','')
            glob_out = ' '.join(glob(f'{data_dir}{prefix}*R{read}*'))
            print(glob_out)
            subprocess.run([f"cat {glob_out} > {data_dir}{merged_dir}{prefix}_LanesMerged_R{read}.fastq.gz"],shell=True)