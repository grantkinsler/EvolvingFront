import pandas as p
import os

data_dir = 'AllRawData/18032-115/'
rename_df = p.read_csv('rename_sequencing_files.csv')

rename_dict = {file_prefix:sample for file_prefix,sample in zip(rename_df['file_prefix'].values,rename_df['sample'].values)}

for file_prefix,sample in rename_dict.items():
    for read in [1,2]:
        if os.path.isfile(f'{data_dir}{file_prefix}_R{read}.fastq.gz'):
            os.system(f'mv {data_dir}{file_prefix}_R{read}.fastq.gz {data_dir}{sample}_R{read}.fastq.gz')
        else:
            print(f'{data_dir}{file_prefix}_R{read}.fastq.gz does not exist')