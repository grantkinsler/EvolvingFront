import pandas as p
import os
from Bio import SeqIO

data_directory = 'raw_data'

conditions = ['Evo2D_IRA1_Mis_lane2','Evo2D_IRA1_Non_lane2','Evo3D_IRA1_Mis_lane2','Evo3D_IRA1_Non_lane2','Evo3D_TCG_lane2']

condition_samples = {c:[] for c in conditions}

sample_names = [['Evo2D_IRA1_Mis_Rep1_T17','Evo3D_IRA1_Mis_Rep1_T14','Evo3D_IRA1_Non_Rep1_T14'],
['Evo2D_IRA1_Mis_Rep1_T21','Evo3D_IRA1_Mis_Rep1_T17','Evo3D_IRA1_Non_Rep1_T17'],
['Evo2D_IRA1_Mis_Rep1_T25','Evo3D_IRA1_Mis_Rep1_T21','Evo3D_IRA1_Non_Rep1_T21'],
['Evo2D_IRA1_Mis_Rep2_T17','Evo3D_IRA1_Mis_Rep2_T14','Evo3D_IRA1_Non_Rep2_T14'],
['Evo2D_IRA1_Mis_Rep2_T21','Evo3D_IRA1_Mis_Rep2_T17','Evo3D_IRA1_Non_Rep2_T17'],
['Evo2D_IRA1_Mis_Rep2_T25','Evo3D_IRA1_Mis_Rep2_T21','Evo3D_IRA1_Non_Rep2_T21'],
['Evo2D_IRA1_Non_Rep1_T17','Evo3D_IRA1_Mis_Rep3_T14','Evo3D_IRA1_Non_Rep3_T14'],
['Evo2D_IRA1_Non_Rep1_T21','Evo3D_IRA1_Mis_Rep3_T17','Evo3D_IRA1_Non_Rep3_T17'],
['Evo2D_IRA1_Non_Rep1_T25','Evo3D_IRA1_Mis_Rep3_T21','Evo3D_IRA1_Non_Rep3_T21'],
['Evo2D_IRA1_Non_Rep2_T17','Evo3D_TCG_Rep1_T14','Evo3D_TCG_Rep2_T14'],
['Evo2D_IRA1_Non_Rep2_T21','Evo3D_TCG_Rep1_T17','Evo3D_TCG_Rep2_T17'],
['Evo2D_IRA1_Non_Rep2_T25','Evo3D_TCG_Rep1_T21','Evo3D_TCG_Rep2_T21']]


# sample_key = p.read_csv('18032-44_samplekey.csv')

# admera_id = 

F_indices = ['F201','F202','F203','F204','F205','F206','F207','F208','F209','F210','F211','F212']
R_indices = ['R301','R302','R303','R304','R305','R306','R307','R308']
N_indices = ['N716','N718','N719','N720','N721','N722','N723','N724','N726','N727','N728','N729']
S_indices = ['S513','S515','S516','S517','S518','S520','S521','S522']

mapped_pairs = {}

for fn,row in enumerate(sample_names):
    for rs,entry in enumerate(row):

        for condition in conditions:
            if condition.replace('_lane2','') in entry:
                condition_samples[condition].append(entry)

        this_sample_admera_id = '092319_NextSeq'

        mapped_pairs[entry] = [this_sample_admera_id,F_indices[fn],R_indices[rs],N_indices[fn],S_indices[rs]]


# for condition,samples in condition_samples.items():
#     print(mapped_pairs[samples[0]][0])
#     gzip_list = ' '.join([f'{data_directory}/{mapped_pairs[sample][0]}_{mapped_pairs[sample][3]}{mapped_pairs[sample][4]}_R1.fastq.gz' for sample in samples])
#     print(gzip_list)
#     os.system(f'cat {gzip_list} > {data_directory}/{condition}_merged_R1.fastq.gz')
#     gzip_list = ' '.join([f'{data_directory}/{mapped_pairs[sample][0]}_{mapped_pairs[sample][3]}{mapped_pairs[sample][4]}_R2.fastq.gz' for sample in samples])
#     os.system(f'cat {gzip_list} > {data_directory}/{condition}_merged_R2.fastq.gz')

constant_regions = SeqIO.parse(open('allCRs.fasta'),'fasta')

F_seqs = {}
R_seqs = {}
for region in constant_regions:
    if region.id in F_indices:
        F_seqs[region.id] = region.seq
    elif region.id in R_indices:
        R_seqs[region.id] = region.seq

with open('evolution_tracking_lane2.inp','w') as slurm_input:
    for condition,samples in condition_samples.items():
        slurm_input.write(f'{condition}_merged_R1.fastq.gz\t{condition}_merged_R2.fastq.gz\t{condition}_multiplex.txt\n')


        with open(f'{condition}_multiplex.txt','w') as condition_multiplex:
            for sample in samples:
                if 'Rep' in sample:
                    rep = int(sample.split('Rep')[1][0])
                else:
                    rep = 0

                time = int(sample.split('_')[-1][1:])
                f_index_seq = F_seqs[mapped_pairs[sample][1]][:6]
                r_index_seq = R_seqs[mapped_pairs[sample][2]][-8:]
                condition_multiplex.write(f'{condition}\t{rep}\t{time}\t{f_index_seq}\t{r_index_seq}\n')







