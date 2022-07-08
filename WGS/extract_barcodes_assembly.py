################# 
### Script to extract barcode regions from Whole Genome Sequencing Data.
### THIS ASSEMBLES ALL THE FOUND BARCODE REGIONS
### Grant Kinsler
################# 

import os
import sys
import re
import numpy as np
import pandas as p

################# 
### Extract the Barcode Region
#################

bc_directory = "barcode_extraction/"
# data_directory = "AllRawData/"

list_file = sys.argv[1]
file_list = []
with open(list_file,'r') as list_f:
    for line in list_f:
        file_list.append(line.split('\n')[0])

bc_dictionary = {}

for file_prefix in file_list:
    this_file = f"{bc_directory}{file_prefix}_extractedBCs.out"
    if os.path.exists(this_file):
        with open(this_file,'r') as barcode_extract:
            for line in barcode_extract:
                bc_dictionary[file_prefix] = line.split('\t')

    else:
        bc_dictionary[file_prefix] = ['','']

with open(f'{bc_directory}allExtractedBCs.csv','w') as out:
    out.write('file_prefix,doubleBC,BC1,BC2\n')
    for file_prefix,bc_list in bc_dictionary.items():

        out.write(f'{file_prefix},{"_".join(bc_list)},{bc_list[0]},{bc_list[1]}\n')





