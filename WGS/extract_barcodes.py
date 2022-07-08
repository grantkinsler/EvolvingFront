################# 
### Script to extract barcode regions from Whole Genome Sequencing Data.
### Grant Kinsler
################# 

import os
import sys
import re
import numpy as np
from Bio.Seq import Seq

################# 
### Extract the Barcode Region
#################

bc_directory = "barcode_extraction/"
data_directory = "AllRawData/"

file_prefix = sys.argv[1]
# file_prefix = "18032FL-96-07-03"

count_threshold = 4
count_ratio = 4

bc1_forward_flanking_seq = "GTCGACGGATCCGATATCGGTACC"
loxP_sequence = "ATAACTTCGTATAATGTATGCTATACGAAGTTAT" ## loxP site between the two barcodes
bc2_reverse_flanking_seq = "GGTACCGATATCAGATCTAAGCTT" 

def reverse_complement(sequence):

    return str(Seq(sequence).reverse_complement())

## this greps the seqeuncing for the loxP site in the middle of the barcode region 
## iterates through both reads and appends to the same file
for read_number in [1,2]:
    # os.system(f"zcat {data_directory}{file_prefix}_R{read_number}.fastq.gz | grep -B 1 -A 2 {bc1_reverse_flanking_seq} >> {bc_directory}{file_prefix}_BC1_extract.out")
    os.system(f"zcat {data_directory}{file_prefix}_R{read_number}.fastq.gz | grep {loxP_sequence} >> {bc_directory}{file_prefix}_BC_forward_reads.out")
    os.system(f"zcat {data_directory}{file_prefix}_R{read_number}.fastq.gz | grep {reverse_complement(loxP_sequence)} >> {bc_directory}{file_prefix}_BC_reverse_reads.out")

### Throw out reads due to amplicon sequencing contamination (barcode region starts at exactly same location as in amplicons)
### Default parameters are for "forward direction" (in orientation of URA3 gene) and accounts for normal length and Ancestor with ApaLI site
def amplicon_contamination(read,seq=loxP_sequence,positions_to_toss=[89,95]):

    flanking_position = read.find(seq)

    if flanking_position in positions_to_toss:
        return True
    else:
        return False

### Figure out what the sequences are... (map to the actual barcode region, extract relevant sequences)

# extract the sequence bewteen the flanking regions using regular expressions.
regex_pattern_bc1 = rf"{bc1_forward_flanking_seq}(.*?){loxP_sequence}"
regex_pattern_bc1_rc = rf"{reverse_complement(loxP_sequence)}(.*?){reverse_complement(bc1_forward_flanking_seq)}"

regex_pattern_bc2 = rf"{loxP_sequence}(.*?){bc2_reverse_flanking_seq}"
regex_pattern_bc2_rc = rf"{reverse_complement(bc2_reverse_flanking_seq)}(.*?){loxP_sequence}"

bc1_sequences = []
bc2_sequences = []
with open(f"{bc_directory}{file_prefix}_BC_forward_reads.out",'r') as forward_reads:
    for read in forward_reads:
        if not amplicon_contamination(read,seq=loxP_sequence,positions_to_toss=[89,95]):
            extracted_bc1s = re.findall(regex_pattern_bc1,read)
            if len(extracted_bc1s) == 1:
                bc1_sequences.append(extracted_bc1s[0])

            extracted_bc2s = re.findall(regex_pattern_bc2,read)
            if len(extracted_bc2s) == 1:
                bc2_sequences.append(extracted_bc2s[0])


# look at files where the barcode region showed up in the reverse complement
with open(f"{bc_directory}{file_prefix}_BC_reverse_reads.out",'r') as reverse_reads:
    for read in reverse_reads:
        if not amplicon_contamination(read,seq=reverse_complement(loxP_sequence),positions_to_toss=[75]):
            extracted_bc1s = re.findall(regex_pattern_bc1_rc,read)
            if len(extracted_bc1s) == 1:
                bc1_sequences.append(reverse_complement(extracted_bc1s[0]))

            extracted_bc2s = re.findall(regex_pattern_bc2_rc,read)
            if len(extracted_bc2s) == 1:
                bc2_sequences.append(reverse_complement(extracted_bc2s[0]))

################
## Identify Barcode Sequences
################
### identify barcodes -> needs to be above count threshold and at least count ratio X the second highest read
def identify_barcode(bc_sequences,count_threshold=count_threshold,count_ratio=count_ratio):
    
    bc_seqs, bc_counts = np.unique(bc_sequences,return_counts=True)

    index_sort = np.argsort(-bc_counts)
    bc_seqs = bc_seqs[index_sort]
    bc_counts = bc_counts[index_sort]

    if len(bc_seqs) > 0:
        if len(bc_seqs) > 1:
            if (bc_counts[0]/bc_counts[1] > count_ratio):
                chosen_bc = bc_seqs[0]
            else:
                chosen_bc = False
                print('no chosen BC')

        elif bc_counts[0] >= count_threshold:
            chosen_bc = bc_seqs[0]
        else:
            chosen_bc = False
            print('no chosen BC')
    else:
        chosen_bc = False
        print('no BC reads...')

    return chosen_bc, bc_seqs, bc_counts




chosen_bc1,bc1_seqs,bc1_counts = identify_barcode(bc1_sequences)

with open(f"{bc_directory}{file_prefix}_BC1_seqs.out",'w') as out:
    for seq,count in zip(bc1_seqs,bc1_counts):
        out.write(f'{seq}\t{count}\n')


chosen_bc2,bc2_seqs,bc2_counts = identify_barcode(bc2_sequences)

with open(f"{bc_directory}{file_prefix}_BC2_seqs.out",'w') as out:
    for seq,count in zip(bc2_seqs,bc2_counts):
        out.write(f'{seq}\t{count}\n')

################
## Output extracted barcodes
################

if chosen_bc1 and chosen_bc2:
    with open(f"{bc_directory}{file_prefix}_extractedBCs.out",'w') as out:
        out.write(f'{chosen_bc1}\t{chosen_bc2}')

