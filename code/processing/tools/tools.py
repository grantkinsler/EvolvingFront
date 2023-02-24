import numpy as np
import pandas as p
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import to_hex
from scipy.stats import pearsonr
from scipy.spatial import distance
from scipy.stats.mstats import gmean
from multiprocessing import Pool
from sklearn.linear_model import LinearRegression
from scipy.spatial.distance import cdist, euclidean
from itertools import combinations
import itertools
import copy


ancestral_mutations = {'IRA1_NON':'II:522427:A:T:IRA1:stop_gained:c.4202T>A:p.Leu1401*',
                      'IRA1_MIS':'II:522697:G:A:IRA1:missense_variant:c.3932C>T:p.Ala1311Val',
                      'CYR1':'X:427906:C:A:CYR1:missense_variant:c.2750C>A:p.Ser917Tyr',
                      'GPB2':'I:40104:T:G:GPB2:stop_gained:c.846T>G:p.Tyr282*',
                      'TOR1':'X:564551:T:G:TOR1:missense_variant:c.5136T>G:p.Phe1712Leu'}

anc_color_map = {'WT':'k',
             'GPB2':'#4daf4a',
             'CYR1':'#e41a1c',
             'TOR1':'#6a51a3',
             'IRA1_MIS':'#02818a',
             'IRA1_NON':'#1f78b4',
             'all_second':'orange',
            }

# anc_color_map = {'WT':'k',
#              'GPB2':'#4daf4a',
#              'CYR1':'#e41a1c',
#              'TOR1':'#6a51a3',
#              'IRA1_MIS':'#02818a',
#              'IRA1_NON':'#1f78b4',
#             }

evo_cond_marker_map = {'Evo1D':'o',
              'Evo2D':'^',
              'Evo3D':'s',
              'Evo5D':'p',
              'Evo1_5D':'D',
              'unknown':'x',
             }

ploidy_marker_map = {'haploid':'o','diploid':'^'}

anc_evo_cond_color_map = {'WT':{'Evo1D':'#cccccc','Evo2D':'#252525','Evo5D':'#969696','Evo1_5D':'#636363'},
             'GPB2':{'Evo1D':'#bae4b3','Evo2D':'#74c476','Evo3D':'#238b45','unknown':'#bae4b3'},
             'CYR1':{'Evo1D':'#fcae91','Evo2D':'#fb6a4a','Evo3D':'#cb181d','unknown':'#fcae91'},
             'TOR1':{'Evo1D':'#cbc9e2','Evo2D':'#9e9ac8','Evo3D':'#6a51a3','unknown':'#cbc9e2'},
             'IRA1_MIS':{'Evo1D':'#bdc9e1','Evo2D':'#67a9cf','Evo3D':'#02818a','unknown':'#bdc9e1'},
             'IRA1_NON':{'Evo1D':'#bdd7e7','Evo2D':'#6baed6','Evo3D':'#2171b5','unknown':'#bdd7e7'}
            }



color_map = {        
#              'WT':{'Evo1D':'#cccccc','Evo2D':'#f7f7f7','Evo5D':'#969696','Evo1_5D':'#636363'},
             'WT':{'Evo1D':'#cccccc','Evo2D':'#252525','Evo5D':'#969696','Evo1_5D':'#636363'},
             'GPB2':{'Evo1D':'#bae4b3','Evo2D':'#74c476','Evo3D':'#238b45','unknown':'#bae4b3'},
             'CYR1':{'Evo1D':'#fcae91','Evo2D':'#fb6a4a','Evo3D':'#cb181d','unknown':'#fcae91'},
             'TOR1':{'Evo1D':'#cbc9e2','Evo2D':'#9e9ac8','Evo3D':'#6a51a3','unknown':'#cbc9e2'},
             'IRA1_MIS':{'Evo1D':'#bdc9e1','Evo2D':'#67a9cf','Evo3D':'#02818a','unknown':'#bdc9e1'},
             'IRA1_NON':{'Evo1D':'#bdd7e7','Evo2D':'#6baed6','Evo3D':'#2171b5','unknown':'#bdd7e7'}
            }

rebarcoding_source_mutants = {
'IRA1_MIS':'CGCTAAAGACATAATGTGGTTTGTTG_CTTCCAACAAAAAATCATTTTTATAC', # BCID 43361 from venkataram 2016
'IRA1_NON':'CGCTAAAGACATAATGTGGTTTGTTG_AGAGTAATCTGCAAGATTCTTTTTCT', # BCID 21967 from venkataram 2016
'CYR1':    'CGCTAAAGACATAATGTGGTTTGTTG_CTCGAAACAGGAAAAGCACTTATCGA', # BCID 43692 from venkataram 2016
'TOR1':    'CGCTAAAGACATAATGTGGTTTGTTG_TAGACAAAATGCAATTGTATTGTCAG' , # BCID 21543 from venkataram 2016
'GPB2':    'CGCTAAAGACATAATGTGGTTTGTTG_TCATGAACGGATAAGCTGGTTGGTTG' } # BCID 7774 from venkataram 2016



long_colors = list(sns.color_palette())+list(sns.color_palette("Set2"))

gene_pathway_map = {
    'Ras/PKA':['IRA1','IRA2','GPB1','GPB2','PDE2','CYR1','GPR1','SHR5'],
    'TOR/Sch9':['TOR1','SCH9','KSP1'],
    'HOG':['HOG1','PBS2','SSK2'],
    'RTG':['RTG2','MKS1','BMH1','BMH2'],
    'TCA cycle':['CIT1','KGD1','MDH1','MAE1','ALD5'],
    'Deadenylation':['PUF3','PAB1','PAN2','PAN3']}

# long_colors += ['gray']*(len(all_genes_sorted)-len(long_colors))

# mutation_color_map = {gene:long_colors[g] for g,gene in enumerate(list(all_genes_sorted[:16]))}
# mutation_color_map = {
#     'KSP1':long_colors[4], # purple (b/c TOR pathway)
#     'PUF3':long_colors[1], # orange
#     'PAB1':long_colors[3], # red
#     'RTG2':long_colors[5], # brown
#     'CIT1':long_colors[11], # pinky orange
#     'ARO80':long_colors[6], # pink
#     'GSH1':long_colors[7], # absolute gray
#     'MKS1':long_colors[13], # light pink
#     'SSK2':long_colors[9], # teal
#     'MKT1':long_colors[8], # greenish yellow
#     'MIT1':long_colors[10], 
#     'GPB2':long_colors[2], # green
#     'KGD1':long_colors[12], 
#     'MAE1':long_colors[14], # light green
#     'MDH1':long_colors[14], # light green
#     'IRA1':long_colors[0], # blue (obviously)
#     'IRA2':long_colors[0], # same color as Ira1
#     'double_mutant':'k'
#     }


mutation_color_map = {
    'KSP1':long_colors[4], # purple (b/c TOR pathway)
    'PUF3':long_colors[1], # orange
    'PAB1':long_colors[3], # red
    
    # browns 
    'RTG2':'#8c510a', # dark brown
    'MKS1':'#bf812d', # medium brown 
    
    ### TCA cycle mutants green
    'CIT1':'#006d2c', # darkest
    'KGD1':'#31a354',
    'MAE1':'#74c476',
    'MDH1':'#bae4b3', # lightest
    'MDH2':'#bae4b3', # lightest


    'ARO80':long_colors[6], # pink
    'GSH1':long_colors[7], # absolute gray

    'SSK2':long_colors[9], # teal
    'MKT1':long_colors[8], # greenish yellow
    'MIT1':long_colors[10], 
    'GPB2':long_colors[2], # green

    'PAN2':'#a50f15',
    'PAN3':'#a50f15',


    'IRA1':long_colors[0], # blue (obviously)
    'IRA2':long_colors[0], # same color as Ira1
    'double_mutant':'k'
    }




labels = {'FerPerHour':'Fermentation per Hour','ResPerHour':'Respiration per Hour','StaPerHour':'Stationary per Hour'}
lims = {'FerPerHour':[-0.01,0.08],
        'ResPerHour':[-0.04,0.12],
        'StaPerHour':[-0.11,0.02],
        'Fit1D_both2%5%_fitness':[-0.05,1.4],
        'Fit2D_early_fitness':[-0.3,3.5],
        'Fit3D_both2%5%_fitness':[-1.0,3.1],
        'Fit5D_both2%5%_fitness':[-4.2,3.2]}




new_lowcomplexity_dict = {
#     'EVO1D_IRAs':[
# #         'CATAAAAAGACTAATCTTATTAATGC', # 74D-2 (Ira1 mis)
#                   'AATGCAATAATGAAATGATTTGAGGA',
# #                  'TGTACAAATCTTAAGAAGATTACAAG',  #YP 12 (IRA1 non)
#                   'GAAATAAACCACAACGACATTCTAAT',
# #                   'CCTACAAATACTAAGGCTATTCCTAT', #YP 10 (IRA1 non)
#                  ]
    
    'EVO1D_IRA1_MIS':['CATAAAAAGACTAATCTTATTAATGC'],
    
    'EVO1D_IRA1_NON':['CCTACAAATACTAAGGCTATTCCTAT','TGTACAAATCTTAAGAAGATTACAAG',
                     'AATGCAATAATGAAATGATTTGAGGA','GAAATAAACCACAACGACATTCTAAT'],
    
#     'EVO1D_TCG':['CCGCCAATCCCGAACCCCGTTTCGCC','GACAGAAAAGCCAAATGGATTTACCG',
#  'ATCAGAAGTTCGAATCAAATTACGAA','CCAACAAAAGGAAACGTATTTATTGA',
#  'TTAAAAATACAAAAAAAGATTTAAGG','AGAACAAAAACTAAACTCATTCATGG',
#  'ACTTAAAAAGCAAACATGATTATTCA','GTATTAAAATTAAAAATAATTGCACA',
#  'CTAGAAATCTCAAAAACTTTTGGCTG','CAGAAAAGCCATAACGCTATTTGAAA'],
    
    'EVO2D_IRA1_NON':['CCAACAAAACACAAATCTGTTGTGTA'],

    'EVO2D_IRA1_MIS':['TGATCAATCTACAAAAATATTTAATG','CATTGAATCACAAAATAGGTTAGATG'],

    'EVO3D_IRA1_NON':['ATCACAATAACTAAACTGATTCTTCA'],

    'EVO3D_IRA1_MIS':['TATCGAAACCCAAAGAGATTTAATCG'],

    'CYR1':['AGAACAAAAACTAAACTCATTCATGG','GACAGAAAAGCCAAATGGATTTACCG',
           'CTAGAAATCTCAAAAACTTTTGGCTG','ACTTAAAAAGCAAACATGATTATTCA',
            'CAGAAAAGCCATAACGCTATTTGAAA'],
    
    'GPB2':['CCGCCAATCCCGAACCCCGTTTCGCC','GTATTAAAATTAAAAATAATTGCACA',
            'CCAACAAAAGGAAACGTATTTATTGA'],
    
    'TOR1':['ATCAGAAGTTCGAATCAAATTACGAA','TTAAAAATACAAAAAAAGATTTAAGG'],
    

}

new_lowcomplexity_bc_to_entry_dict = {}

for key,items in new_lowcomplexity_dict.items():
    for bc in items:
        new_lowcomplexity_bc_to_entry_dict[bc] = key

new_lowcomplexity_barcodes = [item for key,items in new_lowcomplexity_dict.items() for item in items ]

new_lowcomplexity_by_anc = {}
new_lowcomplexity_by_anc['IRA1_MIS'] = []
for key in ['EVO1D_IRA1_MIS','EVO2D_IRA1_MIS','EVO3D_IRA1_MIS']:
    for bc in new_lowcomplexity_dict[key]:
        new_lowcomplexity_by_anc['IRA1_MIS'].append(bc)
    
new_lowcomplexity_by_anc['IRA1_NON'] = []
for key in ['EVO1D_IRA1_NON','EVO2D_IRA1_NON','EVO3D_IRA1_NON']:
    # print(key,new_lowcomplexity_dict[key])
    for bc in new_lowcomplexity_dict[key]:
        new_lowcomplexity_by_anc['IRA1_NON'].append(bc)
        
        
new_lowcomplexity_by_anc['CYR1'] = new_lowcomplexity_dict['CYR1']
new_lowcomplexity_by_anc['GPB2'] = new_lowcomplexity_dict['GPB2']
new_lowcomplexity_by_anc['TOR1'] = new_lowcomplexity_dict['TOR1']


def flatten(list2d):
    return list(itertools.chain.from_iterable(list2d))

def jitter_point(mean,std=0.15):
    return np.random.normal(mean,std)

# def inverse_variance_mean(means,variances,axis=1):
def inverse_variance_mean(means,standard_devs,axis=1):

    variances = np.square(standard_devs)
    # variances = standard_devs

    weighted_mean = np.nansum(means/variances,axis=axis)/(np.nansum(1/variances,axis=axis))

    weighted_standard_devs = np.sqrt((np.nansum(1/variances,axis=axis))**(-1))

    return weighted_mean, weighted_standard_devs

def standard_error(data):

    return np.std(data)/np.sqrt(len(data))

def centroid(arr):
    
    return np.mean(arr,axis=0)


def combine_hex_values(color_list):

    color_list = [to_hex(color) for color in color_list] # convert to hex
  
    weights = [1/len(color_list) for c in color_list]
    tot_weight = sum(weights)

    red = int(sum([int(color[1:3], 16)*weight for color,weight in zip(color_list,weights)])/tot_weight)
    green = int(sum([int(color[3:5], 16)*weight for color,weight in zip(color_list,weights)])/tot_weight)
    blue = int(sum([int(color[5:7], 16)*weight for color,weight in zip(color_list,weights)])/tot_weight)
    zpad = lambda x: x if len(x)==2 else '0' + x
    
    return '#' + zpad(hex(red)[2:]) + zpad(hex(green)[2:]) + zpad(hex(blue)[2:])


def find_mutation_color(gene):

    if gene in mutation_color_map.keys():
        return mutation_color_map[gene]
    else:
        return 'gray'

def hamming_distance(s1,s2):
    assert len(s1) == len(s2)
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1,s2))


def big_extreme(arr):

    return np.abs(arr).max()


