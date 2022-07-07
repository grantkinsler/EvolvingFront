import sys
import os
import pandas as p
import numpy as np
import matplotlib.pyplot as plt 
import seaborn as sns

out_dir = sys.argv[1]
sample_name = sys.argv[2]
window_size = sys.argv[3]
max_coverage = 500

def moving_average(x, w):
    return np.convolve(x, np.ones(w), 'valid') / w

def unfold_array(depth,segment_length):
    unfolded_array = np.zeros(np.sum(segment_length))
    
    running_index= 0
    for d,l in zip(depth,segment_length):
        unfolded_array[running_index:running_index+l] = d
        running_index += l
        
    return unfolded_array

def window_to_median_ratio(sample_coverage):
    
    merged_covs = [y for x in sample_coverage.values() for y in x]
    median_cov = np.median(merged_covs)
    
    ratio = {}
    
    for chrom,depths in sample_coverage.items():
        ratio[chrom] = depths/median_cov
    
    return ratio

def mean_std_coverage_ratio(ratios,chromosome_names=list(chromosome_map.keys())):
    
#     ratios[ratios.key]
    chrom_lengths = [len(x) for x in list(list(ratios.values())[0].values())]
    
    chrom_depths = {}
    for c,chrom in enumerate(list(chromosome_names[:-1])):
        chrom_depths[chrom] = np.zeros((len(ratios.keys()),chrom_lengths[c]))
    
    s = 0
    for sample,chroms in ratios.items():
        for chrom,depth in chroms.items():
            chrom_depths[chrom][s,:] = depth
        s+=1
    
    means = {}
    stds = {}
    z_scores = {}
    for c,chrom in enumerate(list(chromosome_names[:-1])):
        means[chrom] = np.mean(chrom_depths[chrom],axis=0)
        stds[chrom] = np.std(chrom_depths[chrom],axis=0)
        
        means_tiled = np.tile(means[chrom],chrom_depths[chrom].shape[1]).reshape((chrom_depths[chrom].shape[1],len(means[chrome])))
        stds_tiled = np.tile(stds[chrom],chrom_depths[chrom].shape[1]).reshape((chrom_depths[chrom].shape[1],len(stds[chrome])))

        z_scores[chrom] = scipy.stats.zscore(chrom_depths,axis=0)
        
    return means,stds,z_scores,chrom_depths

