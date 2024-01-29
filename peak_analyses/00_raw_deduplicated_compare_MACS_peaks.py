#!/usr/bin/env python
# coding: utf-8

# This notebook extracts MACS peaks that were produced by the snakemake workflow and compares the raw data (no PCR duplicates removed) with data processed to remove PCR duplicates.
# 
# 
# 

# In[ ]:


import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math
import pyranges as pr
import itertools
from collections import defaultdict
from pyfaidx import Fasta
from Bio.Seq import Seq
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler

import pybedtools as pbt

import sys
# sys.path.append("/Genomics/kocherlab/bjarnold/STARRseq/code/notebooks")
import functions as fn


# In[3]:


species = sys.argv[1]
min_biorep_support = int(sys.argv[2]) # minimum number of bioreps to support a peak
qval_threshold = 0 # 1.3 corresponds to pvalue of 0.05, 1 corresponds to pvalue of 0.1
macs_folddiff_threshold = int(sys.argv[3]) # default 2

# these two variables serve as flags as to whether to run the downstream analyses, extracting sequences to use for SEA or deepSTARR
# code related to this is at the very end of this notebook
downstream_sea_anlaysis = False
downstream_deepstarr_analysis = False

# Number of sequences to subsample from each peak set; using the same number to compare across species
number_of_seqs_to_subsample = 500
number_of_subsamples = 100
# size of interval around peak summit to sample sequences from; this is used for grabbing sequences from the reference genome
interval_size = 249
# since this inteval size is used to graab sequences around a peak, if it's odd, shift the lower and upper bound by 1
if interval_size % 2 == 0:
    interval_size_lower, interval_size_upper = int(interval_size//2), int(interval_size//2)
else:
    # if odd, since we cant divide an interval neatly into an equally sized lwoer and upper half, we'll just make the lower half 1 smaller than the upper half
    interval_size_lower, interval_size_upper = int((interval_size-1)//2), int((interval_size+1)//2)

print("These will be the number of bases on either side of the peak summit to grab sequences from:")
print(f"interval_size_lower: {interval_size_lower}")
print(f"interval_size_upper: {interval_size_upper}")

int(interval_size/2)

base_dir = "/Genomics/kocherlab/bjarnold/STARRseq/data"
input_data_dir = f"{base_dir}/peak_calling_snakemake_output/{species}"
# output directory for final, deduplicated peaks
output_dedup_peak_dir = f"{base_dir}/peak_calling_final/min_biorep_support_{min_biorep_support}/qval_{qval_threshold}_folddiff_{macs_folddiff_threshold}/{species}"
# output directory for downstream SEA (simple enrichment analysis) analyses
output_data_sea_dir = f"{base_dir}/meme_suite/min_biorep_support_{min_biorep_support}/qval_{qval_threshold}_folddiff_{macs_folddiff_threshold}/fastas/{species}"
# output directory for downstream deepSTARR analyses
output_data_deepstarr_dir = f"{base_dir}/deepSTARR_input/min_biorep_support_{min_biorep_support}/qval_{qval_threshold}_folddiff_{macs_folddiff_threshold}/{species}"
# output directory for overlapping peaks with ATAC data, for the species that it exists for
output_data_atac_dir = f"{base_dir}/ATAC_MACS_peaks_overlap/MACS_peaks/min_biorep_support_{min_biorep_support}/qval_{qval_threshold}_folddiff_{macs_folddiff_threshold}/{species}"
output_plot_dir = f"{base_dir}/plots_raw_vs_dedup/{species}"


bioreps = ['F1', 'F2', 'F3']
chr1_test = False
min_peak_interval_size = 100 # minimum peak interval size to be considered
frac_overlap = 0.5 # fraction of peak overlap to be considered overlapping

import os
os.system(f"mkdir -p {output_dedup_peak_dir}")
os.system(f"mkdir -p {output_data_sea_dir}")
os.system(f"mkdir -p {output_data_deepstarr_dir}")
os.system(f"mkdir -p {output_data_atac_dir}")
os.system(f"mkdir -p {output_plot_dir}")



# ### Load MACS2 output files, convert to pyranges objects for subsequent overlapping

# In[5]:


# get MACS peak files
macs_raw_files = fn.get_files(f'{input_data_dir}/MACS2', "*_peaks.narrowPeak")
macs_dedup_files = fn.get_files(f'{input_data_dir}/MACS2_dedup', "*_peaks.narrowPeak")
# print(macs_raw_files)
# print(macs_dedup_files)

# load MACS results
macs_raw_dfs = fn.load_peak_caller_results(macs_raw_files, chr1_test) # list of macs dataframes
macs_dedup_dfs = fn.load_peak_caller_results(macs_dedup_files, chr1_test) # list of macs dataframes

# filter MACS results by qval and fold diff
macs_raw_dfs = fn.filter_by_sig_effect_size(macs_raw_dfs, qval_threshold, macs_folddiff_threshold)
macs_dedup_dfs = fn.filter_by_sig_effect_size(macs_dedup_dfs, qval_threshold, macs_folddiff_threshold)

# convert to pyrange objects
macs_raw_prs = [pr.PyRanges(m) for m in macs_raw_dfs]
macs_dedup_prs = [pr.PyRanges(m) for m in macs_dedup_dfs]

print([len(df) for df in macs_raw_dfs])
print([len(df) for df in macs_dedup_dfs])


# load repeatmodeler intervals
repeat_file = f'{base_dir}/repeat_modeler_intervals/{species}_repeat_elements.bed'
# can't read with pandas because variable number of columns; some rows marked with an asterisk if the repeat overlaps with another higher-scoring repeat
# load by reading in all lines instead, then importing into table
# then use pyranges to merge overlapping intervals
# repeat_pr = fn.load_repeat_modeler_intervals(repeat_file)




# ### Load reference genome, which we will use to get sequences within peaks

# In[6]:


# open reference genome ref_genome, count the length of each chromosome, and store in a dictionary
# ref_genome = f'{base_dir}/alternate_references/renamed/{species}.fasta'
ref_genome_locs = '/Genomics/kocherlab/bjarnold/STARRseq/data/reference_genome_locations.txt'
ref_genome = None
with open(ref_genome_locs, 'r') as f:
    for line in f:
        line = line.strip().split()
        if line[0] == species[0:4]:
            ref_genome = line[1]
if ref_genome is None:
    raise Exception(f"Could not find reference genome for {species}")
fasta_sequences = Fasta(ref_genome)

# get scaffold lengths, which we later use to make sure kmers don't overhang edges of scaffolds
scaff_lengths = defaultdict(int)
# Iterate over the sequences and print their names and lengths
for sequence_name in fasta_sequences.keys():
    scaff_lengths[sequence_name] = len(fasta_sequences[sequence_name])


# ## Plot peak interval lengths

# In[7]:


# note for Genrich, peaks thresholded by q value selects for wider peaks with larger q values
fig, axs = plt.subplots(1, 1, figsize=(10, 4), sharex=True, sharey=True)
for df in macs_raw_dfs:
    sns.kdeplot(df['End']-df['Start'], color="blue", alpha=0.5, label="MACS")
    axs.axvline(x=np.mean(df['End']-df['Start']), color="blue", linestyle='--')
for df in macs_dedup_dfs:
    sns.kdeplot(df['End']-df['Start'], color="red", alpha=0.5, label="MACS_dedup")
    axs.axvline(x=np.mean(df['End']-df['Start']), color="red", linestyle='--')
axs.set_xlim(0,3000)
axs.set_xlabel("Peak width (bp)")
axs.set_ylabel("Density")
axs.set_title(f"{species} peak widths")
axs.grid(axis='x', alpha=0.5)  
# add legend
handles, labels = axs.get_legend_handles_labels()
axs.legend(handles, labels)
sns.despine()

plt.savefig(f"{output_plot_dir}/peak_widths.png", dpi=300, bbox_inches='tight')


# ## Plot fold-differences (above background) for raw and deduplicated MACS peaks

# In[8]:


fig, axs = plt.subplots(1, 3, figsize=(15, 4))

for i,(dedup,raw) in enumerate(zip(macs_dedup_prs, macs_raw_prs)):
    tmp = dedup.join(raw, suffix="_raw", how="left", preserve_order=False)
    print("original df length:", len(dedup))
    print("df length after join:", len(tmp))
    sns.scatterplot(data=tmp.df, x="signalValue", y="signalValue_raw", alpha=0.1, ax=axs[i])   
    axs[i].set_xlim(0,7)
    axs[i].set_ylim(0,100)

    axs[i].set_xlabel("deduplicated fold diff")
    axs[i].set_ylabel("raw fold diff")

    # add correlation coefficient to plot
    axs[i].text(0.05, 0.95, f"r={np.corrcoef(tmp.df['signalValue'], tmp.df['signalValue_raw'])[0,1]:.2f}", transform=axs[i].transAxes)
    print(np.corrcoef(tmp.df["signalValue"], tmp.df["signalValue_raw"]))
sns.despine()
# set title
fig.suptitle(f"{species} fold diff correlation between raw and deduplicated peaks")

plt.savefig(f"{output_plot_dir}/fold_diff_correlation.png", dpi=300, bbox_inches='tight')


# ### Overlap peak intervals across biological replicates, for raw and deduplicated data sets

# In[9]:


# look at concordance of peak calling between methods, witihin bioreps
macs_raw_overlaps, macs_dedup_overlaps = fn.get_peak_overlaps_between_method_within_reps(macs_raw_dfs, macs_dedup_dfs, frac_overlap)

print(macs_raw_overlaps)
print(macs_dedup_overlaps)
for i,m in enumerate(macs_raw_overlaps):
    print(f"MACS2 raw: {m[0]} peaks, {m[1]*100:.2f}% overlap with deduplicated peaks")
for i,m in enumerate(macs_dedup_overlaps):
    print(f"MACS2 deduplicated: {m[0]} peaks, {m[1]*100:.2f}% overlap with raw peaks")

# fraction of overlap can seem odd if one caller has smaller peaks that are contained within a larger peak from the other caller


# In[10]:


# look at overlaps within methods between bioreps
macs_raw_overlaps = fn.get_peak_overlaps_within_method_between_reps(macs_raw_dfs, bioreps, frac_overlap)
macs_dedup_overlaps = fn.get_peak_overlaps_within_method_between_reps(macs_dedup_dfs, bioreps, frac_overlap)

print(macs_raw_overlaps)
print(macs_dedup_overlaps)


# In[11]:


biorep_comb = list(itertools.combinations(range(len(bioreps)), 2))
fig, axs = plt.subplots(2, 3, figsize=(18, 8), sharex=False, sharey=False)
for i,comp in enumerate(macs_raw_overlaps):
    rep = biorep_comb[i]
    rep = [int(r)+1 for r in rep]
    l = [f"F{rep[0]}", f"F{rep[1]}", "overlap"]
    axs[0][i].bar(range(len(comp)), comp)
    axs[0][i].set_title(f"MACS RAW F{rep}")
    axs[0][i].set_ylabel("Proportion of peaks overlapping")
    axs[0][i].set_xticks(range(len(l)))
    axs[0][i].set_xticklabels(l, rotation=0)
    axs[0][i].grid(axis='y', alpha=0.5)    
for i,comp in enumerate(macs_dedup_overlaps):
    rep = biorep_comb[i]
    rep = [int(r)+1 for r in rep]
    l = [f"F{rep[0]}", f"F{rep[1]}", "overlap"]
    axs[1][i].bar(range(len(comp)), comp)
    axs[1][i].set_title(f"MACS DEDUP F{rep}")
    axs[1][i].set_ylabel("Proportion of peaks overlapping")
    axs[1][i].set_xticks(range(len(l)))
    axs[1][i].set_xticklabels(l, rotation=0)
    axs[1][i].grid(axis='y', alpha=0.5)    
sns.despine()
# add title
fig.suptitle(f"{species} peak overlap between bioreps")

plt.savefig(f"{output_plot_dir}/peak_overlap_within_method_between_reps.png", dpi=300, bbox_inches='tight')


# ### Categorize peaks by support across bioreps
# For each peak, record how many biological replicates it was detected in

# In[12]:


# convert each df to a pyranges object, count overlaps among biological replicates

macs_raw_overlaps_per_rep = fn.overlaps_per_replicate(macs_raw_dfs)
macs_dedup_overlaps_per_rep = fn.overlaps_per_replicate(macs_dedup_dfs)
# it looks like count_overlaps is counting overlaps b/t all instervals
# if all three replicates have an overlapping interval, pattern could look like this:
#       Chromosome	Start	End	    p1	p2	p3
# 7	    NMEL_chr_1	9079	9297	0	0	0
# 8	    NMEL_chr_1	9297	9309	1	0	0
# 9	    NMEL_chr_1	9309	9649	1	1	0
# 10	NMEL_chr_1	9649	10150	1	1	1
# 11	NMEL_chr_1	10150	10843	0	1	1
# 12	NMEL_chr_1	10843	10914	0	0	1
# 13	NMEL_chr_1	10914	16276	0	0	0
# where index 10 is interval in which all 3 peaks are overlapping, but since start/end points of the peaksare different,
# we can see that replicate 1 has the minimum starting point, and replicate 3 has the maximum ending point

macs_raw_overlaps_per_rep.head()


# In[13]:


print(macs_raw_overlaps_per_rep.combined.unique())
print(macs_dedup_overlaps_per_rep.combined.unique())


# In[14]:


def filter_peaks(overlaps_df, min_interval_size, min_biorep_support, method):
    df = overlaps_df[(overlaps_df['End']-overlaps_df['Start'] > min_interval_size)]
    df.reset_index(inplace=True, drop=True)
    if min_biorep_support == 3:
        df = df[(df['combined'] == '111')]
    elif min_biorep_support == 2:
        df = df[(df['combined'] == '110') | (df['combined'] == '101') | (df['combined'] == '011') | (df['combined'] == '111')]                        
    else:
        raise Exception("min_biorep_support must be 2 or 3")
    df['method'] = method
    return df
    
macs_raw_peaks_df = filter_peaks(macs_raw_overlaps_per_rep, min_peak_interval_size, min_biorep_support, 'macs_raw')
macs_dedup_peaks_df = filter_peaks(macs_dedup_overlaps_per_rep, min_peak_interval_size, min_biorep_support, 'macs_dedup')

print(len(macs_raw_peaks_df), len(macs_dedup_peaks_df))

macs_raw_peaks_df.reset_index(inplace=True, drop=True)
macs_dedup_peaks_df.reset_index(inplace=True, drop=True)


def get_signal_val_and_peak_coord_for_trips_df(trips_df, prs):
    tmp_pr = pr.PyRanges(trips_df)
    for i in range(len(prs)):
        tmp_join = tmp_pr.join(prs[i], suffix="_raw", how="left", preserve_order=True)
        tmp_join = tmp_join.df.reset_index(drop=True, inplace=False)
        # concatenate signalValue as column to macs_raw_peaks_df
        assert len(tmp_join['signalValue']) == len(trips_df)
        trips_df[f'signalValue_{i}'] = tmp_join['signalValue']
        trips_df[f'peak_coord_{i}'] = tmp_join['peak_coord'] # this is the peak summit reported by MACs, which can be outside the peak interval in which 2+ bioreps overlap
    return trips_df

macs_raw_peaks_df = get_signal_val_and_peak_coord_for_trips_df(macs_raw_peaks_df, macs_raw_prs)
macs_dedup_peaks_df = get_signal_val_and_peak_coord_for_trips_df(macs_dedup_peaks_df, macs_dedup_prs)
print(len(macs_raw_peaks_df), len(macs_dedup_peaks_df))


# In[15]:


macs_dedup_peaks_df


# ### Save peak intervals to file

# In[19]:


# select everything but the 'combined' column in macs_dedup_peaks_df
macs_dedup_peaks_df.drop(columns=['combined', 'method']).to_csv(f"{output_dedup_peak_dir}/macs_dedup_peaks.tsv", sep='\t', index=False)


# In[16]:


def remove_bookended_intervals_with_min_support_bioreps_2(df_input):
    # peaks detected in all 3 biological replicates are typically NOT bookended
    # however, peaks detected in only 2 biological replicates are sometimes bookended with peaks detected in all three biological replicates
    # this is because near a given enhancer, support across all three replicates will vary, such that intervals supported by 3 replicates will sometimes be flanked by intervals supported by only 2 replicates
    # since the signal fell below some threshold for one of the replicates
    bookended_indices = []
    # work by chromosome
    for chrom in df_input.Chromosome.unique():
        df = df_input[df_input['Chromosome'] == chrom]
        df.reset_index(inplace=True, drop=False) # don't drop index because we need it to remove rows later
        if len(df) == 1:
            continue
        for i,r in df.iterrows():
            if r['combined'] == '111':
                # only remove intervals with support from less than 3 bioreps
                continue
            if i==0:
                # if i == 0, see if bookended with i+1
                if r['End'] == df.iloc[i+1]['Start']:
                    bookended_indices.append(r['index'])
            elif i==len(df)-1:
                # if i == len(df)-1, see if bookended with i-1
                if r['Start'] == df.iloc[i-1]['End']:
                    bookended_indices.append(r['index'])
            else:
                # else compare to both flanking indices
                if (r['End'] == df.iloc[i+1]['Start']) or (r['Start'] == df.iloc[i-1]['End']):
                    bookended_indices.append(r['index'])
    print(len(bookended_indices))

    return df_input[~df_input.index.isin(bookended_indices)].reset_index(drop=True)


print("MACS raw before removing bookends:", len(macs_raw_peaks_df))
macs_raw_peaks_df = remove_bookended_intervals_with_min_support_bioreps_2(macs_raw_peaks_df)
print("MACS raw after removing bookends:", len(macs_raw_peaks_df))

print("MACS dedup before removing bookends:", len(macs_dedup_peaks_df))
macs_dedup_peaks_df = remove_bookended_intervals_with_min_support_bioreps_2(macs_dedup_peaks_df)
print("MACS dedup after removing bookends:", len(macs_dedup_peaks_df))


# In[17]:


# make pairwise comparisons of signal values columns
fig, axs = plt.subplots(2, 3, figsize=(10, 8), sharex=False, sharey=False)
#get all combinations of integers between 0 and 2
for i, comb in enumerate(list(itertools.combinations(range(3), 2))):
    sns.scatterplot(data=macs_raw_peaks_df, x=f"signalValue_{comb[0]}", y=f"signalValue_{comb[1]}", alpha=0.1, ax=axs[0][i])
    axs[0][i].set_xlim(0,30)
    axs[0][i].set_ylim(0,30)
    # print pearson correlation coefficient
    axs[0][i].text(0.05, 0.95, f"r={np.corrcoef(macs_raw_peaks_df[f'signalValue_{comb[0]}'], macs_raw_peaks_df[f'signalValue_{comb[1]}'])[0,1]:.2f}", transform=axs[0][i].transAxes)


    sns.scatterplot(data=macs_dedup_peaks_df, x=f"signalValue_{comb[0]}", y=f"signalValue_{comb[1]}", alpha=0.1, ax=axs[1][i])
    axs[1][i].set_xlim(0,8)
    axs[1][i].set_ylim(0,8)
    # print pearson correlation coefficient
    axs[1][i].text(0.05, 0.95, f"r={np.corrcoef(macs_dedup_peaks_df[f'signalValue_{comb[0]}'], macs_dedup_peaks_df[f'signalValue_{comb[1]}'])[0,1]:.2f}", transform=axs[1][i].transAxes)
sns.despine()

# add title
fig.suptitle(f"{species} triplet peak signal correlation; row 1: raw peaks, row 2: deduplicated peaks")
plt.savefig(f"{output_plot_dir}/triplet_peak_signal_correlation_across_reps.png", dpi=300, bbox_inches='tight')


# In[18]:


macs_dedup_peaks_df


# ### compute mean signal value, mean peak coordinate

# In[158]:


def compute_mean_value_across_reps(df, col_prefix):
    # we need a more complicated function here to accomodate for peaks that are only supported by 2 bioreps
    mean_values = []
    for i,r in df.iterrows():
        # only take values from bioreps that contained that peak interval
        bioreps_to_avg = [i for i in range(3) if r['combined'][i] == '1']
        mean_values.append( np.mean([r[f'{col_prefix}_{i}'] for i in bioreps_to_avg]) )
    return mean_values

# compute the mean signal value across replicates
macs_raw_peaks_df['mean_signalValue'] = compute_mean_value_across_reps(macs_raw_peaks_df, 'signalValue')
macs_dedup_peaks_df['mean_signalValue'] = compute_mean_value_across_reps(macs_dedup_peaks_df, 'signalValue')

# get midpoint peak_coord and convert to integer
macs_raw_peaks_df['mid_peak_coord'] = (macs_raw_peaks_df['Start'] + macs_raw_peaks_df['End'])//2
macs_dedup_peaks_df['mid_peak_coord'] = (macs_dedup_peaks_df['Start'] + macs_dedup_peaks_df['End'])//2

# get summit peak_coord
macs_raw_peaks_df['mean_summit_coord'] = compute_mean_value_across_reps(macs_raw_peaks_df, 'peak_coord')
macs_dedup_peaks_df['mean_summit_coord'] = compute_mean_value_across_reps(macs_dedup_peaks_df, 'peak_coord')

# convert 'mean_summit_coord' to integer
macs_raw_peaks_df['mean_summit_coord'] = macs_raw_peaks_df['mean_summit_coord'].astype(int)
macs_dedup_peaks_df['mean_summit_coord'] = macs_dedup_peaks_df['mean_summit_coord'].astype(int)

# # compute the mean signal value across replicates
# macs_raw_peaks_df['mean_signalValue'] = macs_raw_peaks_df[[f'signalValue_{i}' for i in range(3)]].mean(axis=1)
# macs_dedup_peaks_df['mean_signalValue'] = macs_dedup_peaks_df[[f'signalValue_{i}' for i in range(3)]].mean(axis=1)

# # get mean peak_coord and convert to integer
# macs_raw_peaks_df['mean_peak_coord'] = macs_raw_peaks_df[[f'peak_coord_{i}' for i in range(3)]].mean(axis=1).astype(int)
# macs_dedup_peaks_df['mean_peak_coord'] = macs_dedup_peaks_df[[f'peak_coord_{i}' for i in range(3)]].mean(axis=1).astype(int)



# ### get interval coordinates around peak mid-point coordinate, these intervals will be used to look for enrichment in TF motifs using SEA

# In[159]:


# extend mid peak coordinate to match interval size
macs_raw_peaks_df['mid_Start'] = macs_raw_peaks_df['mid_peak_coord'] - (interval_size_lower)
macs_raw_peaks_df['mid_End'] = macs_raw_peaks_df['mid_peak_coord'] + (interval_size_upper)

macs_dedup_peaks_df['mid_Start'] = macs_dedup_peaks_df['mid_peak_coord'] - (interval_size_lower)
macs_dedup_peaks_df['mid_End'] = macs_dedup_peaks_df['mid_peak_coord'] + (interval_size_upper)



# ### get interval coordinates around peak summit

# In[160]:


macs_raw_peaks_df['peak_Start'] = macs_raw_peaks_df['mean_summit_coord'] - (interval_size_lower)
macs_raw_peaks_df['peak_End'] = macs_raw_peaks_df['mean_summit_coord'] + (interval_size_upper)

macs_dedup_peaks_df['peak_Start'] = macs_dedup_peaks_df['mean_summit_coord'] - (interval_size_lower)
macs_dedup_peaks_df['peak_End'] = macs_dedup_peaks_df['mean_summit_coord'] + (interval_size_upper)


# In[161]:


macs_dedup_peaks_df


# ### When we created intervals around some coordinate, there's a small chance we went beyond the limits of the chromosome/scaffold

# In[162]:


# get rid of intervals that go beyond boundaries of scaffolds
def get_indices_of_overhangs(df, scaff_lengths, coord_type):
    # get indices of intervals that overhang scaffolds
    overhang_indices = []
    for i,r in df.iterrows():
        if (r[f'{coord_type}_Start'] < 0) or (r[f'{coord_type}_End'] > scaff_lengths[r['Chromosome']]):
            overhang_indices.append(i)
    return overhang_indices

macs_raw_overhang_indices = get_indices_of_overhangs(macs_raw_peaks_df, scaff_lengths, 'peak') 
macs_dedup_overhang_indices = get_indices_of_overhangs(macs_dedup_peaks_df, scaff_lengths, 'peak')

print("macs_raw_overhang_indices overhangs: ", len(macs_raw_overhang_indices))
print("macs_dedup_overhang_indices overhangs: ", len(macs_dedup_overhang_indices))

# get all indices of macs_raw_peaks_df except for those in macs_raw_overhang_indices
macs_raw_peaks_df = macs_raw_peaks_df[~macs_raw_peaks_df.index.isin(macs_raw_overhang_indices)].reset_index(drop=True)

macs_dedup_peaks_df = macs_dedup_peaks_df[~macs_dedup_peaks_df.index.isin(macs_dedup_overhang_indices)].reset_index(drop=True)
 


# ### Record the number of MACS peaks for downstream analyses

# In[163]:


with open(f"./MACS2_summaries/min_biorep_support_{min_biorep_support}/{species}_MACS2_peaks_qval_{qval_threshold}_folddiff_{macs_folddiff_threshold}.txt",'w') as o:
    o.write(" ".join(["species", "raw", "dedup", "\n"]))
    o.write(" ".join([species, str(len(macs_raw_peaks_df)), str(len(macs_dedup_peaks_df)), "\n"]))


# In[164]:


macs_raw_peaks_df


# ## Save intervals as BED file, load into pybedtools object to extract underlying sequences from FASTA file, since this function is bugged in PyRanges
# 
# ### For 'SIMPLE ENRICHMENT ANALYSIS" (SEA): Since raw data likely has many more peaks compared to deduplicated data, downsample raw data to as many observations in deduplicated data, since you may find more TFBS enrichment from sampling random regions of the genome

# In[165]:


# load in reference genome as pybedtools object
ref_genome_pbt = pbt.example_filename(ref_genome)

if downstream_sea_anlaysis:

    # randomly sample rows from macs_raw_peaks_df
    for i in range(number_of_subsamples):
        macs_raw_trips_subsamp_df = macs_raw_peaks_df.sample(n=number_of_seqs_to_subsample, replace=True).reset_index(drop=True) 
        macs_raw_trips_subsamp_df[['Chromosome', 'mid_Start', 'mid_End']].to_csv(f'{output_data_sea_dir}/MACS2_triplet_peaks_raw_rep{i}.bed', sep='\t', header=False, index=False)
        macs_raw_trips_subsamp_pbt = pbt.BedTool(f'{output_data_sea_dir}/MACS2_triplet_peaks_raw_rep{i}.bed')
        macs_raw_trips_subsamp_pbt.sequence(fi=ref_genome_pbt, fo=f"{output_data_sea_dir}/MACS2_triplet_peaks_raw_rep{i}.fa")

    # randomly sample rows from macs_dedup_peaks_df
    for i in range(number_of_subsamples):
        macs_dedup_trips_subsamp_df = macs_dedup_peaks_df.sample(n=number_of_seqs_to_subsample, replace=True).reset_index(drop=True) 
        macs_dedup_trips_subsamp_df[['Chromosome', 'mid_Start', 'mid_End']].to_csv(f'{output_data_sea_dir}/MACS2_triplet_peaks_dedup_rep{i}.bed', sep='\t', header=False, index=False)
        macs_dedup_trips_subsamp_pbt = pbt.BedTool(f'{output_data_sea_dir}/MACS2_triplet_peaks_dedup_rep{i}.bed')
        macs_dedup_trips_subsamp_pbt.sequence(fi=ref_genome_pbt, fo=f"{output_data_sea_dir}/MACS2_triplet_peaks_dedup_rep{i}.fa")

    for i in range(number_of_subsamples):
        x = pbt.BedTool()
        y = x.random(l=interval_size, n=number_of_seqs_to_subsample, g=f"{ref_genome}.fai")
        y.sequence(fi=ref_genome_pbt, fo=f"{output_data_sea_dir}/MACS2_triplet_peaks_random_rep{i}.fa")


# print(open(tmp.seqfn).read())


# In[166]:


if downstream_deepstarr_analysis:

    macs_raw_peaks_df[['Chromosome', 'peak_Start', 'peak_End', 'mean_signalValue']].to_csv(f'{output_data_deepstarr_dir}/MACS2_peaks_raw.bed', sep='\t', header=False, index=False)
    macs_raw_peaks_pbt = pbt.BedTool(f'{output_data_deepstarr_dir}/MACS2_peaks_raw.bed')
    macs_raw_peaks_pbt.sequence(fi=ref_genome_pbt, fo=f"{output_data_deepstarr_dir}/MACS2_peaks_raw.fa")

    macs_dedup_peaks_df[['Chromosome', 'peak_Start', 'peak_End', 'mean_signalValue']].to_csv(f'{output_data_deepstarr_dir}/MACS2_peaks_dedup.bed', sep='\t', header=False, index=False)
    macs_dedup_peaks_pbt = pbt.BedTool(f'{output_data_deepstarr_dir}/MACS2_peaks_dedup.bed')
    macs_dedup_peaks_pbt.sequence(fi=ref_genome_pbt, fo=f"{output_data_deepstarr_dir}/MACS2_peaks_dedup.fa")

    x = pbt.BedTool()
    y = x.random(l=interval_size, n=10000, g=f"{ref_genome}.fai", seed=1)
    y.sequence(fi=ref_genome_pbt, fo=f"{output_data_deepstarr_dir}/random_seqs.fa")


# # Output peak interval files to overlap with ATACseq data (which only exists for a few species)

# In[168]:


macs_raw_peaks_df.to_csv(f'{output_data_atac_dir}/MACS2_peaks_raw.tsv', sep='\t', header=True, index=False)
macs_dedup_peaks_df.to_csv(f'{output_data_atac_dir}/MACS2_peaks_dedup.tsv', sep='\t', header=True, index=False)

