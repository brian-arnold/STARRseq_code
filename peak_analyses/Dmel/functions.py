#!/usr/bin/python -tt

"""
"""

import re
import sys
import os
from collections import defaultdict
import pandas as pd
import pyranges as pr
import itertools
import glob

def get_files(dir, regex):
    files = glob.glob(f'{dir}/{regex}')
    files = [f.split('/')[-1] for f in files]
    # files = [f.split('__')[-1] for f in files]
    return sorted([f'{dir}/{f}' for f in files])

def load_genrich_pileup_files(genrich_pileup_files, chr1_test=False):
    pileup_dfs = [] # list of pileup dataframes
    for f in genrich_pileup_files:
        df = pd.read_csv(f, sep="\t", header=1)
        df.columns = ["Chromosome", "Start", "End", "experimental", "control", "-log(p)"]
        # change experimental, control, and -log(p) columns to numeric
        df[['Start', 'End', 'experimental', 'control', '-log(p)']] = df[['Start', 'End', 'experimental', 'control', '-log(p)']].apply(pd.to_numeric)
        df['fold_diff'] = df['experimental']/df['control']
        df['midpoint'] = (df['End'] - df['Start'])/2 + df['Start']
        df['bases'] = df['End'] - df['Start']
        df = df[['Chromosome', 'Start', 'End', 'bases', 'midpoint', 'experimental', 'control', 'fold_diff', '-log(p)']]
        if chr1_test:
            df = df[df['Chromosome'] == "NMEL_chr_1"]
        # For any region/interval in the exclude list, 0's will be reported in experiment and control columns (making fold_diff NA), and
        # NA will be reported in in the -log(p) column
        # exclude NA values in -log(p) column
        df = df[df['-log(p)'].notna()]
        df = df[df['fold_diff'].notna()]
        pileup_dfs.append(df)
    return pileup_dfs

def load_peak_caller_results(files, chr1_test=False):
    df_list = []
    for f in files:
        df = pd.read_csv(f, sep="\t", header=None)
        df.columns = ["Chromosome", "Start", "End", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak"]
        df = df.sort_values(by=["Chromosome", "Start"])
        df["peak_coord"] = df["peak"] + df["Start"]
        if chr1_test:
            df = df[df['Chromosome'] == "NMEL_chr_1"]
        df_list.append(df)
    return df_list

def load_repeat_modeler_intervals(repeat_file):
    lines=[]
    with open(repeat_file) as f:
        lines = f.readlines()
        [l.strip() for l in lines]
    # convert lines to data frame
    repeat_df = pd.DataFrame([l.split('\t') for l in lines])
    # merge overlapping interval into one superinterval
    repeat_pr = pr.PyRanges(chromosomes=repeat_df[0],
                    starts=repeat_df[1],
                    ends=repeat_df[2])
    return repeat_pr.merge()

def filter_by_sig_effect_size(dfs, qval_threshold, effect_size_threshold):
    for i,m in enumerate(dfs):
        print("before filtering:", len(m))
        dfs[i] = m[(m['qValue'] > qval_threshold) & 
                        (m['signalValue'] > effect_size_threshold)]
        print("after filtering:", len(dfs[i]))
    return dfs
  
def get_summit_overlaps_between_method_within_reps(macs_dfs, genrich_dfs, dist_thresh):
    # here summit refers to the apex of the peak
    macs_summit_overlaps = []
    genrich_summit_overlaps = []
    for macs_df,genrich_df in zip(macs_dfs, genrich_dfs):
        m = pr.PyRanges(chromosomes=macs_df.Chromosome,
                        starts=macs_df.peak_coord - dist_thresh,
                        ends=macs_df.peak_coord + dist_thresh)
        g = pr.PyRanges(chromosomes=genrich_df.Chromosome,
                        starts=genrich_df.peak_coord - dist_thresh,
                        ends=genrich_df.peak_coord + dist_thresh)
        macs_summit_overlaps.append( (len(m), len(m.intersect(g))/len(m)) )
        genrich_summit_overlaps.append( (len(g),len(g.intersect(m))/len(g)) )
    return macs_summit_overlaps, genrich_summit_overlaps


def get_summit_overlaps_within_method_between_reps(df_list, bioreps, dist_thresh):
    # here summit refers to the apex of the peak
    results = []
    biorep_combinations = list(itertools.combinations(range(len(bioreps)), 2))
    for rep in biorep_combinations:
        r1, r2 = rep
        pr1 = pr.PyRanges(chromosomes=df_list[r1].Chromosome,
                        starts=df_list[r1].peak_coord - dist_thresh,
                        ends=df_list[r1].peak_coord + dist_thresh)
        pr2 = pr.PyRanges(chromosomes=df_list[r2].Chromosome,
                            starts=df_list[r2].peak_coord - dist_thresh,
                            ends=df_list[r2].peak_coord + dist_thresh)
        o = pr1.intersect(pr2)
        results.append((len(pr1), len(pr2), len(o)))
    return results 

def get_peak_overlaps_between_method_within_reps(macs_dfs, genrich_dfs, frac_overlap):
    # here a peak refers to the entire interval called by the method, as opposed to the summit/apex of the peak
    print("warning: code currently doesn't account for when multiple overlaps occur for a single peak")
    macs_peaks_overlaps = []
    genrich_peaks_overlaps = []
    for macs_df,genrich_df in zip(macs_dfs, genrich_dfs):
        m = pr.PyRanges(chromosomes=macs_df.Chromosome,
                        starts=macs_df.Start,
                        ends=macs_df.End)
        g = pr.PyRanges(chromosomes=genrich_df.Chromosome,
                        starts=genrich_df.Start,
                        ends=genrich_df.End)
        # use coverage to find overlaps, covert to dataframe to filter those with at least frac_overlap overlap
        m_cov = m.coverage(g, overlap_col="C", fraction_col="F").df
        g_cov = g.coverage(m, overlap_col="C", fraction_col="F").df
      
        m_cov = m_cov[(m_cov['C']>=1) & (m_cov['F']>frac_overlap)]
        g_cov = g_cov[(g_cov['C']>=1) & (g_cov['F']>frac_overlap)]

        macs_peaks_overlaps.append( (len(m), len(m_cov)/len(m)) )
        genrich_peaks_overlaps.append( (len(g), len(g_cov)/len(g)) )
    return macs_peaks_overlaps, genrich_peaks_overlaps


def get_peak_overlaps_within_method_between_reps(df_list, bioreps, frac_overlap):
    # here a peak refers to the entire interval called by the method, as opposed to the summit/apex of the peak
    results = []
    biorep_combinations = list(itertools.combinations(range(len(bioreps)), 2))
    for rep in biorep_combinations:
        r1, r2 = rep
        pr1 = pr.PyRanges(chromosomes=df_list[r1].Chromosome,
                        starts=df_list[r1].Start,
                        ends=df_list[r1].End)
        pr2 = pr.PyRanges(chromosomes=df_list[r2].Chromosome,
                            starts=df_list[r2].Start,
                            ends=df_list[r2].End)
        cov = pr1.coverage(pr2, overlap_col="C", fraction_col="F").df
        cov = cov[(cov['C']>=1) & (cov['F']>frac_overlap)]
        results.append((len(pr1), len(pr2), len(cov)))  
    return results 


def overlaps_per_replicate(dfs):
    # use pyranges 'count_overlaps' to get intervals supported by one, two, or all three biological replicates
    p1 = pr.PyRanges(chromosomes=dfs[0].Chromosome,
                    starts=dfs[0].Start,
                    ends=dfs[0].End)
    p2 = pr.PyRanges(chromosomes=dfs[1].Chromosome,
                    starts=dfs[1].Start,
                    ends=dfs[1].End)

    gr = {"p1":p1, "p2":p2}
    overlaps_per_rep = pr.count_overlaps(gr)
    overlaps_per_rep = overlaps_per_rep.df
    overlaps_per_rep['combined'] = overlaps_per_rep['p1'].astype(str) + overlaps_per_rep['p2'].astype(str)
    return overlaps_per_rep