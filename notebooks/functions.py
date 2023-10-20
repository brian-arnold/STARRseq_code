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

def load_genrich_pileup_files(genrich_pileup_files, chr1_test=False):
  pileup_dfs = [] # list of pileup dataframes
  for f in genrich_pileup_files:
      df = pd.read_csv(f, sep="\t", header=1)
      df.columns = ["Chromosome", "Start", "End", "experimental", "control", "-log(p)"]
      # change experimental, control, and -log(p) columns to numeric
      df[['Start', 'End', 'experimental', 'control', '-log(p)']] = df[['Start', 'End', 'experimental', 'control', '-log(p)']].apply(pd.to_numeric)
      df['fold_diff'] = df['experimental']/df['control']
      df['midpoint'] = (df['End'] - df['Start'])/2 + df['Start']
      if chr1_test:
          df = df[df['Chromosome'] == "NMEL_chr_1"]
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

