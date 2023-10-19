#!/usr/bin/python -tt

"""
"""

import re
import sys
import os
from collections import defaultdict
import pandas as pd

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
  
