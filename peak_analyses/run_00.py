#!/usr/bin/python -tt

"""
"""

import re
import sys
import os
from collections import defaultdict

def main():

  species = ['Nmel', 'Apur', 'Avir', 'Aaur', 'Lvie', 'Bimp', 'BimpMET', 'BimpDMSO', 'Hlig', 'Hqua', 'Lzep', 'Lbal', 'Amel']
  # species = ['Amel_prop0.1', 'Amel_prop0.2', 'Amel_prop0.4', 'Amel_prop0.6', 'Amel_prop0.8']
  # species = ['Amel_prop0.05', 'Amel_prop0.01', 'Amel_prop0.001']

  min_biorep_support = [3,2] # [3,2]
  fold_diffs = [0,2] # [0,2]
  for s in species:
    for m in min_biorep_support:
      for f in fold_diffs:
        print("RUNNING SPECIES: ", s, "MIN_BIOREP_SUPPORT: ", m, "FOLD_DIFF: ", f)
        os.system(f"python 00_raw_deduplicated_compare_MACS_peaks.py {s} {m} {f}")

if __name__ == '__main__':
  main()
