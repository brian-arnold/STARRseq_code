#!/usr/bin/python -tt

"""
"""

import re
import sys
import os
from collections import defaultdict

def main():

  species = ['Nmel', 'Apur', 'Avir', 'Aaur', 'Lvie', 'Bimp', 'BimpMET', 'BimpDMSO', 'Hlig', 'Hqua', 'Lzep', 'Lbal', 'Amel']
  # species = ['Dmel']
  for s in species:
    print("RUNNING SPECIES: ", s)
    os.system(f"python 00_raw_deduplicated_compare_MACS_peaks.py {s}")

if __name__ == '__main__':
  main()
