#!/usr/bin/python -tt

"""
"""

import re
import sys
import os
from collections import defaultdict

def main():

  ref_dir = "/Genomics/kocherlab/bjarnold/STARRseq/data/alternate_references"
  species = ['Nmel','Aaur','Apur','Avir']

  for s in species:
    new_fasta = []
    with open(f'{ref_dir}/{s}.fasta','r') as fasta:
      for line in fasta:
        if line.startswith('>'):
          scaffold = line.strip().split(':')[0]
          scaffold = scaffold.split(' ')[1]
          print(f'>{scaffold}')
          new_fasta.append(f'>{scaffold}')
        else:
          new_fasta.append(line.strip())
          
    with open(f'{ref_dir}/renamed/{s}.fasta','w') as out:
      out.write('\n'.join(new_fasta))

if __name__ == '__main__':
  main()
