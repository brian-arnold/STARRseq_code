#!/usr/bin/python -tt

"""
"""

import re
import sys
import os
from collections import defaultdict

def main():

  species='Lvie'
  fastq_dir = f'/Genomics/kocherlab/bjarnold/STARRseq/data/FASTQ/{species}'
  bioreps = ['F1', 'F2', 'F3']
  library_types = ['RNA', 'input'] # ['RNA', 'input']
  proportions = [0.4, 0.1] # [0.8, 0.6, 0.4, 0.2, 0.1]

  for b in bioreps:
    for l in library_types:
      for p in proportions:
        fastq_infile_1 = f'{fastq_dir}/{species}-{b}-{l}-read-1.fastq.gz'
        fastq_infile_2 = fastq_infile_1.replace('read-1', 'read-3') # make this read-4 for Amel

        fastq_outfile_1 = f'{fastq_dir}_prop{p}/{species}-{b}-{l}-read-1.fastq.gz'
        fastq_outfile_2 = fastq_outfile_1.replace('read-1', 'read-3') # make this read-4 for Amel

        cmd1 = f'seqkit sample --proportion {p} --rand-seed 11 --threads 8 {fastq_infile_1} -o {fastq_outfile_1}'
        cmd2 = f'seqkit sample --proportion {p} --rand-seed 11 --threads 8 {fastq_infile_2} -o {fastq_outfile_2}'
        os.system(cmd1)
        os.system(cmd2)
    
if __name__ == '__main__':
  main()
