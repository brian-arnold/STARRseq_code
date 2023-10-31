#!/usr/bin/python -tt

"""
"""

import re
import sys
import os
from collections import defaultdict

def main():

  var_call_dir = "/Genomics/kocherlab/bjarnold/STARRseq/code/snakemake_snparcher/results"
  species = ['Nmel', 'Aaur', 'Apur', 'Avir']

  for s in species:
    vcf = f'{var_call_dir}/{s}/final_raw.vcf.gz'
    cmd = f'vcftools --gzvcf {vcf} --depth --out {s}'
    os.system(cmd)


if __name__ == '__main__':
  main()
