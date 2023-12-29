#!/usr/bin/python -tt

"""
"""

import re
import sys
import os
from collections import defaultdict

def main():

  var_call_dir = "/Genomics/kocherlab/bjarnold/STARRseq/code/snakemake_snparcher/results"
  # species = ['Nmel', 'Aaur', 'Apur', 'Avir', 'Amel']
  species = ['Bimp']
  out_dir = "/Genomics/kocherlab/bjarnold/STARRseq/data/VCFs/alternate_reference/00_filtered"
  minDP = 4

  for s in species:
    individuals = [f'--indv {s}_RNA_{biorep}' for biorep in ['F1', 'F2', 'F3']]
    individuals = ' '.join(individuals)
    print(individuals)
    vcf = f'{var_call_dir}/{s}/final_raw.vcf.gz'
    # cmd = f'vcftools --gzvcf {vcf} --remove-filtered-all --minDP {minDP} --max-missing-count 1 --min-alleles 2 --max-alleles 2 {individuals} --recode --stdout | gzip -c > {out_dir}/{s}.vcf.gz'
    cmd = f'vcftools --gzvcf {vcf} --minDP {minDP} --max-missing-count 1 --min-alleles 2 --max-alleles 2 {individuals} --recode --stdout | gzip -c > {out_dir}/{s}.vcf.gz'
    print(cmd)
    os.system(cmd)


if __name__ == '__main__':
  main()
