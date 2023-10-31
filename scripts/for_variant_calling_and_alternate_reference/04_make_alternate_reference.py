#!/usr/bin/python -tt

"""
"""

import re
import sys
import os
from collections import defaultdict

def main():

  # use the reference genomes from snpArcher, since these have all of the extra files we need, such as a reference index, dict, etc
  ref_path = "/Genomics/kocherlab/bjarnold/STARRseq/code/snakemake_snparcher/results"
  vcf_dir = "/Genomics/kocherlab/bjarnold/STARRseq/data/VCFs/alternate_reference/01_filtered_within_sample_allele_frequency"
  new_fasta_dir = "/Genomics/kocherlab/bjarnold/STARRseq/data/alternate_references"

  for species in ['Nmel','Aaur','Apur','Avir']:
    vcf = f'{vcf_dir}/{species}.vcf.gz'
    new_fasta = f'{new_fasta_dir}/{species}.fasta'

    cmd = "gatk FastaAlternateReferenceMaker "
    cmd += f"-R {ref_path}/{species}/data/genome/{species}.fna "
    cmd += f"-O {new_fasta} "
    cmd += f"-V {vcf} "
    os.system(cmd)

if __name__ == '__main__':
  main()
