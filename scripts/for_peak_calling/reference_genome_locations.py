#!/usr/bin/python -tt

"""
"""

import re
import sys
import os
from collections import defaultdict
import glob



species_for_variant_calling = ['Nmel', 
                               'Aaur', 
                               'Apur', 
                               'Avir',
                               'Lalb',
                               'Amel',
                               'Hlig',
                               'Lzep',
                               'Hqua',
                               'Lbal',
                               'Sinv',
                               'Lvie',
                               'Bimp']



new_assembly_dir="/Genomics/kocherlab/lab/data/GenomeReleases/official_release_v3.1"
old_assembly_dir="/Genomics/kocherlab/lab/data/GenomeReleases/official_release_v2.1.1"
other_assembly_dir="/Genomics/kocherlab/bjarnold/STARRseq/data/genomes"

def get_ref_genome_locations(species_for_variant_calling, assembly_dir, suffix, exclude_list=[]):
  # GET REFERENCE GENOME LOCATIONS
  assemblies = defaultdict(str)
  # get list of species, which correspond to the subirectories in the new assembly directory
  subdirs = glob.glob(f'{assembly_dir}/*')
  species = [s.split('/')[-1] for s in subdirs]
  for s in species:
    # s is a 4 character string, but convert the last 3 characters to lower case
    s_lower = s[:-3] + s[-3:].lower()
    if s_lower in species_for_variant_calling and s_lower not in exclude_list:
      ref_fastas_paths = glob.glob(f'{assembly_dir}/{s}/*.{suffix}')
      ref_fastas = [f.split('/')[-1] for f in ref_fastas_paths]
      ref_fasta = [f for f in ref_fastas if '.alt.' not in f][0]

      assemblies[s_lower] = f'{assembly_dir}/{s}/{ref_fasta}'

  for s in assemblies:
    print(s, assemblies[s])
  return assemblies


print('### NEW ASSEMBLIES ###')
new_assemblies = get_ref_genome_locations(species_for_variant_calling, new_assembly_dir, "fasta.gz")
print('### OLD ASSEMBLIES ###')
old_assemblies = get_ref_genome_locations(species_for_variant_calling, old_assembly_dir, "fasta.gz", new_assemblies)
print('### OTHER ASSEMBLIES ###')
other_assemblies = get_ref_genome_locations(species_for_variant_calling, other_assembly_dir, "fna.gz", list(new_assemblies.keys()) + list(old_assemblies.keys()))


sys.exit()

# GET REFERENCE GENOME LOCATIONS
new_assemblies = defaultdict(str)
# get list of species, which correspond to the subirectories in the new assembly directory
subdirs = glob.glob(f'{new_assembly_dir}/*')
species = [s.split('/')[-1] for s in subdirs]
for s in species:
  # s is a 4 character string, but convert the last 3 characters to lower case
  s_lower = s[:-3] + s[-3:].lower()
  if s_lower in species_for_variant_calling:
    ref_fastas_paths = glob.glob(f'{new_assembly_dir}/{s}/*.fasta.gz')
    ref_fastas = [f.split('/')[-1] for f in ref_fastas_paths]
    ref_fasta = [f for f in ref_fastas if '.alt.' not in f][0]

    new_assemblies[s_lower] = f'{new_assembly_dir}/{s}/{ref_fasta}'

for s in new_assemblies:
  print(s, new_assemblies[s])

print('###')
old_assemblies = defaultdict(str)
subdirs = glob.glob(f'{old_assembly_dir}/*')
species = [s.split('/')[-1] for s in subdirs]
for s in species:
  # s is a 4 character string, but convert the last 3 characters to lower case
  s_lower = s[:-3] + s[-3:].lower()
  if s_lower in species_for_variant_calling and s_lower not in new_assemblies:
    ref_fastas_paths = glob.glob(f'{old_assembly_dir}/{s}/*_v2.1.0.fasta.gz')
    ref_fastas = [f.split('/')[-1] for f in ref_fastas_paths]
    assert len(ref_fastas) == 1, print(s)
    ref_fasta = ref_fastas[0]

    old_assemblies[s_lower] = f'{old_assembly_dir}/{s}/{ref_fasta}'

for s in old_assemblies:
  print(s, old_assemblies[s])

print('###')
# GET REFERENCE GENOME LOCATIONS
other_assemblies = defaultdict(str)
# get list of species, which correspond to the subirectories in the new assembly directory
subdirs = glob.glob(f'{other_assembly_dir}/*')
species = [s.split('/')[-1] for s in subdirs]
for s in species:
  # s is a 4 character string, but convert the last 3 characters to lower case
  s_lower = s[:-3] + s[-3:].lower()
  if s_lower in species_for_variant_calling:
    ref_fastas_paths = glob.glob(f'{other_assembly_dir}/{s}/*.fna')
    ref_fastas = [f.split('/')[-1] for f in ref_fastas_paths]
    ref_fasta = [f for f in ref_fastas if '.alt.' not in f][0]

    other_assemblies[s_lower] = f'{other_assembly_dir}/{s}/{ref_fasta}'

for s in other_assemblies:
  print(s, other_assemblies[s])