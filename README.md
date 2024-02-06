# STARRseq_code

Calling peaks from DNA input, RNA output FASTQ files

## Before beginning
- In `STARRseq/code/scripts/`
    - Use `bwa_index.py` to create indices of all the genomes you’ll use for peak calling
    - Use `make_repeatmasker_BED.py` to get BED files from repeat modeler output to potentially exclude repeat regions
    - Use `collect_fastqs.py` to copy fastq files for each species into a single location. This script also renames files that have unnecessary strings

## FASTQ to peaks data
- In `STARRseq/code/scripts/for_peak_calling/`
    - Use `compute_effective_genome_size.py` to compute the number of AGCT (non N’s) for effective genome size, this info will get used in the config file for the snakelike workflow. You can use the `reference_genome_locations.py` script to locate the appropriate reference genome to use
    - After you run the snakelike workflow (below), use the `reorganize_files.py` script to copy peak-calling files for all species to a single directory that you can then copy to your local machine for running analyses there
- In ‘/STARRseq/code/snakemake_peak_calling’
    - Put all output in `STARRseq/data/peak_calling_snakemake_output`
    - Run snakemake with `run.sh` script, making sure all relevant parent directories are mounted so the singularity container has access to necessary files


## Exploring peaks, refining output

- The notebook `peak_analyses/00_raw_deduplicated_compare_MACS_peaks.ipynb` performs a variety of tasks
    - Plots many summary statistics of the peak data
    - Uses PyRanges to overlap peaks across biological replicates, e.g. in order to find peaks detected in 2+ or all 3 biological replicates, computing the mean fold-difference across biological replicates, etc.
    - Saves new peak files, where columns `p1`, `p2`, and `p3` indicate whether the peak was present in each biological replicate (possible values here are 0 and 1). There are also fold-difference and peak coordinate (or summit coordinate) columns for each biological replicate that start with `signalValue_` and `peak_coord_`, respectively.
    - Creates FASTA files that contain the sequences underlying peaks. These are used for two downstream applications:
        - the 'Simple Enrichment Analysis' (SEA) tool to study TF motif enrichment
        - DeepSTARR to predict enhancer activity using a convolutional neural network previously trained on Drosophila melanogaster data
- A similarly named script that ends in .py was created from this notebook to automate this analysis across species. It was created using `jupyter nbconvert --to script [YOUR_NOTEBOOK].ipynb`. This python script was run many times, once for each set of conditions, using `run_00.py`.


## Analysing motif enrichment

- files prefixed with `01_` involve running a snakemake workflow that runs SEA on many datasets: 1.) peaks detected in the the "raw" data that contain duplicates, 2.) peaks detected in the deduplicated data and 3.) sequences extracted from random positions across the genome.
- for each data category listed above, there are 500 sequences that were randomly subsampled from a larger number of peaks. Since variable number of peaks were detected in different species and in the "raw" and deduplicated data, I decided to randomly downsample the data to 500 sequences. If there are more sequences, it effects the number of TF motifs that SEA identifies as significantly enriched.
- `02_raw_deduplicated_compare_SEA_output.ipynb` then summarizes the output of these analyses

## Deep Learning to predict enhancer activity

- files prefixed with `05_` involve running DeepSTARR and analyzing the output. This model is a convolutional neural network that has been trained on Drosophila melanogaster data.
- initial attempts to retrain the model on the bee data had limited success, and this code is not included in this repository. Since the Drosphila data was generally higher quality, we decided to just try the pretrained model. After all, this STARRseq data was the result of bee genomes inserted into Drosophila S2 cell lines... so Drosophila transcription factors are driving enhancer activity, so it makes some sense to use a model that recognizes Drosophila TF motifs.

# Deprecated

Below are instructions for how to use a FASTQ to VCF pipeline to create reference genomes that have been updated with polymorphisms found in the specific bees samples. While this approach is being extra careful, I decided to just use the reference genome and ignore polymorphisms detected in the specific samples, assuming a few nucleotide differences here and there wouldn't change the motifs we detect.

Calling short germline variants, from Fastq to VCF

- In `STARRseq/code/scripts/for_variant_calling_and_alternate_reference/`
    - Use `00_make_snparcher_sample_sheet.py` to help make a sample.csv file. NOTE: this also makes a copy of the reference genome that is not compressed, a requirement for snpArcher. These are stored in `ref_genomes_uncompressed` in `snakemake_snparcher`
    - After running snakelike pipeline (below)
        - Run `01_compute_depth_per_sample.py` to compute the depths for each species as a sanity check
        - Run `02_filter_vcfs.py` to filter VCFs for minimum depth, missing genotypes, and balletic sites
        - Run `03_filter_vcfs_within_sample_allele_frequency.py` to create another VCF file, filtered for sites in which 2+ bores have allele frequencies greater than 50%
        - Run `04_make_alternate_reference.py` to create a new reference sequence using the VCF file from the previous step
        - Run `05_rename_scaffolds.py` to revert scaffold names to the original, since these get renamed from the previous step
- In ‘STARRseq/code/snakemake_snparcher’
    - You changed HaplotypeCaller command so that it ignores duplicates flag in BAM
    - In `config/config.yaml`, keep db_scatter_factor set to 1; scattering DB intervals led to an error in which SNPs were duplicated such that merging them back together caused a GATK error
    - Raw VCF files will be stored in the `results` subdirectory