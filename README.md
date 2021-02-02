# Snakemake workflow: dna-seq-pipeline

This Snakemake pipeline implements the DNA-seq pipeline for calling SNV/CNV/SV variants.

## Authors

* Jingyu Yang


## Usage


#### Step 1: Obtain a copy of this workflow

1. Create a new github repository using this workflow [as a template](https://help.github.com/en/articles/creating-a-repository-from-a-template).
2. [Clone](https://help.github.com/en/articles/cloning-a-repository) the newly created repository to your local system, into the place where you want to perform the data analysis.

#### Step 2: Configure workflow

Configure the workflow according to your needs via editing the files `config.yaml`.
Set the raw fastq data located directory at the `input_dir` field.
Make sure there only one sample's fastq data located in your raw fastq data path directory.
Eg.  `sample.R1.fq`, `sample.R2.fq`
#### Step 3: Download reference and control data

1. Follow the intruction of README.md in `Ref_data` folder and download required reference files.
2. Follow the intruction of README.md in `control_bam` folder and download/prepare the required control bam file.

#### Step 4: Install snakemake if not installed

Follow the command below to install snakemake, be sure you have installed the [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/).

    conda install -c conda-forge mamba
    mamba create -c conda-forge -c bioconda -n snakemake snakemake	

#### Step 5: Activate the snakmake environment

Activate the snakemake environment using conda:

    conda activate snakemake

#### Step 6: Run DNA-seq pipeline

1. Dry run pipeline to check the process steps and availabilities:

    snakemake -j all -n --use-conda

2. Run pipeline:

    snakemake -j all --use-conda

#### Step 7: Check results:

All the variants called results will be generated at a subfolder with the raw data sample's name in `Output` folder :

1. SNV final results will be generated at `filterd_calls` subfolder.
2. CNV final results will be generated at `CNV` subfolder.
3. gene fusion(SV) final results will be generated at `lumpySV` subfolder.
