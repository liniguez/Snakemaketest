# Snakemake test for Irimia's Lab
Testing Snakemake in CRG cluster

The only requirement for running these scripts is to have conda installed. More information about conda and how to use it and install it [here](https://bioconda.github.io/)


## Create conda environment

```{bash}
conda install mamba -n base
mamba create -n snakemaketest snakemake
```

## Screen session

```{bash}
screen
conda activate snakemaketest
```

## Run workflow

```{bash}
snakemake -s Snakemake_PE --profile CRG_profile/
```


*Note: you will need to have a folder in the working folder named envs/ with all the environments you will use. Additionally you will be using now a profile where you specify the resources needed for each rule!*

Option -s stands for a specific snakefile script. You might need to change between _PE and _SE depending on the sequencing technique. 

**Important:** samples.txt has the meta annotations of the sample, you can add as many columns as needed but it must have:
Run,download_path,SampleName,Sp,outFolder

The config.yaml has also important information about reference genome and annotation used for mapping. 

