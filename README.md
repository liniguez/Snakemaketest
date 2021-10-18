# Snakemake test for Irimia's Lab
Testing Snakemake in CRG cluster

The only requirement for running these scripts is to have conda isntalled. More information about conda and how to use it and install it [here](https://bioconda.github.io/)


## Option 1: Conda environment for local run

```{bash}
conda install mamba -n base
mamba create -n snakemaketest snakemake=6.9.1 bioconductor-tximport=1.20.0 \
      r-readr=2.0.2 bioconductor-genomicfeatures=1.44.0 salmon=1.5.2 future=0.18.2 \
      parallel-fastq-dump=0.6.6 bowtie=1.3.1
#it might take a while!
```

```{bash}
qlogin -l virtual_free=4G -pe smp 8 -q interactive
git clone https://github.com/LIniguez/Snakemaketest.git
cd Snakemaketest
```

```{bash}
snakemake --keep-going --cores 8 --printshellcmds --reason
```

## Option 2: Conda environment for Cluster/Parallel run

Another posibility is to run snakemake in the login node with a screen. (More info about screen [here](https://linuxize.com/post/how-to-use-linux-screen/)) Using this option snakemake will submit all parallelizable jobs to different nodes. For this option you won't need to install all packages at the beginnig and instead it will use conda environments.

```{bash}
conda install mamba -n base
mamba create -n snakemaketest snakemake=6.9.1
git clone https://github.com/LIniguez/Snakemaketest.git
cd Snakemaketest
snakemake -s Snakefile_cluster --profile CRG_profile/
```

*Note: for this option you will need to have a forder in the workinfolder named envs/ with all the environments you will use. Additionally you will be using now a profile where you specify the resurces needed for each rule!*
Option -s stands for a specific snakefile script. 
