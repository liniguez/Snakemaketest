# Snakemake test for Irimia's Lab
Testing Snakemake

The only requirement for running these scripts is to have conda isntalled. More information about conda and how to use it and install it [here](https://bioconda.github.io/)


## Conda environment

```{bash}
conda install mamba -n base
mamba create -n snakemaketest snakemake=6.9.1 bioconductor-tximport=1.20.0 \
      r-readr=2.0.2 bioconductor-genomicfeatures=1.44.0 salmon=1.5.2 future=0.18.2 \
      parallel-fastq-dump=0.6.6 bowtie=1.3.1
#it might take a while!
```

```{bash}
qlogin -l virtual_free=4G -pe smp 8 -q interactive
conda activate snakemaketest
mkdir snakemaketest & cd snakemaketest
wget https://raw.githubusercontent.com/LIniguez/Snakemaketest/main/Snakefile
wget https://raw.githubusercontent.com/LIniguez/Snakemaketest/main/samples.txt
wget https://raw.githubusercontent.com/LIniguez/Snakemaketest/main/configure.yaml
```

```{bash}
snakemake --keep-going --cores 8 --printshellcmds --reason
```
