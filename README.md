# sexdelim

## installation

### clone this repo

git clone https://github.com/mscharmann/sexdelim

### setup a conda environment

conda create --name delimit_sexregions

conda activate delimit_sexregions

conda install snakemake=5.4 bwa samtools bedtools seqtk vcftools bcftools pixy tabix plink parallel freebayes -y

conda install -c conda-forge -c bioconda -c defaults vcflib

### fix small bug in pixy:
nano $(which pixy)

"OK so it helps to just paste this into the pixy executable script in the conda path:"

import os

os.environ["NUMEXPR_MAX_THREADS"]="272"


### get kmerGO

wget https://github.com/ChnMasterOG/KmerGO/releases/download/v1.5.0/KmerGO_for_linux_x64_cmd.zip

unzip KmerGO_for_linux_x64_cmd.zip

cd KmerGO_for_linux_x64_cmd

chmod +x KmerGO

chmod +x bin/kmc

chmod +x bin/kmc_tools

chmod +x bin/kmc_dump

chmod +x bin/cap3

### get accessory script:
wget https://gist.githubusercontent.com/travc/0c53df2c8eca81c3ebc36616869930ec/raw/eff3032ca7c955ca33bffd8758092e4006949c75/split_ref_by_bai_datasize.py



## run

### modify Snakefile
- adjust the input popmap + the map of samples to fq files + the reference genome file

### local:
snakemake -j24 --keep-going --rerun-incomplete

### on SLURM cluster:
snakemake -j 500 --cluster-config cluster.axiom.json --cluster "sbatch -p {cluster.partition} -t {cluster.time} -c {cluster.CPUs} --mem={cluster.RAM_memory}" --restart-times 3 --keep-going --rerun-incomplete


