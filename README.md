# sexdelim

a pipeline to delimit sex-linked regions in chromosome-scale genome assemblies, using re-sequencing data of males and females

starting from a genome .fasta and fastq reads, ouputs 11 statistics in windows along the genome:

- count of alignments for male-specific kmers (16-mers), aligned with 1 mismatch

- count of alignments for female-specific kmers (16-mers), aligned with 1 mismatch

- log2 ratio of total normalised male over total normalised female read depth

- GWAS / count of significantly sex-associated variants

- average linkage disequilibrium (LD, r2) in the whole population of males and females

- nucleotide diversity pi in the males

- nucleotide diversity pi in the females

- absolute sequence divergence (dxy) between males and females

- Fst between males and females

- sequence divergence (%) between X-like and Y-like "phased" candidate gametologs

- sequence divergence (%) between Z-like and W-like "phased" candidate gametologs

## installation

### clone this repo
```
git clone https://github.com/mscharmann/sexdelim
```
### setup a conda environment

## create stepwise
```
conda create --name delimit_sexregions
conda activate delimit_sexregions
conda install snakemake=5.4 bwa samtools bedtools seqtk vcftools bcftools pixy tabix plink parallel freebayes -y
conda install -c conda-forge -c bioconda -c defaults vcflib -y
conda install -c conda-forge r-ggplot2 r-cowplot -y
```
## OR use this YAML:

### modify prefix of installation path in last line of this file, then
```
conda env create --file delimit_sexregions.2021-02-16.yml
```
### fix small bug in pixy:
nano $(which pixy)

"OK so it helps to just paste this into the pixy executable script in the conda path:"

import os

os.environ["NUMEXPR_MAX_THREADS"]="272"


### get kmerGO
```
wget -P scripts https://github.com/ChnMasterOG/KmerGO/releases/download/v1.5.0/KmerGO_for_linux_x64_cmd.zip

cd scripts
unzip KmerGO_for_linux_x64_cmd.zip
rm KmerGO_for_linux_x64_cmd.zip
cd KmerGO_for_linux_x64_cmd
chmod +x KmerGO
chmod +x bin/kmc
chmod +x bin/kmc_tools
chmod +x bin/kmc_dump
chmod +x bin/cap3
cd ../../
```
### get accessory script:
```
wget -P scripts https://gist.githubusercontent.com/travc/0c53df2c8eca81c3ebc36616869930ec/raw/eff3032ca7c955ca33bffd8758092e4006949c75/split_ref_by_bai_datasize.py
```


## run

### modify Snakefile
- adjust the input popmap + the map of samples to fq files + the reference genome file

### local:
```
snakemake -j24 --keep-going --rerun-incomplete
```
### on SLURM cluster:
```
snakemake -j 500 --cluster-config cluster.axiom.json --cluster "sbatch -p {cluster.partition} -t {cluster.time} -c {cluster.CPUs} --mem={cluster.RAM_memory}" --restart-times 3 --keep-going --rerun-incomplete
```
### on LSF cluster:
```
snakemake -j 500 --cluster-config cluster.EULER.json --cluster "bsub -W {cluster.time} -n {cluster.CPUs} -R {cluster.mem_and_span}" --restart-times 3 --keep-going --rerun-incomplete
```
# Simulate test data: an XY system with the sexchroms and one autosome

```
import numpy as np
import gzip

n_autosomes = 1

autosomes = []
for i in range(n_autosomes):
	autosomes.append( "".join( np.random.choice(["A","C","T","G"], size = 600000, replace = True) ) )

```
the sexchroms are diverged by a 2 MB region between 2Mbp to 4 Mbp

- Y-hemizygous region 1 Mbp
- X-hemizygous region 1 Mbp => Y chrom is SHORTER
- X-Y gametolog region: divergence = 10%
- PARs on both sides.
```
Xchrom = "".join( np.random.choice(["A","C","T","G"], size = 600000, replace = True) )

```
construct the Ychrom: PAR-X_Y_gametolog_region-Y_hemizygous

1	100000	PAR1

100001	200000	X_Y_gametolog_region

200001	300000	Y-hemizygous

300001	500000	PAR2


and therefore the X chrom is: with the PAR2 being 4-6 Mbp on the X, but 3-5 Mb on the Y (because the X-hemiz is larger (2x) than the Y-hemiz region)

1	100000	PAR1

100001	200000	X_Y_gametolog_region

200001	400000	X-hemizygous

400001	600000	PAR2

```
X_Y_gametolog_region = list( Xchrom[100000:200000] )
snpsites = np.random.uniform(0, len(X_Y_gametolog_region), size = int(0.1*len(X_Y_gametolog_region)) )
print len(snpsites)
for s in snpsites:
	s = int(s)
	isnuc = X_Y_gametolog_region[s]
	newnuc = np.random.choice([ x for x in ["A","C","T","G"] if not x == isnuc])
	X_Y_gametolog_region[s] = newnuc 


X_Y_gametolog_region = "".join(X_Y_gametolog_region)
Y_hemiz_region = "".join( np.random.choice(["A","C","T","G"], size = 100000, replace = True) ) 

Ychrom = Xchrom[:100000] + X_Y_gametolog_region + Y_hemiz_region + Xchrom[400000:]


len(Ychrom)
len(Xchrom)


with open("fakegenome.FEMALE.fa", "w") as O:
	cnt = 0
	for chr in autosomes:
		cnt += 1
		O.write(">chrom" + str(cnt) + "\n")
		O.write(chr + "\n")
	O.write(">chrom_X" + "\n")
	O.write(Xchrom + "\n")


with open("fakegenome.FEMALE.diploid.fa", "w") as O:
	cnt = 0
	for chr in autosomes:
		cnt += 1
		O.write(">chrom" + str(cnt) + "_1" + "\n")
		O.write(chr + "\n")
		O.write(">chrom" + str(cnt) + "_2" + "\n")
		O.write(chr + "\n")
	O.write(">chrom_X" + "_1" + "\n")
	O.write(Xchrom + "\n")
	O.write(">chrom_X" + "_2" + "\n")
	O.write(Xchrom + "\n")


with open("fakegenome.MALE.fa", "w") as O:
	cnt = 0
	for chr in autosomes:
		cnt += 1
		O.write(">chrom" + str(cnt) + "\n")
		O.write(chr + "\n")
	O.write(">chrom_Y" + "\n")
	O.write(Ychrom + "\n")


with open("fakegenome.MALE.diploid.fa", "w") as O:
	cnt = 0
	for chr in autosomes:
		cnt += 1
		O.write(">chrom" + str(cnt) + "_1" + "\n")
		O.write(chr + "\n")
		O.write(">chrom" + str(cnt) + "_2" + "\n")
		O.write(chr + "\n")
	O.write(">chrom_Y" + "\n")
	O.write(Ychrom + "\n")
	O.write(">chrom_X" + "\n")
	O.write(Xchrom + "\n")

```

now simulate some WGS read data!

use wgsim from samtools package

```
for i in {1..3}; do
wgsim -N 160000 -1 150 -2 150 fakegenome.MALE.diploid.fa sample_M_${i}.1.fastq sample_M_${i}.2.fastq
done


for i in {1..3}; do
wgsim -N 160000 -1 150 -2 150 fakegenome.FEMALE.diploid.fa sample_F_${i}.1.fastq sample_F_${i}.2.fastq
done

gzip *.fastq
```

