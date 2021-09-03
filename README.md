# sexdelim

a pipeline to delimit sex-linked regions in chromosome-scale genome assemblies, using re-sequencing data of males and females

starting from a genome .fasta and fastq reads, generates 15 statistics in windows along the genome:

- log2 ratio of total normalised male over total normalised female read depth

- count of alignments for female-specific kmers (16-mers), aligned with 1 mismatch

- count of alignments for male-specific kmers (16-mers), aligned with 1 mismatch

- nucleotide diversity pi in the females

- nucleotide diversity pi in the males

- Fst between males and females

- absolute sequence divergence (dxy) between males and females

- average linkage disequilibrium (LD, r2) in the whole population of males and females

- GWAS / count of significantly sex-associated variants

- heterozygosity in males

- heterozygosity in females

- sequence divergence (%) between Z-like and W-like "phased" candidate gametologs

- sequence divergence (%) between X-like and Y-like "phased" candidate gametologs

- net divergence of the sexes = absolute sequence divergence (dxy) - pi in the females

- net divergence of the sexes = absolute sequence divergence (dxy) - pi in the males


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
conda install snakemake=5.4 scipy bwa samtools bedtools seqtk vcftools bcftools tabix plink parallel freebayes -y
conda install -c conda-forge -c bioconda -c defaults vcflib -y
conda install -c conda-forge r-ggplot2 r-cowplot -y
```
## OR use this YAML:

### modify prefix of installation path in last line of this file, then
```
conda env create --file delimit_sexregions.2021-06-30.yml
```

### get kmerGO
```
wget -P scripts https://github.com/ChnMasterOG/KmerGO/releases/download/v1.5.0/KmerGO_for_linux_x64_cmd.zip

cd scripts
unzip KmerGO_for_linux_x64_cmd.zip
rm KmerGO_for_linux_x64_cmd.zip
cd KmerGO_for_linux_x64_cmd
mv KmerGO_for_cmd KmerGO
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

### config.yaml: provide meta-data and set parameters
This file should be mostly self-explanatory and also contains some comments. Further important points are:
- popmap = 2-column tab-sep file listing samples with assignment as 1 (M) or 2 (F)
- samples_reads_map = a table listing the samples and "mapping" to their corresponding fastq files. Both absolute paths and paths relative to the working directory should be tolerated. Both single and paired-end data can be handled, also mixed. CAUTION: currently, each sample may have at most one file (one pair of read files). If the reads are split into several files (pairs of files), these must be concatenated before start.
- genome fasta file path
- regions_for_plot_bed: a BED file with a region of special interest, e.g. the sex chromosome. The pipeline produces statistics for all genomic windows and then extracts a subset (region). Plots will be made for all genomic windows, and for the subset region. This can be useful if the assembly contains many contigs that would make a messy plot.  
- windowsize = length of intervals (windows) in which to split all statistics and plots. Most useful values are between 1000 and 100000.
- VCF_MIN_DEPTH: the default here is 6 (retain all genotypes with at least 6 reads covering the site). I would reduce to as low as 3, but not less.	


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
## post-run
- Inspect the PDF files and .txt tables with statistics in the directory results_processed.
- results_processed are specific for the given windowsize. If a different windowsize is desired, intermediate results in the directory results_raw can be used to generate these relatively quickly, without going back to the read data and variant-calling. Just move the directory results_raw, change the windowsize in config.yaml, and start the pipeline again exactly as before. I usually do this several times for window sizes 1 kb, 10 kb, 25 kb, 50 kb, 100 kb.
- If a different subset region for plotting is desired, there is no need to run the pipeline again, just subset the file calls/allstats.txt and call scripts/plot_allstats.R . For example, to plot only assembly sequences larger than 5 Gbp  
```
seqtk comp data/genome_assembly.fa | awk '{{if($2>5000000) print $1"\t0\t"$2}}' > bigchroms.bed

bedtools intersect -wa -a calls/allstats.txt -b bigchroms.bed > subs.txt

Rscript scripts/plot_allstats.R subs.txt bigchroms.pdf
```
- cleaning up / archiving: Once the final stats for all desired windowsizes are produced, it makes sense to clean up large intermediate files. I suggest to keep only the config.yaml, the popmap and samples_reads_map, the results_processed, or if you have a bit of space, also results_raw. I would delete the large .BAM files and intermediate VCF files; in the worst case these can be re-calculated. To clean up in this way, run 
```
rm -rf mapped_reads FB_chunks FB_chunk_VCFs FB_chunk_VCFs_filtered FB_regionsALL.bed normalization_coefficients.txt
```
- I would also get rid of the hidden .snakemake directory, unless you really need it. It can be quite large/numerous tiny files.




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

