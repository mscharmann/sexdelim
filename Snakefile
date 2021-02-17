popmapfile = "data/popmap.txt"
windowsize = 500
genomefile = "data/fakegenome.MALE.fa"
samples_reads_map = "data/samples_and_readfiles.txt"


def read_popmap(popmapfile):
	
	SAMPLES = []
	with open(popmapfile, "r") as I:
		for line in I:
			if len(line) > 1:
				fields = line.strip("\n").split("\t")
				SAMPLES.append(fields[0])
	return(SAMPLES)


def parse_samples_reads_map (infile):
	
	samples_files_dict = {}
	with open(infile, "r") as I:
		I.readline() # this file has a header!
		for line in I:
			if len(line) > 1:
				fields = line.strip("\n").split("\t")
				s = fields[0]
				fq1 = fields[1]
				try:
					fq2 = fields[2]
					if len(fq2.strip()) == 0:
						fq2 = None
				except IndexError:
					fq2 = None
				try:
					samples_files_dict[s].append( (fq1, fq2) )
				except KeyError:
					samples_files_dict[s] = [ (fq1, fq2) ] 

	return(samples_files_dict)


def get_fastq(wildcards):
	"""Get fastq files of given sample -- DANGER: currently at most one pair of read files per sample is handled!"""
	fqs = samples_files_dict[wildcards.sample][0]
	if not fqs[1]:
		return [fqs[0]]
	return [fqs[0], fqs[1]]



SAMPLES = read_popmap(popmapfile)
print (SAMPLES)
samples_files_dict = parse_samples_reads_map(samples_reads_map)
print (samples_files_dict)


rule all:
	input:
		"calls/all.post_filter.vcf.gz",
		"calls/windows.multicov.F.txt",
		"calls/windows.multicov.M.txt",
		"calls/M_F_norm_coverage_ratio_log2.bed.txt",
		"calls/plink.assoc_results.significant.window.bed.txt",
		"calls/plink.assoc_results.assoc.raw.txt",
		"calls/LD_r2_averaged_per_window.bed.txt",
		"calls/F_specific_kmer.fa",
		"calls/M_specific_kmer.fa",
		"calls/kmers.F_specific.per_window.bed.txt",
		"calls/kmers.M_specific.per_window.bed.txt",
		"calls/gametolog_candidate_alleles.ZW_divergence.windows.bed.txt",
		"calls/gametolog_candidate_alleles.XY_divergence.windows.bed.txt",
		"calls/MvF.dxy.bed.txt",
		"calls/MvF.fst.bed.txt",
		"calls/F.pi.bed.txt",
		"calls/M.pi.bed.txt",
		"calls/allstats.txt",
		"calls/allstats.plots.pdf"
	shell:
		"""
		DIR="FB_chunks"
		if [ -d "$DIR" ]; then
		  rm -r $DIR
		fi
		if [ -f "FB_regionsALL.bed" ]; then
		  rm FB_regionsALL.bed
		fi
		"""


rule bwa_idx:
	input:
	   genomefile
	output:
		"{genomefile}.bwt"
	shell:
		"""
		if [[ ! $( grep ">" {input} ) =~ "|" ]]; then
			bwa index {input}
		else
			echo "refusing to run, fasta headers contain pipe '|' character, dying"
		fi
		"""

rule bwa_map:
	input:
		fa=genomefile,
		gidx=genomefile + ".bwt",
		reads=get_fastq
	output:
		"mapped_reads/{sample}.bam"
	threads: 14
	run:
		if len(input.reads) == 2: # paired-end!
			shell("""
				# filtering alignments to be primary (i.e. each read only once; if multiple locations equally possible than a random one is chosen): 
				# -F 256 == -F 0x0100 == NOT not primary alignment
				# filtering alignments to be NOT supplementary (supplementary: sections of read map to discontinuous coordinates, e.g. across an inversion breakpoint..):
				# -F 2048 == -F 0x800 == NOT supplementary alignment
				# sum of the bit flags: 2304 => filters against BOTH non-primary and supplementary alignments; verified with samtools flagstat
				# filtering alignments to be "properly paired": -f 2
				bwa mem -t {threads} -a {input.fa} {input.reads[0]} {input.reads[1]} -R "@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tPL:Illumina" | samtools view -F 2304 -f 2 -b -@ 2 - > {output} 
				""")
		else: # single-end
			shell("""
				# filtering alignments to be primary (i.e. each read only once; if multiple locations equally possible than a random one is chosen): 
				# -F 256 == -F 0x0100 == NOT not primary alignment
				# filtering alignments to be NOT supplementary (supplementary: sections of read map to discontinuous coordinates, e.g. across an inversion breakpoint..):
				# -F 2048 == -F 0x800 == NOT supplementary alignment
				# -F 4 read unmapped (0x4)
				# sum of the bit flags: 2308 => filters against non-primary and supplementary alignments and unmapped
				bwa mem -t {threads} -a {input.fa} {input.reads[0]} -R "@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tPL:Illumina" | samtools view -F 2308 -b -@ 2 - > {output} 
				""")


rule samtools_sort:
	input:
		"mapped_reads/{sample}.bam"
	output:
		"mapped_reads/{sample}.sorted.bam"
	shell:
		"""
		samtools sort -T mapped_reads/{wildcards.sample} -O bam {input} > {output}
		rm {input}
		"""
		

rule samtools_index:
	input:
		"mapped_reads/{sample}.sorted.bam"
	output:
		"mapped_reads/{sample}.sorted.bam.bai"
	shell:
		"samtools index {input}"


# A checkpoint that shall trigger re-evaluation of the DAG
# https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#data-dependent-conditional-execution
# we need this because the number of output files (i.e. chunks / region files for Freebayes) is not known before execution!
checkpoint split_ref_for_freebayes:
	input:
		ref=genomefile,
		indexes=expand("mapped_reads/{sample}.sorted.bam.bai", sample=SAMPLES),
		samples=expand("mapped_reads/{sample}.sorted.bam", sample=SAMPLES)
	output:
		directory("FB_chunks")
	params:
		chunksize_Mb = 100000
	shell:
		"""
		samtools faidx {input.ref}
		
		echo {input.samples} | sed 's/ /\\n/g' > bamlistforsplit
		# wget https://gist.githubusercontent.com/travc/0c53df2c8eca81c3ebc36616869930ec/raw/eff3032ca7c955ca33bffd8758092e4006949c75/split_ref_by_bai_datasize.py
		python scripts/split_ref_by_bai_datasize.py -r {input.ref}.fai -L bamlistforsplit --target-data-size {params.chunksize_Mb} > FB_regionsALL.bed
		rm bamlistforsplit
		split --lines=1 FB_regionsALL.bed FB_regions_ --numeric-suffixes --suffix-length=4 --additional-suffix=.bed
		mkdir {output}
		for i in FB_regions_* ; do mv $i {output} ; done		
		"""

rule freebayes:
	input:
		ref=genomefile,
		regions="FB_chunks/{i}.bed"
	output:
		"FB_chunk_VCFs/{i}.bed.vcf.gz"
	params: 
		skip_cov = 6750 # set this to 10 * N_samples * expected average coverage (= bases sequenced per sample / genomesize); e.g. 20*10*(15/0.444)
	shell:
		"""
		freebayes -L <( ls mapped_reads/*.bam ) -t {input.regions} -f {input.ref} --min-coverage 6 --min-mapping-quality 5 \
			--min-base-quality 10 --max-complex-gap 3 --min-repeat-entropy 1 --binomial-obs-priors-off --no-population-priors \
			--use-best-n-alleles 4 --report-monomorphic --skip-coverage {params.skip_cov} | \
			vcfallelicprimitives -kg | \
			tr "|" "/" | \
			bgzip -c > {output}
			
		# unfortunately, the use of vcfallelicprimitives -kg can result in a known and unresolved bug:
		# https://github.com/vcflib/vcflib/issues/225
		
		# ULTRA-important freebayes options for efficiency:
		# -g --skip-coverage N                                                          
        #           Skip processing of alignments overlapping positions with coverage >N.
        #           This filters sites above this coverage, but will also reduce data nearby.
        #           default: no limit                                               
		# https://groups.google.com/g/freebayes/c/5L657BepjpY
		# Erik Garrison writes: "This will jump over sites that have total coverage in all samples >N. This is important when processing genomes with collapsed repeats, 
		# where coverage may go very high. In many genomes, it makes sense to set this to a multiple of your expected read coverage (e.g. 5-10x, or 
		# alternatively several standard deviations above the mean). Variant calls made in these regions of collapsed sequence are unlikely to be "correct" 
		# in that the actual local copy number can be very high, and so the genotype won't be the ploidy that you've set."
		#
		#  WHAT EG DOES NOT SAY: skip-coverage evaluates the TOTAL SUM OVER ALL SAMPLES
		#
		#  --limit-coverage N                                                    
        #           Downsample per-sample coverage to this level if greater than this coverage.
        #           default: no limit
        #
        # EG: "So this will pseudorandomly downsample each biosample until its max coverage is N. If you have a lot of samples, the coverage will still be higher.
        # I think it may not adjust DP. Look at the per-sample DP to see if it changes."
        #
        # WHAT EG DOES NOT SAY: limit-coverage does not help or only small reduction in running time + memory load; 
        # apparently all the reads still get loaded to memory before they get downsampled!?                                              
		"""


def merge_vcfs_input(wildcards):
	checkpoint_output = checkpoints.split_ref_for_freebayes.get(**wildcards).output[0]
	# print(checkpoint_output)
	these_files = sorted( expand("FB_chunk_VCFs/{i}.bed.vcf.gz", i=glob_wildcards(os.path.join(checkpoint_output, "{i}.bed")).i) )
	# print (thef)
	return these_files
           
           

rule merge_vcfs:
	input:
		region_vcfs=merge_vcfs_input # refers to the function above which evaluates the checkpoint
	output:
		"calls/all.pre_filter.vcf.gz"
	threads: 3
	shell:
		"""
		# MUST NOT USE bcftools concat: it cannot resolve POS that are non-monotonically icreasing (which ca happen at the interval boundaries)
		# the code below is only slightly modified from freebayes-parallel script: https://github.com/freebayes/freebayes/blob/master/scripts/freebayes-parallel
		# zcat {input.region_vcfs} | python2.7 $(which vcffirstheader) | vcfstreamsort -w 10000 | vcfuniq | bgzip -c > {output}
		# zcat alone may complain about too many arguments, so better use find -exec :
		find FB_chunk_VCFs/*.bed.vcf.gz -type f -exec zcat {{}} \\; | python2.7 $(which vcffirstheader) | vcfstreamsort -w 10000 | vcfuniq | bgzip -c > {output}
		"""


rule VCF_get_mean_depth_per_site:
	input:
		"calls/all.pre_filter.vcf.gz"
	output:
		"calls/meandepthpersite.txt"
	shell:
		"""
		# get a subset of 0.1% of the VCF records, which is used to estimate mean depth per site
		bcftools view {input} | vcfrandomsample -r 0.001 | bgzip -c > subset.vcf.gz

		# Calculate mean depth per site
		vcftools --gzvcf subset.vcf.gz --site-mean-depth --stdout | cut -f3 > {output}
		
		rm subset.vcf.gz
		"""

	
rule VCF_filter_variants:
	input:
		mdps="calls/meandepthpersite.txt",
		gzvcf="calls/all.pre_filter.vcf.gz"
	output:
		"calls/variant_sites.vcf.gz"
	params:
		MISS=0.1,
		QUAL=20,
		MIN_DEPTH=6
	shell:
		"""
		## https://bcbio.wordpress.com/2013/10/21/updated-comparison-of-variant-detection-methods-ensemble-freebayes-and-minimal-bam-preparation-pipelines/#comment-1469
		## => "We use ... a very vanilla filter with only depth and quality (QUAL < 20, DP < 5)"
		## => impose a minimum of QUAL, and a minimum of DEPTH
		## 	for variants, I want those two, AND a maximum DPETH cutoff at 95% of the mean depth over all samples
		## 	for invariants, we can only impose a minimum DEPTH and a maximum DEPTH.
		## use code from dDocent to calcualte a mean depth histogram, then find the 95th percentile: no. too complicated; there are easier ways to do this.
		## also good info:
		## 	https://speciationgenomics.github.io/filtering_vcfs/

		# set filters
		# MAX_DEPTH is specific to the dataset; we take the 95th percentile of the mean depths per site:
		MAX_DEPTH=$( sort -n {input.mdps} | awk 'BEGIN{{c=0}} length($0){{a[c]=$0;c++}}END{{p5=(c/100*5); p5=p5%1?int(p5)+1:p5; print a[c-p5-1]}}' )

		echo $MAX_DEPTH

		# variants:
		vcftools --gzvcf {input.gzvcf} --mac 1 --minQ {params.QUAL} --min-meanDP {params.MIN_DEPTH} --minDP {params.MIN_DEPTH} \
		--max-meanDP $MAX_DEPTH --max-missing {params.MISS} --recode --stdout | bgzip -c > {output}
		tabix {output}
		"""

rule VCF_filter_invariants:
	input:
		mdps="calls/meandepthpersite.txt",
		gzvcf="calls/all.pre_filter.vcf.gz"
	output:
		"calls/invariant_sites.vcf.gz"
	params:
		MISS=0.1,
		QUAL=20,
		MIN_DEPTH=6
	shell:
		"""
		## https://bcbio.wordpress.com/2013/10/21/updated-comparison-of-variant-detection-methods-ensemble-freebayes-and-minimal-bam-preparation-pipelines/#comment-1469
		## => "We use ... a very vanilla filter with only depth and quality (QUAL < 20, DP < 5)"
		## => impose a minimum of QUAL, and a minimum of DEPTH
		## 	for variants, I want those two, AND a maximum DPETH cutoff at 95% of the mean depth over all samples
		## 	for invariants, we can only impose a minimum DEPTH and a maximum DEPTH.
		## use code from dDocent to calcualte a mean depth histogram, then find the 95th percentile: no. too complicated; there are easier ways to do this.
		## also good info:
		## 	https://speciationgenomics.github.io/filtering_vcfs/

		# set filters
		# MAX_DEPTH is specific to the dataset; we take the 95th percentile of the mean depths per site:
		MAX_DEPTH=$( sort -n {input.mdps} | awk 'BEGIN{{c=0}} length($0){{a[c]=$0;c++}}END{{p5=(c/100*5); p5=p5%1?int(p5)+1:p5; print a[c-p5-1]}}' )

		echo $MAX_DEPTH

		# variants:
		vcftools --gzvcf {input.gzvcf} --max-maf 0 --min-meanDP {params.MIN_DEPTH} --minDP {params.MIN_DEPTH} \
		--max-meanDP $MAX_DEPTH --max-missing {params.MISS} --recode --stdout | bgzip -c > {output}		
		tabix {output}
		"""

rule VCF_merge_post_filter:
	input:
		invar="calls/invariant_sites.vcf.gz",
		var="calls/variant_sites.vcf.gz"
	output:
		"calls/all.post_filter.vcf.gz"
	shell:
		"""
		# combine the two VCFs using bcftools concat
		bcftools concat --allow-overlaps {input.var} {input.invar} -O z -o {output}
		"""
		
rule get_coverage_in_windows:
	input:
		fa=genomefile,
		popmap={popmapfile},
		bam=expand("mapped_reads/{sample}.sorted.bam", sample=SAMPLES),
		bai=expand("mapped_reads/{sample}.sorted.bam.bai", sample=SAMPLES)
	output:
		fcov="calls/windows.multicov.F.txt",
		mcov="calls/windows.multicov.M.txt",
		log2ratio="calls/M_F_norm_coverage_ratio_log2.bed.txt"
	threads: 2
	shell:
		"""
		seqtk comp {input.fa} | awk '{{print $1"\t"$2}}' > genomefile.txt
		bedtools makewindows -w {windowsize} -g genomefile.txt > windows.bed
		rm genomefile.txt
		
		# separate calls for M and F	
		cat {input.popmap} | awk 'NF {{if($2==1) print $1}}' | tr "\\n" "\\t" | awk 'BEGIN {{ OFS="\\t"}}; {{$1=$1; print "chr","start","stop",$0}}' > {output.mcov}
		cat {input.popmap} | awk 'NF {{if($2==2) print $1}}' | tr "\\n" "\\t" | awk 'BEGIN {{ OFS="\\t"}}; {{$1=$1; print "chr","start","stop",$0}}' > {output.fcov}
		( bedtools multicov -bams $( cat {input.popmap} | awk 'BEGIN {{ ORS=" "}}; NF {{if($2==1) print "mapped_reads/"$1".sorted.bam"}}' ) -bed windows.bed >> {output.mcov} )&
		( bedtools multicov -bams $( cat {input.popmap} | awk 'BEGIN {{ ORS=" "}}; NF {{if($2==2) print "mapped_reads/"$1".sorted.bam"}}' ) -bed windows.bed >> {output.fcov} )&
		wait
		rm windows.bed
		
		################ get M-F norm read coverage ratio log2 (we add 1e-6 to the read counts to avoid zero division)
		awk '{{ for(i=4; i<=NF;i++) j+=$i; print j; j=0 }}' {output.mcov} > total.M
		totm=$( awk '{{sum+=$1;}} END{{print sum;}}' total.M )
		awk -v totalsum="$totm" '{{print ($1*1000000)/totalsum}}' total.M > total.norm.M

		awk '{{ for(i=4; i<=NF;i++) j+=$i; print j; j=0 }}' {output.fcov} > total.F
		totf=$( awk '{{sum+=$1;}} END{{print sum;}}' total.F )
		awk -v totalsum="$totf" '{{print ($1*1000000)/totalsum}}' total.F > total.norm.F

		paste total.norm.M total.norm.F | awk '{{print log((($1+1e-6)/($2+1e-6)))/log(2)}}' > log2_ratio

		paste <(cut -f1,2,3 {output.fcov} ) log2_ratio | tail -n +2 > {output.log2ratio}

		rm total.M total.F total.norm.M total.norm.F log2_ratio 	
		"""

rule pixy_stats:
	input:
		popmap={popmapfile},
		gzvcf="calls/all.post_filter.vcf.gz"
	output:
		dxy="calls/MvF.dxy.full_report.txt",
		dxy_bed="calls/MvF.dxy.bed.txt",
		fst="calls/MvF.fst.full_report.txt",
		fst_bed="calls/MvF.fst.bed.txt",
		pi="calls/M_F.pi.full_report.txt",
		pi_f_bed="calls/F.pi.bed.txt",
		pi_m_bed="calls/M.pi.bed.txt"		
	shell:
		"""
		# this sometimes throws an error:
		# raise ValueError('values must be monotonically increasing')
		# ValueError: values must be monotonically increasing
		
		
		# first, an ugly workaround to avoid an issue caused by vcfallelicprimitives duplicating an INFO field..
		bcftools view {input.gzvcf} | grep -v "##INFO=<ID=NUMALT" | bgzip -c > vcf_noheader.vcf.gz

		DIR="tmp_zarr"
		if [ -d "$DIR" ]; then
		  rm -r $DIR
		fi
		
		mkdir tmp_zarr
		pixy --stats pi fst dxy \
		--vcf vcf_noheader.vcf.gz \
		--zarr_path tmp_zarr \
		--window_size {windowsize} \
		--populations {input.popmap} \
		--bypass_filtration yes
		rm -r tmp_zarr
		rm vcf_noheader.vcf.gz
		
		mv pixy_output_dxy.txt {output.dxy}
		mv pixy_output_fst.txt {output.fst}
		mv pixy_output_pi.txt {output.pi}
		
		cat {output.dxy} | awk 'OFS="\\t" {{print $3,$4,$5,$6}}' | tail -n +2 > {output.dxy_bed}
		cat {output.fst} | awk 'OFS="\\t" {{print $3,$4,$5,$6}}' | tail -n +2 > {output.fst_bed}
	
		cat {output.pi} | awk 'OFS="\t" {{if($1 == 1) print $2,$3,$4,$5}}' | tail -n +2 > {output.pi_m_bed}
		cat {output.pi} | awk 'OFS="\t" {{if($1 == 2) print $2,$3,$4,$5}}' | tail -n +2 > {output.pi_f_bed}
		
		"""

	
rule GWAS_plink:
	input:
		popmap={popmapfile},
		gzvcf="calls/variant_sites.vcf.gz",
		fa=genomefile
	output:
		plinkassoc="calls/plink.assoc_results.significant.window.bed.txt",
		rawplink="calls/plink.assoc_results.assoc.raw.txt"
	shell:
		"""
		cat {input.popmap} | awk 'NF {{print $1"\\t"$1"\\t"$2}}' > thepheno
		plink --vcf {input.gzvcf} --double-id --allow-extra-chr --pheno thepheno --allow-no-sex --assoc --out assoc_results	
		rm thepheno
		
		cat assoc_results.assoc | tr -s " " "\\t" | awk '{{if($9 <= 0.05) print}}' > plink.assoc_results.significant.txt
		mv assoc_results.assoc {output.rawplink}
		rm assoc_results*

		cat plink.assoc_results.significant.txt | cut -f1,3,9 | awk '{{print $1"\\t"$2"\\t"$2}}' > plink.assoc_results.significant.bed
		rm plink.assoc_results.significant.txt

		seqtk comp {input.fa} | awk '{{print $1"\\t"$2}}' > genomefile.assoc.txt
		bedtools makewindows -w {windowsize} -g genomefile.assoc.txt > windows.assoc.bed
		rm genomefile.assoc.txt 
		bedtools coverage -counts -a windows.assoc.bed -b plink.assoc_results.significant.bed > {output.plinkassoc}
		rm windows.assoc.bed plink.assoc_results.significant.bed
		"""

rule LD_plink:
	input:
		gzvcf="calls/variant_sites.vcf.gz",
		fa=genomefile
	output:
		"calls/LD_r2_averaged_per_window.bed.txt"
	shell:
		"""
		# apply a moderate MAF filter before calculating LD as rare alleles imply very high LD by definition; they are thus useless for our purpose.
		# vcftools --gzvcf {input.gzvcf} --maf 0.1 --recode --stdout > forplinkld.vcf
		# plink --vcf forplinkld.vcf --double-id --allow-extra-chr --r2 --ld-window-r2 0.0 --ld-window 5 --ld-window-kb 2 
		# rm forplinkld.vcf 
		
		plink --vcf {input.gzvcf} --double-id --allow-extra-chr --maf 0.1 --r2 --ld-window-r2 0.0 --ld-window 5 --ld-window-kb 2 
		
		# --ld-window X	compute LD only for pairs that are at most X SNPs apart (default 10)
		# --ld-window-kb X	compute LD only for pairs that are at most X kb apart (default 1000 kb)
		# --ld-window-r2 X	minimum r2 to report, else omit from output. Default = 0.2

		rm plink.log plink.nosex

		cat plink.ld | tr -s ' ' '\\t' | cut -f1,2,4,5,7 | tail -n +2 | awk '{{ print $1"\\t"$2"\\t"$2"\\t"$5"\\n"$3"\\t"$4"\\t"$4"\\t"$5  }}' > ld_clean.bed

		bedtools sort -i ld_clean.bed > ld_clean.sorted.bed

		seqtk comp {input.fa} | awk '{{print $1"\\t"$2}}' > genomefile.ld.txt
		bedtools makewindows -w {windowsize} -g genomefile.ld.txt > windows.ld.bed
		rm genomefile.ld.txt 
		bedtools map -a windows.ld.bed -b ld_clean.sorted.bed -c 4 -o mean > {output}
		rm ld_clean.sorted.bed windows.ld.bed ld_clean.bed plink.ld
		"""
	
	
rule kmerGO:
	input:
		popm={popmapfile},
		SamReMap={samples_reads_map}
	output:
		fkmer="calls/F_specific_kmer.fa",
		mkmer="calls/M_specific_kmer.fa",
		fcont="calls/F_specific.fa.cap.contigs",
		mcont="calls/M_specific.fa.cap.contigs"
	threads: 14
	shell:
		"""
		# for kmerGO to recognise input read files, there must be no "." characters in the filename except in the suffix; 
		# allowed suffixes are: FASTQ_suffix = ['.fq', '.fastq', '.fq.gz', '.fastq.gz']
		# also, we must concatenate all read files to a single file per sample... 
		cat {input.popm} | awk 'NF {{if($2==1) print $1",M"; else print $1",F"}}' > traits_sex_for_kmerGO.txt
		tail -n +2 {input.SamReMap} > nohead
		mkdir read_files
		cd read_files
		while read line; do
			sid=$( echo $line | awk '{{print $1}}' )
			f1=$( echo $line | awk '{{print $2}}')
			f2=$( echo $line | awk '{{print $3}}' )
			echo $sid $f1 $f2
			if [ -z "$f2" ]
			then
				cat ../$f1 >> ./$sid.fastq.gz 
			else
				cat ../$f1 >> ./$sid.fastq.gz
				cat ../$f2 >> ./$sid.fastq.gz	
			fi
		done < ../nohead
		cd ../
		
		rm nohead
		
		# 16-mers
		# -ci 	/ 	No minimal K-mer occurring times (default: 2)
		scripts/KmerGO_for_linux_x64_cmd/KmerGO -n {threads} -k 16 -ci 3 -i read_files -t traits_sex_for_kmerGO.txt -p 0.01 -assl 0.8 -assn 0.8
		rm -r kmer_features kmer_matrix read_files kmer_countings traits_sex_for_kmerGO.txt
		mv contig_result/F_specific_kmer.fa {output.fkmer}
		mv contig_result/M_specific_kmer.fa {output.mkmer}
		mv contig_result/F_specific.fa.cap.contigs {output.fcont}
		mv contig_result/M_specific.fa.cap.contigs {output.mcont}
		rm -r contig_result
		
		"""
	
rule match_kmers_to_genome:
	input:
		popm={popmapfile},
		fa=genomefile,
		gidx=genomefile +".bwt",
		fkmers="calls/F_specific_kmer.fa",
		mkmers="calls/M_specific_kmer.fa"
	output:
		fsp="calls/kmers.F_specific.per_window.bed.txt",
		msp="calls/kmers.M_specific.per_window.bed.txt"
	shell:
		"""
		# searching for perfect matches of 21-mers: https://bioinformatics.stackexchange.com/questions/7298/mapping-heteryzygous-kmers-on-a-genome
		# The -k 21 says to use a minimum seed length of 21, which forces exact matches. The -T 21 requires a minimum score of 21, which also enforces an exact match. The -a parameter reports all matches, since a "best" match doesn't make sense in this situation. The -c parameter limits how many matches are reported, which may need to be adjusted depending on how repetitive the 21-mer is.
		# we match IMPERFECTLY = allowing 1 mismatch in 16 bp kmers. => 6.25% divergence
		
		# use samtools view -F 4 : removes unmapped
	
		bwa mem -t 1 -k 15 -T 15 -a -c 5000 {input.fa} {input.fkmers} | samtools view -F 4 -b - > F_specific_kmer.bam
		samtools sort -T F_specific_kmer -O bam F_specific_kmer.bam > F_specific_kmer.sorted.bam
		samtools index F_specific_kmer.sorted.bam
		rm F_specific_kmer.bam

		bwa mem -t 1 -k 15 -T 15 -a -c 5000 {input.fa} {input.mkmers} | samtools view -F 4 -b - > M_specific_kmer.bam
		samtools sort -T M_specific_kmer -O bam M_specific_kmer.bam > M_specific_kmer.sorted.bam
		samtools index M_specific_kmer.sorted.bam
		rm M_specific_kmer.bam

		seqtk comp {input.fa} | awk '{{print $1"\\t"$2}}' > genomefile.kmer.txt
		bedtools makewindows -w {windowsize} -g genomefile.kmer.txt > windows.kmer.bed
		rm genomefile.kmer.txt
		bedtools multicov -bams F_specific_kmer.sorted.bam -bed windows.kmer.bed > {output.fsp}
		bedtools multicov -bams M_specific_kmer.sorted.bam -bed windows.kmer.bed > {output.msp}
		rm windows.kmer.bed M_specific_kmer.sorted.bam M_specific_kmer.sorted.bam.bai F_specific_kmer.sorted.bam F_specific_kmer.sorted.bam.bai
		"""
	
	
rule calculate_freqs:
	input:
		popm={popmapfile},
		gzvcf="calls/all.post_filter.vcf.gz"
	output:
		mfreq="calls/freq_males.txt.gz",
		ffreq="calls/freq_females.txt.gz"
	shell:
		"""
		vcftools --gzvcf {input.gzvcf} --keep <( cat {input.popm} | awk '{{if($2==1) print $1}}' ) --freq --stdout | gzip > {output.mfreq}
		vcftools --gzvcf {input.gzvcf} --keep <( cat {input.popm} | awk '{{if($2==2) print $1}}' ) --freq --stdout | gzip > {output.ffreq}
		"""
	
rule pseudo_phase_gametologs:
	input:
		mfreq="calls/freq_males.txt.gz",
		ffreq="calls/freq_females.txt.gz",
		gzvcf="calls/all.post_filter.vcf.gz"
	output:
		"calls/gametolog_candidate_alleles.allsites.vcf.gz"
	shell:
		"""
		python scripts/pseudo_phase_gametologs.py {input.mfreq} {input.ffreq} {input.gzvcf} | bgzip -c > {output}
		"""
	
rule score_XY_gametolog_divergence:
	input:
		"calls/gametolog_candidate_alleles.allsites.vcf.gz"
	output:
		full="calls/gametolog_candidate_alleles.XY_divergence.windows.full_report.txt",
		bedformat="calls/gametolog_candidate_alleles.XY_divergence.windows.bed.txt"
	shell:
		"""
		echo -e "XY_Y_like\\tXY_Y_like" > xypopm
		echo -e "XY_X_like\\tXY_X_like" >> xypopm

		DIR="tmp_zarr_xy"
		if [ -d "$DIR" ]; then
		  rm -r $DIR
		fi
		
		mkdir tmp_zarr_xy
		pixy --stats dxy \
		--vcf {input} \
		--zarr_path tmp_zarr_xy \
		--window_size {windowsize} \
		--populations xypopm \
		--bypass_filtration yes \
		--outfile_prefix ./XY
		rm -r tmp_zarr_xy xypopm

		mv XY_dxy.txt {output.full}
		cat {output.full} | awk 'OFS="\\t" {{print $3,$4,$5,$6}}' | tail -n +2 > {output.bedformat}	
		"""

rule score_ZW_gametolog_divergence:
	input:
		"calls/gametolog_candidate_alleles.allsites.vcf.gz"
	output:
		full="calls/gametolog_candidate_alleles.ZW_divergence.windows.full_report.txt",
		bedformat="calls/gametolog_candidate_alleles.ZW_divergence.windows.bed.txt"
	shell:
		"""
		echo -e "ZW_W_like\\tZW_W_like" > zwpopm
		echo -e "ZW_Z_like\\tZW_Z_like" >> zwpopm

		DIR="tmp_zarr_zw"
		if [ -d "$DIR" ]; then
		  rm -r $DIR
		fi
		
		mkdir tmp_zarr_zw
		pixy --stats dxy \
		--vcf {input} \
		--zarr_path tmp_zarr_zw \
		--window_size {windowsize} \
		--populations zwpopm \
		--bypass_filtration yes \
		--outfile_prefix ./ZW
		rm -r tmp_zarr_zw zwpopm

		mv ZW_dxy.txt {output.full}
		cat {output.full} | awk 'OFS="\\t" {{print $3,$4,$5,$6}}' | tail -n +2 > {output.bedformat}	
		"""

rule plot_all:
	input:
		a="calls/M_F_norm_coverage_ratio_log2.bed.txt",
		b="calls/kmers.F_specific.per_window.bed.txt",
		c="calls/kmers.M_specific.per_window.bed.txt",
		d="calls/F.pi.bed.txt",
		e="calls/M.pi.bed.txt",
		f="calls/MvF.fst.bed.txt",
		g="calls/MvF.dxy.bed.txt",
		h="calls/LD_r2_averaged_per_window.bed.txt",
		i="calls/plink.assoc_results.significant.window.bed.txt",
		j="calls/gametolog_candidate_alleles.ZW_divergence.windows.bed.txt",
		k="calls/gametolog_candidate_alleles.XY_divergence.windows.bed.txt"
	output:
		stats="calls/allstats.txt",
		plot="calls/allstats.plots.pdf"
	shell:
		"""
		cat {input.a} > {output.stats}

		for i in {input.b} {input.c} {input.d} {input.e} {input.f} {input.g} {input.h} {input.i} {input.j} {input.k} ; do
			paste {output.stats} <(cut -f4 $i ) > tmpst
			mv tmpst {output.stats}
		done
		
		Rscript scripts/plot_allstats.R {output.stats} {output.plot}	
		"""


	
