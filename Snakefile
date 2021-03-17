configfile: "config.yaml"

popmapfile = config["popmapfile"]
windowsize = int( config["windowsize"] )
genomefile = config["genomefile"]
samples_reads_map = config["samples_reads_map"]
regions_for_plot_bed = config["regions_for_plot_bed"]


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
		"results_raw/all.post_filter.vcf.gz",
		"results_processed/windows.multicov.F.txt",
		"results_processed/windows.multicov.M.txt",
		"results_processed/M_F_norm_coverage_ratio_log2.bed.txt",
		"results_processed/plink.assoc_results.significant.window.bed.txt",
		"results_raw/plink.assoc_results.assoc.raw.txt",
		"results_processed/LD_r2_averaged_per_window.bed.txt",
		"results_raw/F_specific_kmer.fa",
		"results_raw/M_specific_kmer.fa",
		"results_processed/kmers.F_specific.per_window.bed.txt",
		"results_processed/kmers.M_specific.per_window.bed.txt",
		"results_processed/gametolog_candidate_alleles.ZW_divergence.windows.bed.txt",
		"results_processed/gametolog_candidate_alleles.XY_divergence.windows.bed.txt",
		"results_processed/MvF.dxy.bed.txt",
		"results_processed/MvF.fst.bed.txt",
		"results_processed/F.pi.bed.txt",
		"results_processed/M.pi.bed.txt",
		"results_processed/MvF.netdiv_M.bed.txt",
		"results_processed/MvF.netdiv_F.bed.txt",
		"results_processed/het_males.bed.txt",
		"results_processed/het_females.bed.txt",
		"results_processed/allstats.txt",
		"results_processed/allstats.plots.pdf",
		"results_processed/region.stats.plots.pdf"



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
		chunksize_Mb = config["freebayes_chunksize"]
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
		skip_cov = config["freebayes_skip_cov"] # set this to 10 * N_samples * expected average coverage (= bases sequenced per sample / genomesize); e.g. 20*10*(15/0.444)
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
		"results_raw/all.pre_filter.vcf.gz"
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
		"results_raw/all.pre_filter.vcf.gz"
	output:
		"results_raw/meandepthpersite.txt"
	threads: 2
	shell:
		"""
		# index VCF
		( tabix {input} )&
		
		# get a subset of 10k to 100k of the raw VCF records; iteratively keeping only 10% of lines in file
		( zcat {input} | cut -f1,2 > all_pos_in_vcf.txt 
		
		cp all_pos_in_vcf.txt tmp1
		nlines=$( cat tmp1 | wc -l )
		echo $nlines
		while [ "$nlines " -gt 100000 ] ; do
		awk 'NR == 1 || NR % 10 == 0' tmp1 > tmp2
		mv tmp2 tmp1
		nlines=$( cat tmp1 | wc -l )
		echo $nlines
		done
		)&
		wait
		
		grep -v "#" tmp1 > postokeep.txt
	
		# subset the VCF
		tabix -h -R postokeep.txt {input} > subset.vcf
		
		# Calculate mean depth per site
		vcftools --vcf subset.vcf --site-mean-depth --stdout | cut -f3 > {output}
		
		rm subset.vcf tmp1 all_pos_in_vcf.txt postokeep.txt 
		"""

	
rule VCF_filter_variants:
	input:
		mdps="results_raw/meandepthpersite.txt",
		gzvcf="results_raw/all.pre_filter.vcf.gz"
	output:
		"results_raw/variant_sites.vcf.gz"
	params:
		MISS=config["VCF_MISS"],
		QUAL=config["VCF_QUAL"],
		MIN_DEPTH=config["VCF_MIN_DEPTH"]
	shell:
		"""
		## https://bcbio.wordpress.com/2013/10/21/updated-comparison-of-variant-detection-methods-ensemble-freebayes-and-minimal-bam-preparation-pipelines/#comment-1469
		## => "We use ... a very vanilla filter with only depth and quality (QUAL < 20, DP < 5)"
		## => impose a minimum of QUAL, and a minimum of DEPTH
		## 	for variants, I want those two, AND a maximum DPETH cutoff 
		## 	for invariants, we can only impose a minimum DEPTH and a maximum DEPTH.
		## use code from dDocent to calcualte a mean depth histogram, then find the 95th percentile: no. too complicated; there are easier ways to do this.
		## also good info:
		## 	https://speciationgenomics.github.io/filtering_vcfs/

		# set filters
		# MAX_DEPTH is specific to the dataset; we take the 98th percentile of the mean depths per site:
		MAX_DEPTH=$( sort -n {input.mdps} | awk 'BEGIN{{c=0}} length($0){{a[c]=$0;c++}}END{{p2=(c/100*2); p2=p2%1?int(p2)+1:p2; print a[c-p2-1]}}' )

		echo $MAX_DEPTH

		# variants:
		vcftools --gzvcf {input.gzvcf} --mac 1 --minQ {params.QUAL} --min-meanDP {params.MIN_DEPTH} --minDP {params.MIN_DEPTH} \
		--max-meanDP $MAX_DEPTH --max-missing {params.MISS} --recode --stdout | bgzip -c > {output}
		tabix {output}
		"""

rule VCF_filter_invariants:
	input:
		mdps="results_raw/meandepthpersite.txt",
		gzvcf="results_raw/all.pre_filter.vcf.gz"
	output:
		"results_raw/invariant_sites.vcf.gz"
	params:
		MISS=config["VCF_MISS"],
		QUAL=config["VCF_QUAL"],
		MIN_DEPTH=config["VCF_MIN_DEPTH"]
	shell:
		"""
		## https://bcbio.wordpress.com/2013/10/21/updated-comparison-of-variant-detection-methods-ensemble-freebayes-and-minimal-bam-preparation-pipelines/#comment-1469
		## => "We use ... a very vanilla filter with only depth and quality (QUAL < 20, DP < 5)"
		## => impose a minimum of QUAL, and a minimum of DEPTH
		## 	for variants, I want those two, AND a maximum DPETH cutoff 
		## 	for invariants, we can only impose a minimum DEPTH and a maximum DEPTH.
		## use code from dDocent to calcualte a mean depth histogram, then find the 95th percentile: no. too complicated; there are easier ways to do this.
		## also good info:
		## 	https://speciationgenomics.github.io/filtering_vcfs/

		# set filters
		# MAX_DEPTH is specific to the dataset; we take the 98th percentile of the mean depths per site:
		MAX_DEPTH=$( sort -n {input.mdps} | awk 'BEGIN{{c=0}} length($0){{a[c]=$0;c++}}END{{p2=(c/100*2); p2=p2%1?int(p2)+1:p2; print a[c-p2-1]}}' )

		echo $MAX_DEPTH

		# variants:
		vcftools --gzvcf {input.gzvcf} --max-maf 0 --min-meanDP {params.MIN_DEPTH} --minDP {params.MIN_DEPTH} \
		--max-meanDP $MAX_DEPTH --max-missing {params.MISS} --recode --stdout | bgzip -c > {output}		
		tabix {output}
		"""

rule VCF_merge_post_filter:
	input:
		invar="results_raw/invariant_sites.vcf.gz",
		var="results_raw/variant_sites.vcf.gz"
	output:
		"results_raw/all.post_filter.vcf.gz"
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
		fcov="results_processed/windows.multicov.F.txt",
		mcov="results_processed/windows.multicov.M.txt",
		log2ratio="results_processed/M_F_norm_coverage_ratio_log2.bed.txt"
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

rule pi_rawstats_M:
	input:
		popmap={popmapfile},
		gzvcf="results_raw/all.post_filter.vcf.gz"
	output:
		"results_raw/pi_raw.M.txt"		
	shell:
		"""
		cat {input.popmap} | awk '{{if($2==1) print $1}}' > mpop_pi
		vcftools --gzvcf {input.gzvcf} --keep mpop_pi --max-missing 0.01 --recode --stdout | python2.7 scripts/make_denominator_and_numerator_for_pi.py > {output}
		rm mpop_pi
		"""

rule pi_rawstats_F:
	input:
		popmap={popmapfile},
		gzvcf="results_raw/all.post_filter.vcf.gz"
	output:
		"results_raw/pi_raw.F.txt"		
	shell:
		"""
		cat {input.popmap} | awk '{{if($2==2) print $1}}' > fpop_pi
		vcftools --gzvcf {input.gzvcf} --keep fpop_pi --max-missing 0.01 --recode --stdout | python2.7 scripts/make_denominator_and_numerator_for_pi.py > {output}
		rm fpop_pi
		"""

rule pi_windowed:
	input:
		pim="results_raw/pi_raw.M.txt",
		pif="results_raw/pi_raw.F.txt",
		fa=genomefile
	output:
		pi_f_bed="results_processed/F.pi.bed.txt",
		pi_m_bed="results_processed/M.pi.bed.txt"		
	threads: 4
	shell:
		"""
		seqtk comp {input.fa} | awk '{{print $1"\\t"$2}}' > genomefile.pi.txt
		bedtools makewindows -w {windowsize} -g genomefile.pi.txt > windows.pi.bed

		# ugly construct to handle the "bug" that chromosomes in the VCF are not sorted in 
		# same order as in the fasta genome file; its unclear how this can happen but it DOES
		cut -f1 {input.pim} | uniq | awk '{{print $1"\\t0\\t9999999999999999"}}' > pi_sort_order_in_scores.txt
		# find chroms NOT YET in the list
		comm -13 <(cut -f1 pi_sort_order_in_scores.txt | sort | uniq ) <(cut -f1 genomefile.pi.txt | sort | uniq) | awk '{{print $1"\\t0\\t9999999999999999"}}' >> pi_sort_order_in_scores.txt
			 
		bedtools sort -g pi_sort_order_in_scores.txt -i windows.pi.bed > windows.pi.sort_order_in_scores.bed
		
		# MALE
		(bedtools map -a windows.pi.sort_order_in_scores.bed -b {input.pim} -c 4 -o mean > pi_num_M )& 
		(bedtools map -a windows.pi.sort_order_in_scores.bed -b {input.pim} -c 5 -o mean > pi_denom_M )&

		# FEMALE
		(bedtools map -a windows.pi.sort_order_in_scores.bed -b {input.pif} -c 4 -o mean > pi_num_F )&
		(bedtools map -a windows.pi.sort_order_in_scores.bed -b {input.pif} -c 5 -o mean > pi_denom_F )&
		
		wait
		
		paste pi_num_M pi_denom_M > pi_both_M		
		awk '{{ if($8>0) print $1"\\t"$2"\\t"$3"\\t"$4/$8 ; else print $1"\\t"$2"\\t"$3"\\tNA" }}' pi_both_M > pi_tmp 			
		bedtools sort -g genomefile.pi.txt -i pi_tmp > {output.pi_m_bed}

		paste pi_num_F pi_denom_F > pi_both_F		
		awk '{{ if($8>0) print $1"\\t"$2"\\t"$3"\\t"$4/$8 ; else print $1"\\t"$2"\\t"$3"\\tNA" }}' pi_both_F > pi_tmp		
		bedtools sort -g genomefile.pi.txt -i pi_tmp > {output.pi_f_bed}
		
		rm genomefile.pi.txt windows.pi.bed pi_num_M pi_denom_M pi_num_F pi_denom_F pi_both_M pi_both_F pi_tmp windows.pi.sort_order_in_scores.bed pi_sort_order_in_scores.txt
		"""

rule dxy_rawstats:
	input:
		popmap={popmapfile},
		gzvcf="results_raw/all.post_filter.vcf.gz"
	output:
		"results_raw/dxy_raw.txt"		
	shell:
		"""
		zcat {input.gzvcf} | python2.7 scripts/make_denominator_and_numerator_for_dxy.py --popmap {input.popmap} > {output}
		"""

rule dxy_windowed:
	input:
		raw="results_raw/dxy_raw.txt",
		fa=genomefile
	output:
		"results_processed/MvF.dxy.bed.txt"
	threads: 2
	shell:
		"""
		seqtk comp {input.fa} | awk '{{print $1"\\t"$2}}' > genomefile.dxy.txt
		bedtools makewindows -w {windowsize} -g genomefile.dxy.txt > windows.dxy.bed
		
		# ugly construct to handle the "bug" that chromosomes in the VCF are not sorted in 
		# same order as in the fasta genome file; its unclear how this can happen but it DOES
		cut -f1 {input.raw} | uniq | awk '{{print $1"\\t0\\t9999999999999999"}}' > dxy_sort_order_in_scores.txt
		# find chroms NOT YET in the list
		comm -13 <(cut -f1 dxy_sort_order_in_scores.txt | sort | uniq ) <(cut -f1 genomefile.dxy.txt | sort | uniq) | awk '{{print $1"\\t0\\t9999999999999999"}}' >> dxy_sort_order_in_scores.txt

		bedtools sort -g dxy_sort_order_in_scores.txt -i windows.dxy.bed > windows.dxy.sort_order_in_scores.bed

		# correct
		( bedtools map -a windows.dxy.sort_order_in_scores.bed -b {input.raw} -c 4 -o mean > dxy_num )&
		( bedtools map -a windows.dxy.sort_order_in_scores.bed -b {input.raw} -c 5 -o mean > dxy_denom )&
		wait

		paste dxy_num dxy_denom > dxy_both
		
		awk '{{ if($8>0) print $1"\\t"$2"\\t"$3"\\t"$4/$8 ; else print $1"\\t"$2"\\t"$3"\\tNA" }}' dxy_both > dxy_tmp
			
		bedtools sort -g genomefile.dxy.txt -i dxy_tmp > {output}
		
		rm genomefile.dxy.txt windows.dxy.bed dxy_num dxy_denom dxy_both dxy_tmp windows.dxy.sort_order_in_scores.bed dxy_sort_order_in_scores.txt
		"""


rule netdiv_windowed:
	input:
		dxy="results_processed/MvF.dxy.bed.txt",
		pi_f_bed="results_processed/F.pi.bed.txt",
		pi_m_bed="results_processed/M.pi.bed.txt",		
		fa=genomefile
	output:
		netdiv_M="results_processed/MvF.netdiv_M.bed.txt",
		netdiv_F="results_processed/MvF.netdiv_F.bed.txt"
	shell:
		"""
		seqtk comp {input.fa} | awk '{{print $1"\\t"$2}}' > genomefile.netdiv.txt
		bedtools makewindows -w {windowsize} -g genomefile.netdiv.txt > windows.netdiv.bed
				
		paste {input.dxy} {input.pi_m_bed} > netdiv_prep		
		awk '{{ if($8!="NA" && $4!="NA") print $1"\\t"$2"\\t"$3"\\t"$4-$8 ; else print $1"\\t"$2"\\t"$3"\\tNA" }}' netdiv_prep > {output.netdiv_M}

		paste {input.dxy} {input.pi_f_bed} > netdiv_prep		
		awk '{{ if($8!="NA" && $4!="NA") print $1"\\t"$2"\\t"$3"\\t"$4-$8 ; else print $1"\\t"$2"\\t"$3"\\tNA" }}' netdiv_prep > {output.netdiv_F}
		
		rm genomefile.netdiv.txt windows.netdiv.bed netdiv_prep
		"""


rule fst_rawstats:
	input:
		popmap={popmapfile},
		gzvcf="results_raw/variant_sites.vcf.gz"
	output:
		"results_raw/fst_raw.txt"		
	shell:
		"""
		cat {input.popmap} | awk '{{if($2==1) print $1}}' > mpop_fst
		cat {input.popmap} | awk '{{if($2==2) print $1}}' > fpop_fst

		vcftools --gzvcf {input.gzvcf} --weir-fst-pop mpop_fst --weir-fst-pop fpop_fst --stdout | awk '{{print $1"\\t"$2"\\t"$2"\\t"$3}}'> {output}
		rm mpop_fst fpop_fst
		"""

rule fst_windowed:
	input:
		fraw="results_raw/fst_raw.txt",
		fa=genomefile
	output:
		"results_processed/MvF.fst.bed.txt"
	shell:
		"""
		seqtk comp {input.fa} | awk '{{print $1"\\t"$2}}' > genomefile.fst.txt
		bedtools makewindows -w {windowsize} -g genomefile.fst.txt > windows.fst.bed
		
		bedtools map -a windows.fst.bed -b {input.fraw} -g genomefile.fst.txt -c 4 -o mean > {output}
		
		rm genomefile.fst.txt windows.fst.bed		
		"""

	
rule GWAS_plink_raw:
	input:
		popmap={popmapfile},
		gzvcf="results_raw/variant_sites.vcf.gz",
		fa=genomefile
	output:
		"results_raw/plink.assoc_results.assoc.raw.txt"
	resources:
		mem_mb=20000,
		cpus=4
	shell:
		"""
		# PLINK has a few problematic behaviours: 
		# - VCF gets converted to temporary plink format files with hard-coded names.. => RACE CONDITON when multiple instances of plink run
		# - it will by default use ALL-1 CPUs => must use argument --threads
		# - By default, PLINK 1.9 tries to reserve half of your system's RAM for its main workspace. => --memory <main workspace size, in MB>
		
		DIR="tmpdir_plink_GWAS"
		if [ -d "$DIR" ]; then
		  rm -r $DIR
		fi
		mkdir tmpdir_plink_GWAS
		cd tmpdir_plink_GWAS

		cat ../{input.popmap} | awk 'NF {{print $1"\\t"$1"\\t"$2}}' > thepheno
		plink --vcf ../{input.gzvcf} --double-id --allow-extra-chr --pheno thepheno --allow-no-sex --assoc --out assoc_results	--threads {resources.cpus} --memory {resources.mem_mb} 
		
		mv assoc_results.assoc ../{output}
		cd ../
		rm -r tmpdir_plink_GWAS	
		"""


rule GWAS_plink_windows:
	input:
		fa=genomefile,
		rawgwas="results_raw/plink.assoc_results.assoc.raw.txt"
	output:
		"results_processed/plink.assoc_results.significant.window.bed.txt",
	shell:
		"""
		cat {input.rawgwas} | tr -s " " "\\t" | awk '{{if($9 <= 0.05) print}}' > plink.assoc_results.significant.txt

		cat plink.assoc_results.significant.txt | cut -f1,3,9 | awk '{{print $1"\\t"$2"\\t"$2}}' > plink.assoc_results.significant.bed

		seqtk comp {input.fa} | awk '{{print $1"\\t"$2}}' > genomefile.assoc.txt
		bedtools makewindows -w {windowsize} -g genomefile.assoc.txt > windows.assoc.bed
		bedtools coverage -counts -a windows.assoc.bed -b plink.assoc_results.significant.bed > {output}

		rm plink.assoc_results.significant.txt plink.assoc_results.significant.bed genomefile.assoc.txt windows.assoc.bed 
		"""

rule LD_plink_raw:
	input:
		gzvcf="results_raw/variant_sites.vcf.gz",
		fa=genomefile
	output:
		"results_raw/ld_clean.sorted.bed.gz"
#		"results_processed/LD_r2_averaged_per_window.bed.txt"
	resources:
		mem_mb=20000,
		cpus=4
	shell:
		"""
		# PLINK has a few problematic behaviours: 
		# - VCF gets converted to temporary plink format files with hard-coded names.. => RACE CONDITON when multiple instances of plink run
		# - it will by default use ALL-1 CPUs => must use argument --threads
		# - By default, PLINK 1.9 tries to reserve half of your system's RAM for its main workspace. => --memory <main workspace size, in MB>
		
		DIR="tmpdir_plink_LD"
		if [ -d "$DIR" ]; then
		  rm -r $DIR
		fi
		mkdir tmpdir_plink_LD
		cd tmpdir_plink_LD

		# apply a moderate MAF filter before calculating LD as rare alleles imply very high LD by definition; they are thus useless for our purpose.
		# vcftools --gzvcf {input.gzvcf} --maf 0.1 --recode --stdout > forplinkld.vcf
		# plink --vcf forplinkld.vcf --double-id --allow-extra-chr --r2 --ld-window-r2 0.0 --ld-window 5 --ld-window-kb 2 
		# rm forplinkld.vcf 
		
		plink --vcf ../{input.gzvcf} --double-id --allow-extra-chr --maf 0.1 --r2 --ld-window-r2 0.0 --ld-window 5 --ld-window-kb 2 --threads {resources.cpus} --memory {resources.mem_mb} 
		
		# --ld-window X	compute LD only for pairs that are at most X SNPs apart (default 10)
		# --ld-window-kb X	compute LD only for pairs that are at most X kb apart (default 1000 kb)
		# --ld-window-r2 X	minimum r2 to report, else omit from output. Default = 0.2

		cat plink.ld | tr -s ' ' '\\t' | cut -f1,2,4,5,7 | tail -n +2 | awk '{{ print $1"\\t"$2"\\t"$2"\\t"$5"\\n"$3"\\t"$4"\\t"$4"\\t"$5  }}' | sort --parallel {resources.cpus} -S 2G -T . -k1,1 -k2,2n | gzip -c > ../{output}
		cd ../
		rm -r tmpdir_plink_LD
		"""
		
rule LD_plink_windows:
	input:
		ldraw="results_raw/ld_clean.sorted.bed.gz",
		fa=genomefile
	output:
		"results_processed/LD_r2_averaged_per_window.bed.txt"
	resources:
		mem_mb=20000,
		cpus=4
	shell:
		"""
		# create windows, sort chroms lexicographically (same as the LD output)
		seqtk comp {input.fa} | awk '{{print $1"\\t"$2}}' > genomefile.ld.txt
		bedtools makewindows -w {windowsize} -g genomefile.ld.txt | sort --parallel {resources.cpus} -S 2G -T . -k1,1 -k2,2n > windows.ld.bed
		
		gunzip --stdout {input.ldraw} > ld_clean.sorted.bed
		
		# get the mean LD per window
		bedtools map -a windows.ld.bed -b ld_clean.sorted.bed -c 4 -o mean > tmpoutperwindow.txt

		# sort chroms back to the same order as the genome file:
		bedtools sort -g genomefile.ld.txt -i tmpoutperwindow.txt > {output}

		rm ld_clean.sorted.bed genomefile.ld.txt tmpoutperwindow.txt windows.ld.bed
		"""
	
	
rule kmerGO:
	input:
		popm={popmapfile},
		SamReMap={samples_reads_map}
	output:
		fkmer="results_raw/F_specific_kmer.fa",
		mkmer="results_raw/M_specific_kmer.fa",
		fcont="results_raw/F_specific.fa.cap.contigs",
		mcont="results_raw/M_specific.fa.cap.contigs"
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
		fkmers="results_raw/F_specific_kmer.fa",
		mkmers="results_raw/M_specific_kmer.fa"
	output:
		fsp="results_processed/kmers.F_specific.per_window.bed.txt",
		msp="results_processed/kmers.M_specific.per_window.bed.txt"
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
	
rule calc_indiv_het:
	input:
		popm={popmapfile},
		gzvcf="results_raw/variant_sites.vcf.gz"
	output:
		mhet="results_raw/het_males.txt",
		fhet="results_raw/het_females.txt"
	shell:
		"""
		vcftools --gzvcf {input.gzvcf} --keep <( cat {input.popm} | awk '{{if($2==1) print $1}}' ) --mac 1 --hardy --stdout | awk -F'[\\t/]' 'NR>1 {{print $1"\\t"$2"\\t"$2"\\t"$4/($3+$4+$5)}}' > {output.mhet}
		vcftools --gzvcf {input.gzvcf} --keep <( cat {input.popm} | awk '{{if($2==2) print $1}}' ) --mac 1 --hardy --stdout | awk -F'[\\t/]' 'NR>1 {{print $1"\\t"$2"\\t"$2"\\t"$4/($3+$4+$5)}}' > {output.fhet}
		"""

rule indiv_het_per_windows:
	input:
		fa=genomefile,
		mhet="results_raw/het_males.txt",
		fhet="results_raw/het_females.txt"		
	output:
		mhetwindows="results_processed/het_males.bed.txt",
		fhetwindows="results_processed/het_females.bed.txt"
	shell:
		"""	
		# create windows, sort chroms lexicographically (same as the LD output)
		seqtk comp {input.fa} | awk '{{print $1"\\t"$2}}' > genomefile.het.txt
		bedtools makewindows -w {windowsize} -g genomefile.het.txt > windows.het.bed

		# get the mean HET per window
		bedtools map -a windows.het.bed -b {input.mhet} -c 4 -o mean > {output.mhetwindows}
		bedtools map -a windows.het.bed -b {input.fhet} -c 4 -o mean > {output.fhetwindows}
		
		rm genomefile.het.txt windows.het.bed		
		"""
	
rule on:
	input:
		gzvcf="results_raw/all.post_filter.vcf.gz",
		popm={popmapfile}
	output:
		"results_raw/gametolog_candidate_alleles.allsites.vcf.gz"
	shell:
		"""
		python scripts/pseudo_phase_gametologs.py {input.popm} {input.gzvcf} | bgzip -c > {output}
		"""
	
rule score_XY_gametolog_divergence:
	input:
		"results_raw/gametolog_candidate_alleles.allsites.vcf.gz"
	output:
		"results_raw/gametolog_candidate_alleles.XY_divergence.rawstats.txt"
	shell:
		"""
		echo -e "XY_Y_like\\tXY_Y_like" > xypopm
		echo -e "XY_X_like\\tXY_X_like" >> xypopm
		
		echo -e "XY_X_like" > xypop
		echo -e "XY_Y_like" >> xypop

		vcftools --gzvcf {input} --keep xypop --recode --stdout | python2.7 scripts/make_denominator_and_numerator_for_dxy.py --popmap xypopm > {output}
		rm xypopm xypop
		"""

rule score_ZW_gametolog_divergence:
	input:
		"results_raw/gametolog_candidate_alleles.allsites.vcf.gz"
	output:
		"results_raw/gametolog_candidate_alleles.ZW_divergence.rawstats.txt"
	shell:
		"""
		echo -e "ZW_W_like\\tZW_W_like" > zwpopm
		echo -e "ZW_Z_like\\tZW_Z_like" >> zwpopm
		
		echo -e "ZW_W_like" > zwpop
		echo -e "ZW_Z_like" >> zwpop

		vcftools --gzvcf {input} --keep zwpop --recode --stdout | python2.7 scripts/make_denominator_and_numerator_for_dxy.py --popmap zwpopm > {output}
		rm zwpopm zwpop
		"""

rule XY_div_windows:
	input:
		raw="results_raw/gametolog_candidate_alleles.XY_divergence.rawstats.txt",
		fa=genomefile
	output:
		"results_processed/gametolog_candidate_alleles.XY_divergence.windows.bed.txt"
	threads: 2
	shell:
		"""
		seqtk comp {input.fa} | awk '{{print $1"\\t"$2}}' > genomefile.xygametologs.txt
		bedtools makewindows -w {windowsize} -g genomefile.xygametologs.txt > windows.xygametologs.bed

		(bedtools map -a windows.xygametologs.bed -b {input.raw} -g genomefile.xygametologs.txt -c 4 -o mean > XY_dxy_num )&
		(bedtools map -a windows.xygametologs.bed -b {input.raw} -g genomefile.xygametologs.txt -c 5 -o mean > XY_dxy_denom )&
		wait
		
		paste XY_dxy_num XY_dxy_denom > XY_dxy_both
		
		awk '{{ if($8>0) print $1"\\t"$2"\\t"$3"\\t"$4/$8 ; else print $1"\\t"$2"\\t"$3"\\tNA" }}' XY_dxy_both > {output}
		
		rm genomefile.xygametologs.txt windows.xygametologs.bed XY_dxy_denom XY_dxy_num XY_dxy_both		
		"""

rule ZW_div_windows:
	input:
		raw="results_raw/gametolog_candidate_alleles.ZW_divergence.rawstats.txt",
		fa=genomefile
	output:
		"results_processed/gametolog_candidate_alleles.ZW_divergence.windows.bed.txt"
	threads: 2
	shell:
		"""
		seqtk comp {input.fa} | awk '{{print $1"\\t"$2}}' > genomefile.zwgametologs.txt
		bedtools makewindows -w {windowsize} -g genomefile.zwgametologs.txt > windows.zwgametologs.bed

		( bedtools map -a windows.zwgametologs.bed -b {input.raw} -g genomefile.zwgametologs.txt -c 4 -o mean > ZW_dxy_num )&
		( bedtools map -a windows.zwgametologs.bed -b {input.raw} -g genomefile.zwgametologs.txt -c 5 -o mean > ZW_dxy_denom )&
		wait
				
		paste ZW_dxy_num ZW_dxy_denom > ZW_dxy_both
		
		awk '{{ if($8>0) print $1"\\t"$2"\\t"$3"\\t"$4/$8 ; else print $1"\\t"$2"\\t"$3"\\tNA" }}' ZW_dxy_both > {output}
		
		rm genomefile.zwgametologs.txt windows.zwgametologs.bed ZW_dxy_denom ZW_dxy_num ZW_dxy_both		
		"""


rule plot_all:
	input:
		regions=regions_for_plot_bed,
		a="results_processed/M_F_norm_coverage_ratio_log2.bed.txt",
		b="results_processed/kmers.F_specific.per_window.bed.txt",
		c="results_processed/kmers.M_specific.per_window.bed.txt",
		d="results_processed/F.pi.bed.txt",
		e="results_processed/M.pi.bed.txt",
		f="results_processed/MvF.fst.bed.txt",
		g="results_processed/MvF.dxy.bed.txt",
		h="results_processed/LD_r2_averaged_per_window.bed.txt",
		i="results_processed/plink.assoc_results.significant.window.bed.txt",
		j="results_processed/gametolog_candidate_alleles.ZW_divergence.windows.bed.txt",
		k="results_processed/gametolog_candidate_alleles.XY_divergence.windows.bed.txt",
		l="results_processed/het_males.bed.txt",
		m="results_processed/het_females.bed.txt",
		n="results_processed/MvF.netdiv_F.bed.txt",
		o="results_processed/MvF.netdiv_M.bed.txt"
	output:
		stats="results_processed/allstats.txt",
		regionstats="results_processed/region.stats.txt",
		plot="results_processed/allstats.plots.pdf",
		regionplot="results_processed/region.stats.plots.pdf"
	shell:
		"""
		cat {input.a} > {output.stats}

		for i in {input.b} {input.c} {input.d} {input.e} {input.f} {input.g} {input.h} {input.i} {input.j} {input.k} {input.l} {input.m} {input.n} {input.o} ; do
			paste {output.stats} <(cut -f4 $i ) > tmpst
			mv tmpst {output.stats}
		done
		
		bedtools intersect -wa -a {output.stats} -b {input.regions} > {output.regionstats}
		
		Rscript scripts/plot_allstats.R {output.stats} {output.plot}
		Rscript scripts/plot_allstats.R {output.regionstats} {output.regionplot}	
		"""


	
