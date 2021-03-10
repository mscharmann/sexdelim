## pseudo_phase_gametologs.py

import gzip, sys

def read_popmap (infile):
	
	popdict = {}
	with open(infile, "r") as INF:
		for line in INF:
			if len( line ) > 2:
				fields = line.strip("\n").split("\t")
				try:
					popdict[ fields[1] ].append( fields[0] )
				except KeyError:
					popdict[ fields[1] ] = [ fields[0] ]
	return popdict				


def get_major_allele (gts_p1, gts_p2):
	nomiss = list( (gts_p1 + gts_p2).replace(".", "") )
	return max(set(nomiss), key = nomiss.count) 
	


popmapfile = sys.argv[1]
gzvcffile = sys.argv[2]

popdict = read_popmap(popmapfile)

vcf_infile = gzip.open(gzvcffile, "rt")
for line in vcf_infile: # discard the header!
	if line.startswith("#CHROM"):
		header_line = line.lstrip("#").strip("\n").split("\t")
		samples = sorted(header_line[9:])
		popdict_vcf_idx = {}
		for pop in popdict.keys():
			for sample in  popdict[pop]:
				vcf_idx = header_line.index(sample)
				try:
					popdict_vcf_idx[pop].append(vcf_idx)
				except KeyError:
					popdict_vcf_idx[pop] = [vcf_idx]
		pops = sorted(popdict.keys())
		break

vcf_header = ["##fileformat=VCFv4.2",
"""##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">""",
"\t".join(["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","XY_Y_like","XY_X_like","ZW_W_like","ZW_Z_like"])]

# vcf_outfile = open("gametolog_candidate_alleles.allsites.vcf", "w")
# vcf_outfile.write("\n".join(vcf_header) + "\n")

sys.stdout.write( "\n".join(vcf_header) + "\n" )


#print(popdict_vcf_idx)


for line in vcf_infile :
	fields = line.strip("\n").split("\t")
	chrom = fields[0]
	pos = fields[1]
	refallele = fields[3]
	altallele = fields[4]
	gt_fields = fields[9:]
	gts_p1 = "".join( [x.split(":")[0] for x in [fields[i] for i in popdict_vcf_idx[pops[0]]]] ).replace("/","")		
	gts_p2 = "".join( [x.split(":")[0] for x in [fields[i] for i in popdict_vcf_idx[pops[1]]]] ).replace("/","")
	pres_M = float((len(gts_p1)-gts_p1.count(".")))
	pres_F = float((len(gts_p2)-gts_p2.count(".")))
	
	major_allele = get_major_allele (gts_p1, gts_p2)
	hom_maj = major_allele + "/" + major_allele
	
	try:
		alleles = set(gts_p1)
		alleles = alleles.union( set(gts_p2) )
		alleles.discard(".")
		n_alleles = len(alleles)
	except TypeError: 
		alleles = set()
		n_alleles = 3 # dummy

	# continue only for bi-allelic or fixed sites, NOT MISSING SITES
	if n_alleles == 1: # fixed site, may or may not have missing data
		alleles = list(alleles)
		if pres_M == 0:
			if pres_F == 0:
				# both absent; in a correctly fitlered input VCF this case shoud never occur.
				sys.stdout.write( "\t".join([chrom, pos, ".", refallele, altallele, "60", ".", ".", "GT", "./.", "./.", "./.", "./." ]) + "\n" )
			else:
				# F present, M absent: can not have Y-like or Z-like; X-like and W-like are fine!
				# they are assigned the major allele
				XY_Y_like = "./."
				ZW_Z_like = "./."
				XY_X_like = hom_maj
				ZW_W_like = hom_maj
				sys.stdout.write( "\t".join([chrom, pos, ".", refallele, altallele, "60", ".", ".", "GT", XY_Y_like, XY_X_like, ZW_W_like, ZW_Z_like ]) + "\n" )
		elif pres_F == 0:
			# F absent, M present: : can not have X-like or W-like; Y-like and Z-like are fine!
			# they are assigned the major allele
			XY_Y_like = hom_maj
			ZW_Z_like = hom_maj
			XY_X_like = "./."
			ZW_W_like = "./."
			sys.stdout.write( "\t".join([chrom, pos, ".", refallele, altallele, "60", ".", ".", "GT", XY_Y_like, XY_X_like, ZW_W_like, ZW_Z_like ]) + "\n" )
		else:
			# both present	
			sys.stdout.write( "\t".join([chrom, pos, ".", refallele, altallele, "60", ".", ".", "GT", hom_maj, hom_maj, hom_maj, hom_maj ]) + "\n" )

	elif n_alleles == 2: #  bi-allelic site; may or may not have missing data
		alleles = list(alleles)
		if pres_M == 0:
			if pres_F == 0:
				# both absent; in a correctly fitlered input VCF this case shoud never occur.
				XY_Y_like = "./."
				ZW_Z_like = "./."
				XY_X_like = "./."
				ZW_W_like = "./."
			else:
				# F present, M absent: can not have Y-like or Z-like; X-like and W-like are fine
				# they are assigned the major allele
				XY_Y_like = "./."
				ZW_Z_like = "./."
				XY_X_like = hom_maj
				ZW_W_like = hom_maj

		elif pres_F == 0:
			# F absent, M present: : can not have X-like or W-like; Y-like and Z-like are fine; 
			# they are assigned the major allele
			XY_X_like = "./."
			ZW_W_like = "./."
			XY_Y_like = hom_maj
			ZW_Z_like = hom_maj

		else:
			# both present; now investigate heterozygosity and frequency in either sex	
			freq_0_p1 = gts_p1.count(alleles[0]) / float((len(gts_p1)-gts_p1.count(".")))
			freq_1_p1 = gts_p1.count(alleles[1]) / float((len(gts_p1)-gts_p1.count(".")))
			freq_0_p2 = gts_p2.count(alleles[0]) / float((len(gts_p2)-gts_p2.count(".")))
			freq_1_p2 = gts_p2.count(alleles[1]) / float((len(gts_p2)-gts_p2.count(".")))
			gts_p1 = [y.split("/") for y in [x.split(":")[0] for x in [fields[i] for i in popdict_vcf_idx["1"]]]]
			gts_p2 = [y.split("/") for y in [x.split(":")[0] for x in [fields[i] for i in popdict_vcf_idx["2"]]]]
			het_males = len([x for x in gts_p1 if x[0] != x[1]]) / float( len(gts_p1) )
			het_females = len([x for x in gts_p2 if x[0] != x[1]]) / float( len(gts_p2) )
			if het_males == 1.0:
				if freq_0_p2 == 1:
					#print (alleles[1], "Y-like")
					XY_Y_like = alleles[1] + "/" + alleles[1]
					XY_X_like = alleles[0] + "/" + alleles[0]
					ZW_W_like = hom_maj
					ZW_Z_like = hom_maj
				elif freq_1_p2 == 1:
					#print (alleles[0], "Y-like")
					XY_Y_like = alleles[0] + "/" + alleles[0]
					XY_X_like = alleles[1] + "/" + alleles[1] 
					ZW_W_like = hom_maj
					ZW_Z_like = hom_maj
				else:
					# not gametolog-like; therefore all categories get the major allele 
					XY_X_like = hom_maj
					ZW_W_like = hom_maj
					XY_Y_like = hom_maj
					ZW_Z_like = hom_maj

			elif het_females == 1.0:
				if freq_0_p1 == 1:
					XY_Y_like = hom_maj
					XY_X_like = hom_maj
					ZW_W_like = alleles[1] + "/" + alleles[1]
					ZW_Z_like = alleles[0] + "/" + alleles[0]
				elif freq_1_p1 == 1:
					XY_Y_like = hom_maj
					XY_X_like = hom_maj
					ZW_W_like = alleles[0] + "/" + alleles[0]
					ZW_Z_like = alleles[1] + "/" + alleles[1]
				else:
					# not gametolog-like; therefore all categories get the major allele 
					XY_X_like = hom_maj
					ZW_W_like = hom_maj
					XY_Y_like = hom_maj
					ZW_Z_like = hom_maj

			else:
				# not gametolog-like; therefore all categories get the major allele
				XY_X_like = hom_maj
				ZW_W_like = hom_maj
				XY_Y_like = hom_maj
				ZW_Z_like = hom_maj
					
		sys.stdout.write( "\t".join([chrom, pos, ".", refallele, altallele, "60", ".", ".", "GT", XY_Y_like, XY_X_like, ZW_W_like, ZW_Z_like ]) + "\n" )

