## pseudo_phase_gametologs.py

import random, gzip, sys, numpy as np

mfreqfile = sys.argv[1]
ffreqfile = sys.argv[2]
gzvcffile = sys.argv[3]

m_vcf_freqs_file = gzip.open(mfreqfile, "rt")
f_vcf_freqs_file = gzip.open(ffreqfile, "rt")

m_vcf_freqs_file.readline()
f_vcf_freqs_file.readline()

vcf_infile = gzip.open(gzvcffile, "rt")
for line in vcf_infile: # discard the header!
	if line.startswith("#CHROM"):
		break

vcf_header = ["##fileformat=VCFv4.2",
"""##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">""",
"\t".join(["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","XY_Y_like","XY_X_like","ZW_W_like","ZW_Z_like"])]

# vcf_outfile = open("gametolog_candidate_alleles.allsites.vcf", "w")
# vcf_outfile.write("\n".join(vcf_header) + "\n")

sys.stdout.write( "\n".join(vcf_header) + "\n" )

for vcfline, mline, fline in zip(vcf_infile, m_vcf_freqs_file, f_vcf_freqs_file) :
	fields = vcfline.strip("\n").split("\t")
	chrom = fields[0]
	pos = fields[1]
	refallele = fields[3]
	altallele = fields[4]
	mfields = mline.strip("\n").split("\t")
	ffields = fline.strip("\n").split("\t")
	mfreq_p = float(mfields[4].split(":")[1])
	ffreq_p = float(ffields[4].split(":")[1])
	if altallele == ".": # fixed sites 
		if np.isnan(mfreq_p):
			if np.isnan(ffreq_p):
				# both absent; in a correctly fitlered input VCF this case shoud never occur.
				sys.stdout.write( "\t".join([chrom, pos, ".", refallele, altallele, "60", ".", ".", "GT", "./.", "./.", "./.", "./." ]) + "\n" )
			else:
				# F present, M absent: can not have Y-like or Z-like; X-like and W-like are fine!
				XY_Y_like = "./."
				ZW_Z_like = "./."
				XY_X_like = "0/0"
				ZW_W_like = "0/0"
				sys.stdout.write( "\t".join([chrom, pos, ".", refallele, altallele, "60", ".", ".", "GT", XY_Y_like, XY_X_like, ZW_W_like, ZW_Z_like ]) + "\n" )
		elif np.isnan(ffreq_p):
			# F absent, M present: : can not have X-like or W-like; Y-like and Z-like are fine!
			XY_Y_like = "0/0"
			ZW_Z_like = "0/0"
			XY_X_like = "./."
			ZW_W_like = "./."
			sys.stdout.write( "\t".join([chrom, pos, ".", refallele, altallele, "60", ".", ".", "GT", XY_Y_like, XY_X_like, ZW_W_like, ZW_Z_like ]) + "\n" )
		else:
			# both present	
			sys.stdout.write( "\t".join([chrom, pos, ".", refallele, altallele, "60", ".", ".", "GT", "0/0", "0/0", "0/0", "0/0" ]) + "\n" )
	else: # variable sites
		mfreq_q = float(mfields[5].split(":")[1])
		ffreq_q = float(ffields[5].split(":")[1])
		m_p = mfields[4].split(":")[0]
		m_q = mfields[5].split(":")[0]
		f_p = ffields[4].split(":")[0]
		f_q = ffields[5].split(":")[0]
		if np.isnan(mfreq_p):
			if np.isnan(ffreq_p):
				# both absent; in a correctly fitlered input VCF this case shoud never occur.
				sys.stdout.write( "\t".join([chrom, pos, ".", refallele, altallele, "60", ".", ".", "GT", "./.", "./.", "./.", "./." ]) + "\n" )
			else:
				# F present, M absent: can not have Y-like or Z-like; X-like and W-like are fine!
				XY_Y_like = "./."
				ZW_Z_like = "./."
				if ffreq_p >= ffreq_q:
					XY_X_like = f_p
				else:
					XY_X_like = f_q
				if ffreq_p >= 0.4 and mfreq_p == 0.0:
					ZW_W_like = f_p
				elif ffreq_q >= 0.4 and mfreq_q == 0.0:
					ZW_W_like = f_q
				else:
					ZW_W_like = random.choice([f_p,f_q])
				if XY_X_like == refallele:
					s2 = "0/0"
				else:
					s2 = "1/1"
				if ZW_W_like == refallele:
					s3 = "0/0"
				else:
					s3 = "1/1"
				sys.stdout.write( "\t".join([chrom, pos, ".", refallele, altallele, "60", ".", ".", "GT", XY_Y_like, s2, s3, ZW_Z_like ]) + "\n" )
		elif np.isnan(ffreq_p):
			# F absent, M present: : can not have X-like or W-like; Y-like and Z-like are fine!
			XY_X_like = "./."
			ZW_W_like = "./."
			if mfreq_p >= 0.4 and ffreq_p == 0.0:
				XY_Y_like = m_p
			elif mfreq_q >= 0.4 and ffreq_q == 0.0:
				XY_Y_like = m_q
			else:
				XY_Y_like = random.choice([m_p,m_q])
			if mfreq_p >= mfreq_q:
				ZW_Z_like = m_p
			else:
				ZW_Z_like = m_q
			if XY_Y_like == refallele:
				s1 = "0/0"
			else:
				s1 = "1/1"
			if ZW_Z_like == refallele:
				s4 = "0/0"
			else:
				s4 = "1/1"
			sys.stdout.write( "\t".join([chrom, pos, ".", refallele, altallele, "60", ".", ".", "GT", s1, XY_X_like, ZW_W_like, s4 ]) + "\n" )
		else:
			# both present	
			if ffreq_p >= ffreq_q:
				XY_X_like = f_p
			else:
				XY_X_like = f_q
			if mfreq_p >= 0.4 and ffreq_p == 0.0:
				XY_Y_like = m_p
			elif mfreq_q >= 0.4 and ffreq_q == 0.0:
				XY_Y_like = m_q
			else:
				XY_Y_like = random.choice([m_p,m_q])
			# print (XY_Y_like, XY_X_like)
			if ffreq_p >= 0.4 and mfreq_p == 0.0:
				ZW_W_like = f_p
			elif ffreq_q >= 0.4 and mfreq_q == 0.0:
				ZW_W_like = f_q
			else:
				ZW_W_like = random.choice([f_p,f_q])
			if mfreq_p >= mfreq_q:
				ZW_Z_like = m_p
			else:
				ZW_Z_like = m_q
			# print (ZW_W_like, ZW_Z_like)
			# order of "samples" in output vcf: "XY_Y_like","XY_X_like","ZW_W_like","ZW_Z_like"
			if XY_Y_like == refallele:
				s1 = "0/0"
			else:
				s1 = "1/1"
			if XY_X_like == refallele:
				s2 = "0/0"
			else:
				s2 = "1/1"
			if ZW_W_like == refallele:
				s3 = "0/0"
			else:
				s3 = "1/1"
			if ZW_Z_like == refallele:
				s4 = "0/0"
			else:
				s4 = "1/1"
			# vcf_outfile.write( "\t".join([chrom, pos, ".", refallele, altallele, "60", ".", ".", "GT", s1, s2, s3, s4 ]) + "\n")
			sys.stdout.write( "\t".join([chrom, pos, ".", refallele, altallele, "60", ".", ".", "GT", s1, s2, s3, s4 ]) + "\n" )


# vcf_outfile.close()
vcf_infile.close()
m_vcf_freqs_file.close()
f_vcf_freqs_file.close()

