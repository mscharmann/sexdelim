#!/usr/local/bin/python
# Python 2
# 
# 
# Mathias Scharmann


# usage example
# 

"""
dxy
"""


#import scipy.special
import sys
import argparse

########################## HEAD

# parses command line arguments
def get_commandline_arguments ():
	
	parser = argparse.ArgumentParser()
	
	parser.add_argument("--popmap", required=True, help="name/path of popmap file", metavar="FILE")
	
	args = parser.parse_args()

	return args

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


def alldiffs (inlist): 
	
	seen = set()
	diff = 0
	for i in range(len(inlist)):
		for j in range(len(inlist)):
			if i != j:
				pairstring = "-".join(sorted([str(x) for x in [i,j]]))
				if not pairstring in seen:
					if inlist[i] != inlist[j]:
						diff += 1
					seen.add( pairstring )	
	print len(seen) # this is correct!
	return diff	


################################## MAIN

args = get_commandline_arguments ()

popdict= read_popmap(args.popmap)


for line in sys.stdin:
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
		
		continue		
	elif line.startswith("##"):
		continue
	elif len(line) < 2: # empty lines or so
		continue
	fields = line.strip("\n").split("\t")
	outl = fields[:2]
	gt_fields = fields[9:]
	gts_p1 = "".join( [x.split(":")[0] for x in [fields[i] for i in popdict_vcf_idx[pops[0]]]] ).replace("/","")		
	gts_p2 = "".join( [x.split(":")[0] for x in [fields[i] for i in popdict_vcf_idx[pops[1]]]] ).replace("/","")
	try:
		alleles = set(gts_p1)
		alleles.add( gts_p2 )
		alleles.discard(".")
		n_alleles = len(alleles)
	except TypeError:
		n_alleles = 3 # dummy
	if n_alleles <= 2:	# only bi-allelic or fixed sites
		n_pairs = len(gts_p1) * len(gts_p2) ## this is the site-specific denominator
		count_p1 = gts_p1.count(alleles.pop())  # count of a random allele , incl fixed sites!
		count_p2 = gts_p2.count(alleles.pop())  # count THE OTHER ALLELE
		# print count_p
		countsproduct = ( count_p1*count_p2)
		# pi_combin will return the same result, but takes MUCH longer!
		# pi_combin = alldiffs( [x for x in gts if x != "."] ) / scipy.special.binom(nt, 2)
		outl += [str(countsproduct), str(n_pairs)]
		sys.stdout.write("\t".join(outl) + "\n")
						

