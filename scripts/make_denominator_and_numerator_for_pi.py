#!/usr/local/bin/python
# Python 3

# Mathias Scharmann




#import scipy.special
import sys

########################## HEAD
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
	print (len(seen)) # this is correct!
	return diff	

################################## MAIN

for line in sys.stdin:
	if line.startswith("#"):
		continue
	if len(line) < 2: # empty lines or so
		continue
	fields = line.strip("\n").split("\t")
	outl = [fields[0],str(int(fields[1])-1), fields[1]] 
	# convert from 1-based, closed [start, end] Variant Call Format v4 (VCF)
	# to sorted, 0-based, half-open [start-1, end) extended BED data
	gt_fields = fields[9:]
	gts = "".join( [x.split(":")[0] for x in gt_fields] ).replace("/","")		
	try:
		alleles = set(gts)
		alleles.discard(".")
		n_alleles = len(alleles)
	except TypeError:
		n_alleles = 3 # dummy
	if n_alleles in set([1,2]):	# only bi-allelic or fixed sites, NOT MISSING SITES
		totlen = len(gts)
		nt = totlen - gts.count(".")
		n_pairs = (nt*(nt-1))/2 ## this is the site-specific denominator
		count_p = gts.count(alleles.pop())  # count of a random allele , incl fixed sites!
		# print count_p
		countsproduct = ( count_p*( nt- count_p) )
		# pi_combin will return the same result, but takes MUCH longer!
		# pi_combin = alldiffs( [x for x in gts if x != "."] ) / scipy.special.binom(nt, 2)
		outl += [str(countsproduct), str(n_pairs)]
		sys.stdout.write("\t".join(outl) + "\n")
						

