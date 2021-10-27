#!/usr/local/bin/python
# Python 3

# Mathias Scharmann

"""
for pooled data, we define pi based on the count of alleles in the READS instead of the SAMPLES/CHROMOSOMES:

pi = ( count_p_in_reads *( count_all_reads - count_p_in_reads) ) / ( (count_all_reads*(count_all_reads-1))/2.0 )

"""


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
	gt_field = fields[9]
	gts = gt_field.split(":")[0].replace("/","")
	try:
		alleles = set(gts)
		alleles.discard(".")
		n_alleles = len(alleles)
	except TypeError:
		n_alleles = 3 # dummy
	if n_alleles in set([1,2]):	# only bi-allelic or fixed sites, NOT MISSING SITES: sites with more than 2 alleles are NOT utilised!!
		gt_fields = fields[9].split(":")
		count_all_reads = float( gt_fields[3] )
		n_pairs = (count_all_reads*(count_all_reads-1.0))/2.0 ## this is the site-specific denominator
		count_p_in_reads = float( gt_fields[1] )
		countsproduct = ( count_p_in_reads*( count_all_reads - count_p_in_reads) )
		# pi_combin will return the same result, but takes MUCH longer!
		# pi_combin = alldiffs( [x for x in gts if x != "."] ) / scipy.special.binom(nt, 2)
		outl += [str(countsproduct), str(n_pairs)]
		sys.stdout.write("\t".join(outl) + "\n")
						

