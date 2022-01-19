#!/usr/local/bin/python
# Python 3
# 
# 
# Mathias Scharmann


# usage example
# 

"""
dxy

	Nei & Li (1979): unnumbered eq., between eqs. 24 and 25 
	Nei 1987 eq 10.20 
	for a single biallelic SNP (1,2) in two pops it simplifies to: dxy = pop1_freq1 * ( 1.0 - pop2_freq1 ) + pop2_freq1 * ( 1.0 - pop1_freq1 )


- this script will only work for a poolSNP VCF with exactly TWO populations in it.

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
	if line.startswith("#CHROM"):
		header_line = line.lstrip("#").strip("\n").split("\t")
		samples = sorted(header_line[9:])
	
		continue		
	elif line.startswith("##"):
		continue
	elif len(line) < 2: # empty lines or so
		continue
	fields = line.strip("\n").split("\t")
	outl = [fields[0],str(int(fields[1])-1), fields[1]]
	# convert from 1-based, closed [start, end] Variant Call Format v4 (VCF)
	# to sorted, 0-based, half-open [start-1, end) extended BED data

	gts_p1 = fields[9].split(":")[0].replace("/","")		
	gts_p2 = fields[10].split(":")[0].replace("/","")
	try:
		alleles = set(gts_p1)
		alleles = alleles.union( set(gts_p2) )
		alleles.discard(".")
		n_alleles = len(alleles)
	except TypeError: 
		alleles = set()
		n_alleles = 3 # dummy
	if n_alleles in set([1,2]):	# only bi-allelic or fixed sites, NOT MISSING SITES: sites with more than 2 alleles are NOT utilised!!
		gt_fields_p1 = fields[9].split(":")
		gt_fields_p2 = fields[10].split(":")
		count_all_reads_p1 = float( gt_fields_p1[3] )
		count_all_reads_p2 = float( gt_fields_p2[3] )
		n_pairs = count_all_reads_p1 * count_all_reads_p2 ## this is the site-specific denominator; susbtract missing genotypes
		try:
			count_p_p1 = float( gt_fields_p1[1])
			count_q_p1 = float( gt_fields_p1[1])
		except ValueError: # gt_fields_xx is of the form: '1/2:0:16,2:18:0.89,0.11'. This means we had 2 ALT alleles, while the REF allele did NOT occur in the populations. It is a bi-allelic site nevertheless and will be used. 
			count_p_p1 = float( gt_fields_p1[1].split(",")[0] )
			count_q_p1 = float( gt_fields_p1[2].split(",")[1] )
		try:
			count_p_p2 = float( gt_fields_p2[1])
			count_q_p2 = float( gt_fields_p2[2])
		except ValueError:
			count_p_p2 = float( gt_fields_p2[1].split(",")[0] )
			count_q_p2 = float( gt_fields_p2[2].split(",")[1] )
			
		countsproduct_sum = ( count_p_p1*count_q_p2 ) + ( count_q_p1*count_p_p2 ) 
		outl += [str(countsproduct_sum), str(n_pairs)]
		sys.stdout.write("\t".join(outl) + "\n")
						

