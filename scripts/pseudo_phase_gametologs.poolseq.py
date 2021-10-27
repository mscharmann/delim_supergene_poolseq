## pseudo_phase_gametologs.py

"""
FREQ in a poolSNP vcf is ALWAYS the frequency of the ALT allele.

"""


import gzip, sys


def get_major_allele (ALT_freq_p1, ALT_freq_p2):
	if ALT_freq_p1 == ".":
		if float(ALT_freq_p2) > 0.5:
			return "ALT"
		else:
			return "REF"
	if ALT_freq_p2 == ".":
		if float(ALT_freq_p1) > 0.5:
			return "ALT"
		else:
			return "REF"
	avg_freq = (float(ALT_freq_p1) + float(ALT_freq_p2))/2.0
	if avg_freq < 0.5:
		major = "REF"
	else:
		major = "ALT"
	return major
	

vcf_header = ["##fileformat=VCFv4.2",
"""##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">""",
"\t".join(["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","XY_Y_like","XY_X_like","ZW_W_like","ZW_Z_like"])]

# vcf_outfile = open("gametolog_candidate_alleles.allsites.vcf", "w")
# vcf_outfile.write("\n".join(vcf_header) + "\n")

sys.stdout.write( "\n".join(vcf_header) + "\n" )


#print(popdict_vcf_idx)


for line in sys.stdin:
	if line.startswith("#"):
		if line.startswith("#CHROM"):
			header_line = line.lstrip("#").strip("\n").split("\t")
			samples = sorted(header_line[9:])
			continue
		continue
	if len(line) < 2: # empty lines or so
		continue
	fields = line.strip("\n").split("\t")
	chrom = fields[0]
	pos = fields[1]
	refallele = fields[3]
	altallele = fields[4]
	
	# p1 are the 'males', p2 are the 'females'
	
	gts_p1 = fields[9].split(":")[0].replace("/","")		
	gts_p2 = fields[10].split(":")[0].replace("/","")

	gt_fields = fields[9:]

	pres_M = float((len(gts_p1)-gts_p1.count(".")))
	pres_F = float((len(gts_p2)-gts_p2.count(".")))
	
	major_allele = get_major_allele (fields[9].split(":")[-1] , fields[10].split(":")[-1])
	if major_allele == "REF":
		hom_maj = "0/0"
		hom_min = "1/1"
	else:
		hom_maj = "1/1"
		hom_min = "0/0"
			
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
				#continue
			else:
				# F present, M absent: can not have Y-like or Z-like; X-like and W-like are fine!
				# they are assigned the major allele
				XY_Y_like = "./."
				ZW_Z_like = "./."
				XY_X_like = hom_maj
				ZW_W_like = hom_maj
				sys.stdout.write( "\t".join([chrom, pos, ".", refallele, altallele, "60", ".", ".", "GT", XY_Y_like, XY_X_like, ZW_W_like, ZW_Z_like ]) + "\n" )
				#continue
		elif pres_F == 0:
			# F absent, M present: : can not have X-like or W-like; Y-like and Z-like are fine!
			# they are assigned the major allele
			XY_Y_like = hom_maj
			ZW_Z_like = hom_maj
			XY_X_like = "./."
			ZW_W_like = "./."
			sys.stdout.write( "\t".join([chrom, pos, ".", refallele, altallele, "60", ".", ".", "GT", XY_Y_like, XY_X_like, ZW_W_like, ZW_Z_like ]) + "\n" )
			#continue
		else:
			# both present	
			sys.stdout.write( "\t".join([chrom, pos, ".", refallele, altallele, "60", ".", ".", "GT", hom_maj, hom_maj, hom_maj, hom_maj ]) + "\n" )
			#continue

			
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
			# both present; now investigate frequency in either sex	
			not_gametolog_like = False
			freq_0_M = float( fields[9].split(":")[4] )
			freq_0_F = float( fields[10].split(":")[4] )

			if freq_0_M > 0.35 and freq_0_M < 0.65: # around freq 0.5 in M, but fixed in F -> X-Y patterned. because site present in both M and F, ALT must be the Y-like allele.
				if freq_0_F == 0 or freq_0_F == 1.0:
					XY_Y_like = hom_min
					XY_X_like = hom_maj
					ZW_W_like = hom_maj
					ZW_Z_like = hom_maj
					#print(fields, "X-Y like")
				else:
					not_gametolog_like = True

			elif freq_0_F > 0.35 and freq_0_F < 0.65: # around freq 0.5 in F, but fixed in M -> Z-W patterned. because site present in both M and F, ALT must be the W-like allele.
				if freq_0_M == 0 or freq_0_M == 1.0:
					XY_Y_like = hom_maj
					XY_X_like = hom_maj
					ZW_W_like = hom_min
					ZW_Z_like = hom_maj
					#print(fields, "Z-W like")
				else:
					not_gametolog_like = True
			else:
				not_gametolog_like = True			
			if not_gametolog_like:
				# not gametolog-like; therefore all categories get the major allele
				XY_X_like = hom_maj
				ZW_W_like = hom_maj
				XY_Y_like = hom_maj
				ZW_Z_like = hom_maj


		sys.stdout.write( "\t".join([chrom, pos, ".", refallele, altallele, "60", ".", ".", "GT", XY_Y_like, XY_X_like, ZW_W_like, ZW_Z_like ]) + "\n" )













