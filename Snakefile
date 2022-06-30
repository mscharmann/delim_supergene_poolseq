configfile: "config.yaml"

windowsize = int( config["windowsize"] )
genomefile = config["genomefile"]
samples_units_fqs_map = config["samples_units_fqs_map"]
regions_for_plot_bed = config["regions_for_plot_bed"]

avg_samples_per_pool = int( config["avg_samples_per_pool"] )
kmer_covs = [int(0.5*avg_samples_per_pool), 1*avg_samples_per_pool, 2*avg_samples_per_pool, 5*avg_samples_per_pool, 10*avg_samples_per_pool, 100*avg_samples_per_pool]

import pandas as pd

samples_units_fqs = pd.read_table(samples_units_fqs_map, dtype=str).set_index(
	["sample", "unit", "fq1", "fq2"], drop=False)

SAMPLES = list( set(samples_units_fqs["sample"]) ) 
print(samples_units_fqs)


def get_sample_bams(wildcards):
	"""Get all aligned reads of given sample."""
	return expand(
		"mapped_reads_per_unit/{sample}-{unit}.sorted.bam",
		sample=wildcards.sample,
		unit=samples_units_fqs.loc[wildcards.sample].unit,
	)


def get_fastq_sample_unit(wildcards):
	"""Get fastq files of given sample-unit."""
	fastqs = samples_units_fqs.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]]
	if fastqs.fq2.isnull().values.any():
		return [ fastqs.fq1.item() ]
	return [ fastqs.fq1.item(), fastqs.fq2.item() ]


def get_fastq_sample_ALL(wildcards):
	"""Get list of fastq files of given sample including ALL units."""
	fastqs_pd = samples_units_fqs.loc[(wildcards.sample), ["fq1", "fq2"]]
	fastqs = set( fastqs_pd.fq1.tolist() + fastqs_pd.fq2.tolist() )
	fastqs_clean = list( {x for x in fastqs if pd.notna(x)} )
	return fastqs_clean



rule all:
	input:
		"results_raw/poolsnp.vcf.gz",
		expand("results_processed/window_cov.{sample}.txt", sample=SAMPLES),
		"results_processed/norm_coverage_ratio_log2.bed.txt",
		"results_raw/poolsnp.vcf.gz",
		expand("results_processed/pi_per_window.{sample}.bed.txt", sample=SAMPLES),
		"results_processed/dxy_per_window.bed.txt",
		expand("results_processed/netDiv_per_window.{sample}.bed.txt", sample=SAMPLES),
		"results_processed/Fst_Hudson_per_window.bed.txt",
		"results_processed/gametolog_candidate_alleles.ZW_divergence.windows.bed.txt",
		"results_processed/gametolog_candidate_alleles.XY_divergence.windows.bed.txt",
		"results_raw/kmers_diff_ci_evaluation.txt",
		expand("results_processed/kmers.{sample}_specific.per_window.bed.txt", sample=SAMPLES),
		"results_processed/allstats.txt",
		"results_processed/region.stats.txt",
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
		reads=get_fastq_sample_unit
	output:
		temp("mapped_reads_per_unit/{sample}-{unit}.bam")
	threads: 24
	run:
		if len(input.reads) == 2: # paired-end!
			shell("""
				# filtering alignments to be primary (i.e. each read only once; if multiple locations equally possible than a random one is chosen): 
				# -F 256 == -F 0x0100 == NOT not primary alignment
				# filtering alignments to be NOT supplementary (supplementary: sections of read map to discontinuous coordinates, e.g. across an inversion breakpoint..):
				# -F 2048 == -F 0x800 == NOT supplementary alignment
				# sum of the bit flags: 2304 => filters against BOTH non-primary and supplementary alignments; verified with samtools flagstat
				# filtering alignments to be "properly paired": -f 2
				# filtering against multi-mapping alignments (which have MAPQ=0): --min-MQ 1
				bwa mem -t {threads} -a {input.fa} {input.reads[0]} {input.reads[1]} -R "@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tPL:Illumina" | samtools view -F 2304 -f 2 --min-MQ 1 -b -@ 2 - > {output} 
				""")
		else: # single-end
			shell("""
				# filtering alignments to be primary (i.e. each read only once; if multiple locations equally possible than a random one is chosen): 
				# -F 256 == -F 0x0100 == NOT not primary alignment
				# filtering alignments to be NOT supplementary (supplementary: sections of read map to discontinuous coordinates, e.g. across an inversion breakpoint..):
				# -F 2048 == -F 0x800 == NOT supplementary alignment
				# -F 4 read unmapped (0x4)
				# sum of the bit flags: 2308 => filters against non-primary and supplementary alignments and unmapped
				# filtering against multi-mapping alignments (which have MAPQ=0): --min-MQ 1
				bwa mem -t {threads} -a {input.fa} {input.reads[0]} -R "@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tPL:Illumina" | samtools view -F 2308 --min-MQ 1 -b -@ 2 - > {output} 
				""")


rule samtools_sort:
	input:
		"mapped_reads_per_unit/{sample}-{unit}.bam"
	output:
		temp( "mapped_reads_per_unit/{sample}-{unit}.sorted.bam" )
	shell:
		"""
		samtools sort -T mapped_reads_per_unit/{wildcards.sample}.{wildcards.unit} -O bam {input} > {output}
		"""
		
rule merge_bams_per_sample:
	input:
		bams=get_sample_bams
	output:
		"mapped_reads/{sample}.sorted.bam"
	threads: 4
	shell:
		"""
		samtools merge --threads {threads} {output} {input.bams}
		"""


rule samtools_index:
	input:
		"mapped_reads/{sample}.sorted.bam"
	output:
		"mapped_reads/{sample}.sorted.bam.bai"
	shell:
		"samtools index {input}"


rule samtools_mpileup:
	input:
		ref=genomefile,
		bamfiles=expand("mapped_reads/{sample}.sorted.bam", sample=SAMPLES)
	output:
		"results_raw/pools_mpileup.gz"
	threads: 2
	shell:
		"""
 		echo {input.bamfiles} | sed "s/mapped_reads\///g" | sed "s/.sorted.bam//g" > results_raw/order_of_samples_in_mpileup.txt
		
		# -d 500: use at most 500 reads per input BAM file, apparently these are sampled in order of appearance.
		# --min-MQ 15: minimum mapping quality of an alignment, otherwise skip
		# --no-BAQ : do NOT re-calculate mapping quality (which involves re-aligning). Instead, will use MAPQ as stored in BAM file.
		# --min-BQ INT        skip bases with baseQ/BAQ smaller than INT [13]
		samtools mpileup -d 500 -f {input.ref} --min-MQ 20 --min-BQ 15 --no-BAQ {input.bamfiles} | gzip -c > {output}
		"""


rule run_PoolSNP:
	input:
		mp="results_raw/pools_mpileup.gz",
		ref=genomefile
	output:
		vcf="results_raw/poolsnp.vcf.gz"
	threads: 14
	params:
		MAX_DEPTH_PROPORTION=config["VCF_MAX_DEPTH_PROPORTION"],
		MIN_DEPTH=config["VCF_MIN_DEPTH"],
		MIN_COUNT=config["VCF_MIN_COUNT"]
	shell:
		## miss-frac=0.5 is hardcoded here: we need genotypes even for sites present in only one of the two pools!!
		"""
		thesampleorder=$( cat results_raw/order_of_samples_in_mpileup.txt | tr " " "," )
		
		cd PoolSNP
		bash PoolSNP.sh   \
		mpileup=../{input.mp} \
		reference=../{input.ref} \
		names=$thesampleorder \
		max-cov={params.MAX_DEPTH_PROPORTION} \
		min-cov={params.MIN_DEPTH} \
		min-count={params.MIN_COUNT} \
		min-freq=0.01 \
		miss-frac=0.5 \
		jobs={threads} \
		BS=0 \
		output=poolsnp
		cd ..
				
		mv PoolSNP/poolsnp.vcf.gz {output.vcf}
		tabix {output.vcf} 	
		"""
		



rule get_coverage_in_windows:
	input:
		fa=genomefile,
		bam="mapped_reads/{sample}.sorted.bam",
		bai="mapped_reads/{sample}.sorted.bam.bai"
	output:
		"results_processed/window_cov.{sample}.txt"
	threads: 2
	shell:
		"""
		seqtk comp {input.fa} | awk '{{print $1"\t"$2}}' > genomefile.cov.{wildcards.sample}.txt
		bedtools makewindows -w {windowsize} -g genomefile.cov.{wildcards.sample}.txt > windows.cov.{wildcards.sample}.bed
		rm genomefile.cov.{wildcards.sample}.txt
		
		bedtools multicov -bams {input.bam} -bed windows.cov.{wildcards.sample}.bed > {output} 
		rm windows.cov.{wildcards.sample}.bed
		"""



rule get_cov_ratio:
	input:
		expand("results_processed/window_cov.{sample}.txt", sample=SAMPLES)
	output:
		"results_processed/norm_coverage_ratio_log2.bed.txt"
	shell:
		"""	
		################ get M-F norm read coverage ratio log2 (we add 1e-6 to the read counts to avoid zero division)
		awk '{{ for(i=4; i<=NF;i++) j+=$i; print j; j=0 }}' {input[0]} > total.M
		totm=$( awk '{{sum+=$1;}} END{{print sum;}}' total.M )
		awk -v totalsum="$totm" '{{print ($1*1000000)/totalsum}}' total.M > total.norm.M

		awk '{{ for(i=4; i<=NF;i++) j+=$i; print j; j=0 }}' {input[1]} > total.F
		totf=$( awk '{{sum+=$1;}} END{{print sum;}}' total.F )
		awk -v totalsum="$totf" '{{print ($1*1000000)/totalsum}}' total.F > total.norm.F

		paste total.norm.M total.norm.F | awk '{{print log((($1+1e-6)/($2+1e-6)))/log(2)}}' > log2_ratio

		paste <(cut -f1,2,3 {input[0]} ) log2_ratio > {output}

		rm total.M total.F total.norm.M total.norm.F log2_ratio 	
		"""


rule pi_rawstats:
	input:
		gzvcf="results_raw/poolsnp.vcf.gz"
	output:
		"results_raw/pi_raw.{sample}.txt.gz"		
	shell:
		"""
		# keep only one pool/sample, exclude missing genotypes
		bcftools view -Ov --genotype ^miss -s {wildcards.sample} {input.gzvcf} | python scripts/make_denominator_and_numerator_for_pi.poolseq.py | gzip -c > {output}
		"""


rule pi_windowed:
	input:
		pir="results_raw/pi_raw.{sample}.txt.gz",
		fa=genomefile
	output:
		pi_bed="results_processed/pi_per_window.{sample}.bed.txt"
	threads: 2
	shell:
		"""
		# unzip raw
		gunzip -c {input.pir} > tmp.pi_raw.{wildcards.sample}.txt
		
		seqtk comp {input.fa} | awk '{{print $1"\\t"$2}}' > genomefile.pi.{wildcards.sample}.txt
		bedtools makewindows -w {windowsize} -g genomefile.pi.{wildcards.sample}.txt > windows.pi.{wildcards.sample}.bed

		# ugly construct to handle the "bug" that chromosomes in the VCF are not sorted in 
		# same order as in the fasta genome file; its unclear how this can happen but it DOES
		cut -f1 tmp.pi_raw.{wildcards.sample}.txt | uniq | awk '{{print $1"\\t0\\t9999999999999999"}}' > pi_sort_order_in_scores.{wildcards.sample}.txt
		# find chroms NOT YET in the list
		comm -13 <(cut -f1 pi_sort_order_in_scores.{wildcards.sample}.txt | sort | uniq ) <(cut -f1 genomefile.pi.{wildcards.sample}.txt | sort | uniq) | awk '{{print $1"\\t0\\t9999999999999999"}}' >> pi_sort_order_in_scores.{wildcards.sample}.txt
			 
		bedtools sort -g pi_sort_order_in_scores.{wildcards.sample}.txt -i windows.pi.{wildcards.sample}.bed > windows.pi.sort_order_in_scores.{wildcards.sample}.bed
		
		# make averages of num and denom per window
		(bedtools map -a windows.pi.sort_order_in_scores.{wildcards.sample}.bed -b tmp.pi_raw.{wildcards.sample}.txt -c 4 -o mean > pi_num_{wildcards.sample}.txt )& 
		(bedtools map -a windows.pi.sort_order_in_scores.{wildcards.sample}.bed -b tmp.pi_raw.{wildcards.sample}.txt -c 5 -o mean > pi_denom_{wildcards.sample}.txt )&
		wait
		
		# calc actual pi per window
		paste pi_num_{wildcards.sample}.txt pi_denom_{wildcards.sample}.txt > pi_both_{wildcards.sample}.txt	
		awk '{{ if($8>0) print $1"\\t"$2"\\t"$3"\\t"$4/$8 ; else print $1"\\t"$2"\\t"$3"\\tNA" }}' pi_both_{wildcards.sample}.txt > pi_tmp.{wildcards.sample}.txt 			
		bedtools sort -g genomefile.pi.{wildcards.sample}.txt -i pi_tmp.{wildcards.sample}.txt > {output.pi_bed}
		
		# cleanup temp files
		rm tmp.pi_raw.{wildcards.sample}.txt genomefile.pi.{wildcards.sample}.txt windows.pi.{wildcards.sample}.bed pi_sort_order_in_scores.{wildcards.sample}.txt windows.pi.sort_order_in_scores.{wildcards.sample}.bed pi_num_{wildcards.sample}.txt pi_denom_{wildcards.sample}.txt pi_both_{wildcards.sample}.txt pi_tmp.{wildcards.sample}.txt

		"""


rule dxy_rawstats:
	input:
		gzvcf="results_raw/poolsnp.vcf.gz"
	output:
		"results_raw/dxy_raw.txt.gz"		
	shell:
		"""
		bcftools view -Ov -M 2 --genotype ^miss {input.gzvcf} | python scripts/make_denominator_and_numerator_for_dxy.two_pools.py  | gzip -c > {output}
		"""


rule dxy_windowed:
	input:
		raw="results_raw/dxy_raw.txt.gz",
		fa=genomefile
	output:
		"results_processed/dxy_per_window.bed.txt"
	threads: 2
	shell:
		"""
		# unzip raw file
		gunzip -c {input.raw} > tmp.dxy_raw.txt
		
		seqtk comp {input.fa} | awk '{{print $1"\\t"$2}}' > genomefile.dxy.txt
		bedtools makewindows -w {windowsize} -g genomefile.dxy.txt > windows.dxy.bed
		
		# ugly construct to handle the "bug" that chromosomes in the VCF are not sorted in 
		# same order as in the fasta genome file; its unclear how this can happen but it DOES
		cut -f1 tmp.dxy_raw.txt | uniq | awk '{{print $1"\\t0\\t9999999999999999"}}' > dxy_sort_order_in_scores.txt
		# find chroms NOT YET in the list
		comm -13 <(cut -f1 dxy_sort_order_in_scores.txt | sort | uniq ) <(cut -f1 genomefile.dxy.txt | sort | uniq) | awk '{{print $1"\\t0\\t9999999999999999"}}' >> dxy_sort_order_in_scores.txt

		bedtools sort -g dxy_sort_order_in_scores.txt -i windows.dxy.bed > windows.dxy.sort_order_in_scores.bed

		# correct
		( bedtools map -a windows.dxy.sort_order_in_scores.bed -b tmp.dxy_raw.txt -c 4 -o mean > dxy_num )&
		( bedtools map -a windows.dxy.sort_order_in_scores.bed -b tmp.dxy_raw.txt -c 5 -o mean > dxy_denom )&
		wait

		paste dxy_num dxy_denom > dxy_both
		
		awk '{{ if($8>0) print $1"\\t"$2"\\t"$3"\\t"$4/$8 ; else print $1"\\t"$2"\\t"$3"\\tNA" }}' dxy_both > dxy_tmp
			
		bedtools sort -g genomefile.dxy.txt -i dxy_tmp > {output}
		
		rm genomefile.dxy.txt windows.dxy.bed dxy_num dxy_denom dxy_both dxy_tmp windows.dxy.sort_order_in_scores.bed dxy_sort_order_in_scores.txt tmp.dxy_raw.txt
		"""


rule netdiv_windowed:
	input:
		dxy="results_processed/dxy_per_window.bed.txt",
		pi_bed="results_processed/pi_per_window.{sample}.bed.txt",
		fa=genomefile
	output:
		"results_processed/netDiv_per_window.{sample}.bed.txt"
	shell:
		"""
		seqtk comp {input.fa} | awk '{{print $1"\\t"$2}}' > genomefile.netdiv.{wildcards.sample}.txt
		bedtools makewindows -w {windowsize} -g genomefile.netdiv.{wildcards.sample}.txt > windows.netdiv.{wildcards.sample}.bed
				
		paste {input.dxy} {input.pi_bed} > netdiv_prep.{wildcards.sample}.txt
		awk '{{ if($8!="NA" && $4!="NA") print $1"\\t"$2"\\t"$3"\\t"$4-$8 ; else print $1"\\t"$2"\\t"$3"\\tNA" }}' netdiv_prep.{wildcards.sample}.txt > {output}

		rm genomefile.netdiv.{wildcards.sample}.txt windows.netdiv.{wildcards.sample}.bed netdiv_prep.{wildcards.sample}.txt
		"""


rule fst_rawstats:
	input:
		"results_raw/poolsnp.vcf.gz"
	output:
		"results_raw/fst_raw.txt.gz"		
	shell:
		"""
		bcftools view -g ^miss -m2 -M2 -Ov {input} | python scripts/make_denominator_and_numerator_for_Hudson_Fst.two_pools.py | gzip -c > {output} 
		"""



rule fst_windowed:
	input:
		fraw="results_raw/fst_raw.txt.gz",
		fa=genomefile
	output:
		"results_processed/Fst_Hudson_per_window.bed.txt"
	shell:
		"""
		## as found to be more robust to presence of many rare variants, we get a "ratio of averages" rather than an "average of ratios"..
		## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3759727/
		
		# unzip raw
		gunzip -c {input.fraw} > tmp.fst_raw.txt
		
		seqtk comp {input.fa} | awk '{{print $1"\\t"$2}}' > genomefile.fst.txt
		bedtools makewindows -w {windowsize} -g genomefile.fst.txt > windows.fst.bed
		
		fstlines=$(cat tmp.fst_raw.txt | wc -l)
		if [ "$fstlines" -gt 2 ]; then
			bedtools map -a windows.fst.bed -b tmp.fst_raw.txt -g genomefile.fst.txt -c 4 -o mean > fst_num
			bedtools map -a windows.fst.bed -b tmp.fst_raw.txt -g genomefile.fst.txt -c 5 -o mean > fst_denom

			paste fst_num fst_denom > fst_both
		
			awk '{{ if($8>0) print $1"\\t"$2"\\t"$3"\\t"1-($4/$8) ; else print $1"\\t"$2"\\t"$3"\\tNA" }}' fst_both > {output}
			
			
			
		else
			# fst file was empty, because all values in genome were NA. Make an NA windowed-file to allow downstream to proceed anyway.
			cat windows.fst.bed | awk '{{print $0"\\tNA"}}' > {output}
		fi
		
		rm genomefile.fst.txt windows.fst.bed tmp.fst_raw.txt fst_num fst_denom fst_both 		
		"""

	
rule pseudo_phase_gametologs:
	input:
		"results_raw/poolsnp.vcf.gz"
	output:
		"results_raw/gametolog_candidate_alleles.allsites.vcf.gz"
	shell:
		"""
		bcftools view -M2 -Ov {input} | python scripts/pseudo_phase_gametologs.poolseq.py | bgzip -c > {output}
		"""


	
rule score_XY_gametolog_divergence:
	input:
		"results_raw/gametolog_candidate_alleles.allsites.vcf.gz"
	output:
		"results_raw/gametolog_candidate_alleles.XY_divergence.rawstats.txt.gz"
	shell:
		"""
		echo -e "XY_Y_like\\tXY_Y_like" > xypopm
		echo -e "XY_X_like\\tXY_X_like" >> xypopm
		
		echo -e "XY_X_like" > xypop
		echo -e "XY_Y_like" >> xypop

		vcftools --gzvcf {input} --keep xypop --recode --stdout | python scripts/make_denominator_and_numerator_for_dxy.py --popmap xypopm | gzip -c > {output}
		rm xypopm xypop
		"""

rule score_ZW_gametolog_divergence:
	input:
		"results_raw/gametolog_candidate_alleles.allsites.vcf.gz"
	output:
		"results_raw/gametolog_candidate_alleles.ZW_divergence.rawstats.txt.gz"
	shell:
		"""
		echo -e "ZW_W_like\\tZW_W_like" > zwpopm
		echo -e "ZW_Z_like\\tZW_Z_like" >> zwpopm
		
		echo -e "ZW_W_like" > zwpop
		echo -e "ZW_Z_like" >> zwpop

		vcftools --gzvcf {input} --keep zwpop --recode --stdout | python scripts/make_denominator_and_numerator_for_dxy.py --popmap zwpopm | gzip -c > {output}
		rm zwpopm zwpop
		"""

rule XY_div_windows:
	input:
		raw="results_raw/gametolog_candidate_alleles.XY_divergence.rawstats.txt.gz",
		fa=genomefile
	output:
		"results_processed/gametolog_candidate_alleles.XY_divergence.windows.bed.txt"
	threads: 2
	shell:
		"""
		# unzip raw
		gunzip -c {input.raw} > tmp.gametolog_candidate_alleles.XY_divergence.rawstats.txt
		
		seqtk comp {input.fa} | awk '{{print $1"\\t"$2}}' > genomefile.xygametologs.txt
		bedtools makewindows -w {windowsize} -g genomefile.xygametologs.txt > windows.xygametologs.bed

		(bedtools map -a windows.xygametologs.bed -b tmp.gametolog_candidate_alleles.XY_divergence.rawstats.txt -g genomefile.xygametologs.txt -c 4 -o mean > XY_dxy_num )&
		(bedtools map -a windows.xygametologs.bed -b tmp.gametolog_candidate_alleles.XY_divergence.rawstats.txt -g genomefile.xygametologs.txt -c 5 -o mean > XY_dxy_denom )&
		wait
		
		paste XY_dxy_num XY_dxy_denom > XY_dxy_both
		
		awk '{{ if($8>0) print $1"\\t"$2"\\t"$3"\\t"$4/$8 ; else print $1"\\t"$2"\\t"$3"\\tNA" }}' XY_dxy_both > {output}
		
		rm genomefile.xygametologs.txt windows.xygametologs.bed XY_dxy_denom XY_dxy_num XY_dxy_both	tmp.gametolog_candidate_alleles.XY_divergence.rawstats.txt	
		"""

rule ZW_div_windows:
	input:
		raw="results_raw/gametolog_candidate_alleles.ZW_divergence.rawstats.txt.gz",
		fa=genomefile
	output:
		"results_processed/gametolog_candidate_alleles.ZW_divergence.windows.bed.txt"
	threads: 2
	shell:
		"""
		# unzip raw
		gunzip -c {input.raw} > tmp.gametolog_candidate_alleles.ZW_divergence.rawstats.txt

		seqtk comp {input.fa} | awk '{{print $1"\\t"$2}}' > genomefile.zwgametologs.txt
		bedtools makewindows -w {windowsize} -g genomefile.zwgametologs.txt > windows.zwgametologs.bed

		( bedtools map -a windows.zwgametologs.bed -b tmp.gametolog_candidate_alleles.ZW_divergence.rawstats.txt -g genomefile.zwgametologs.txt -c 4 -o mean > ZW_dxy_num )&
		( bedtools map -a windows.zwgametologs.bed -b tmp.gametolog_candidate_alleles.ZW_divergence.rawstats.txt -g genomefile.zwgametologs.txt -c 5 -o mean > ZW_dxy_denom )&
		wait
				
		paste ZW_dxy_num ZW_dxy_denom > ZW_dxy_both
		
		awk '{{ if($8>0) print $1"\\t"$2"\\t"$3"\\t"$4/$8 ; else print $1"\\t"$2"\\t"$3"\\tNA" }}' ZW_dxy_both > {output}
		
		rm genomefile.zwgametologs.txt windows.zwgametologs.bed ZW_dxy_denom ZW_dxy_num ZW_dxy_both	tmp.gametolog_candidate_alleles.ZW_divergence.rawstats.txt	
		"""


rule KMC_count:
	input:
		get_fastq_sample_ALL
	output:
		pre="results_raw/kmer_{sample}.kmc_pre",
		suf="results_raw/kmer_{sample}.kmc_suf"
	threads: 14
	shell:
		"""
		# count kmers that occur at least 2 times (-ci), and count them up to 1 million (-cs)
		echo {input} | tr " " "\n" >input_file_names.{wildcards.sample}
		mkdir -p ./kmc_tempdir.{wildcards.sample}
		kmc -m5 -sm -k16 -fq -ci2 -cs1000000 -t{threads} @input_file_names.{wildcards.sample} results_raw/kmer_{wildcards.sample} ./kmc_tempdir.{wildcards.sample}
		rm input_file_names.{wildcards.sample}
		rm -r ./kmc_tempdir.{wildcards.sample}
		"""


rule KMC_substract:
	input:
		expand("results_raw/kmer_{sample}.kmc_suf", sample=SAMPLES)
	output:
		diff1="results_raw/kmer_diff_ci{cov}_" + SAMPLES[0] + "_minus_" + SAMPLES[1] + ".kmc_suf",
		diff2="results_raw/kmer_diff_ci{cov}_" + SAMPLES[1] + "_minus_" + SAMPLES[0] + ".kmc_suf"
	threads: 14
	shell:
		"""
		kmc_tools -t{threads} simple results_raw/kmer_{SAMPLES[0]} -ci{wildcards.cov} results_raw/kmer_{SAMPLES[1]} kmers_subtract results_raw/kmer_diff_ci{wildcards.cov}_{SAMPLES[0]}_minus_{SAMPLES[1]}
		kmc_tools -t{threads} simple results_raw/kmer_{SAMPLES[1]} -ci{wildcards.cov} results_raw/kmer_{SAMPLES[0]} kmers_subtract results_raw/kmer_diff_ci{wildcards.cov}_{SAMPLES[1]}_minus_{SAMPLES[0]}
		"""

rule evaluate_KMC_ci:
	input:
		expand("results_raw/kmer_diff_ci{cov}_" + SAMPLES[0] + "_minus_" + SAMPLES[1] + ".kmc_suf", cov=kmer_covs),
		expand("results_raw/kmer_diff_ci{cov}_" + SAMPLES[1] + "_minus_" + SAMPLES[0] + ".kmc_suf", cov=kmer_covs)	
	output:
		"results_raw/kmers_diff_ci_evaluation.txt"
	run:
		shell("echo ci {SAMPLES[0]}_minus_{SAMPLES[1]} {SAMPLES[1]}_minus_{SAMPLES[0]} > {output}")
		for civalue in kmer_covs:
			print(civalue)
			shell("echo {civalue} > kk")
			shell("kmc_tools transform results_raw/kmer_diff_ci{civalue}_{SAMPLES[0]}_minus_{SAMPLES[1]} dump /dev/stdout | wc -l > evtemp1")
			shell("kmc_tools transform results_raw/kmer_diff_ci{civalue}_{SAMPLES[1]}_minus_{SAMPLES[0]} dump /dev/stdout | wc -l > evtemp2")
			shell("paste kk evtemp1 evtemp2 >> {output}")
			shell("rm kk evtemp1 evtemp2")   
		
		# find the ci under which one of the pools/samples has greatest number of private kmers, "best ci"
		# this is the maximum of the absolutes of the log2 ratios of the numbers of kmers private to pool1 over those private to pool2
		# it should minimise false-privates and maximise true-privates
		import math
		abs_log2ratios = []
		cis = []
		with open("results_raw/kmers_diff_ci_evaluation.txt", "r") as I:
			I.readline()
			for line in I:
				fields = line.strip().split()
				cis.append(fields[0])
				abs_log2ratio = abs( math.log( (float(fields[1])+1e-6) / (float(fields[2])+1e-6), 2 ) )
				abs_log2ratios.append( abs_log2ratio )		
		best_ci = cis[ abs_log2ratios.index(max(abs_log2ratios)) ]
		print("best ci = ", best_ci)
		
		# now we get these "best ci" kmers in fasta format
		shell("""
		kmc_tools transform results_raw/kmer_diff_ci{best_ci}_{SAMPLES[0]}_minus_{SAMPLES[1]} dump /dev/stdout | awk '{{print ">mer\\n"$1}}' > results_raw/kmers_best_ci.{SAMPLES[0]}_minus_{SAMPLES[1]}.fasta
 		kmc_tools transform results_raw/kmer_diff_ci{best_ci}_{SAMPLES[1]}_minus_{SAMPLES[0]} dump /dev/stdout | awk '{{print ">mer\\n"$1}}' > results_raw/kmers_best_ci.{SAMPLES[1]}_minus_{SAMPLES[0]}.fasta
 		""")
	

	
rule match_kmers_to_genome:
	input:
		fa=genomefile,
		gidx=genomefile +".bwt",
		ci_eval="results_raw/kmers_diff_ci_evaluation.txt"
	output:
		expand("results_processed/kmers.{sample}_specific.per_window.bed.txt", sample=SAMPLES)	
	shell:
		"""
		# searching for perfect matches of 21-mers: https://bioinformatics.stackexchange.com/questions/7298/mapping-heteryzygous-kmers-on-a-genome
		# The -k 21 says to use a minimum seed length of 21, which forces exact matches. The -T 21 requires a minimum score of 21, which also enforces an exact match. The -a parameter reports all matches, since a "best" match doesn't make sense in this situation. The -c parameter limits how many matches are reported, which may need to be adjusted depending on how repetitive the 21-mer is.
		# we match IMPERFECTLY = allowing 1 mismatch in 16 bp kmers. => 6.25% divergence
		
		# use samtools view -F 4 : removes unmapped
	
		bwa mem -t 1 -k 15 -T 15 -a -c 5000 {input.fa} results_raw/kmers_best_ci.{SAMPLES[1]}_minus_{SAMPLES[0]}.fasta | samtools view -F 4 -b - > specific_to_{SAMPLES[1]}_kmer.bam
		samtools sort -T specific_to_{SAMPLES[1]}_kmer -O bam specific_to_{SAMPLES[1]}_kmer.bam > specific_to_{SAMPLES[1]}_kmer.sorted.bam
		samtools index specific_to_{SAMPLES[1]}_kmer.sorted.bam
		rm specific_to_{SAMPLES[1]}_kmer.bam

		bwa mem -t 1 -k 15 -T 15 -a -c 5000 {input.fa} results_raw/kmers_best_ci.{SAMPLES[0]}_minus_{SAMPLES[1]}.fasta | samtools view -F 4 -b - > specific_to_{SAMPLES[0]}_kmer.bam
		samtools sort -T specific_to_{SAMPLES[0]}_kmer -O bam specific_to_{SAMPLES[0]}_kmer.bam > specific_to_{SAMPLES[0]}_kmer.sorted.bam
		samtools index specific_to_{SAMPLES[0]}_kmer.sorted.bam
		rm specific_to_{SAMPLES[0]}_kmer.bam

		seqtk comp {input.fa} | awk '{{print $1"\\t"$2}}' > genomefile.kmer.txt
		bedtools makewindows -w {windowsize} -g genomefile.kmer.txt > windows.kmer.bed
		rm genomefile.kmer.txt
		bedtools multicov -bams specific_to_{SAMPLES[0]}_kmer.sorted.bam -bed windows.kmer.bed > results_processed/kmers.{SAMPLES[0]}_specific.per_window.bed.txt
		bedtools multicov -bams specific_to_{SAMPLES[1]}_kmer.sorted.bam -bed windows.kmer.bed > results_processed/kmers.{SAMPLES[1]}_specific.per_window.bed.txt
		rm windows.kmer.bed specific_to_{SAMPLES[0]}_kmer.sorted.bam specific_to_{SAMPLES[0]}_kmer.sorted.bam.bai specific_to_{SAMPLES[1]}_kmer.sorted.bam specific_to_{SAMPLES[1]}_kmer.sorted.bam.bai
		"""



rule plot_all:
	input:
		regions=regions_for_plot_bed,
		a="results_processed/norm_coverage_ratio_log2.bed.txt",
		b="results_processed/kmers."+SAMPLES[0]+"_specific.per_window.bed.txt", 
		c="results_processed/kmers."+SAMPLES[1]+"_specific.per_window.bed.txt",
		d="results_processed/pi_per_window."+SAMPLES[0]+".bed.txt",
		e="results_processed/pi_per_window."+SAMPLES[1]+".bed.txt",
		f="results_processed/dxy_per_window.bed.txt",
		g="results_processed/netDiv_per_window."+SAMPLES[0]+".bed.txt",
		h="results_processed/netDiv_per_window."+SAMPLES[1]+".bed.txt",
		i="results_processed/Fst_Hudson_per_window.bed.txt",
		j="results_processed/gametolog_candidate_alleles.XY_divergence.windows.bed.txt",
		k="results_processed/gametolog_candidate_alleles.ZW_divergence.windows.bed.txt"
	output:
		stats="results_processed/allstats.txt",
		regionstats="results_processed/region.stats.txt",
		plot="results_processed/allstats.plots.pdf",
		regionplot="results_processed/region.stats.plots.pdf"
	shell:
		"""
		cat {input.a} > {output.stats}

		for i in {input.b} {input.c} {input.d} {input.e} {input.f} {input.g} {input.h} {input.i} {input.j} {input.k} ; do
			paste {output.stats} <(cut -f4 $i ) > tmpst
			mv tmpst {output.stats}
		done
		
		bedtools intersect -wa -a {output.stats} -b {input.regions} > {output.regionstats}

		Rscript scripts/plot_allstats.R {output.stats} {output.plot}
		Rscript scripts/plot_allstats.R {output.regionstats} {output.regionplot}	
		
		"""


	
