# delim_supergene_poolseq

a pipeline to delimit supergene-like regions in chromosome-scale genome assemblies, using POOLED re-sequencing data of two groups 

starting from a genome .fasta and fastq reads, generates statistics in windows along the genome:

- group difference in read coverage

- group-specific kmers (16-mers), aligned with 1 mismatch

- nucleotide diversity per group

- Fst between groups

- abs divergence (dxy) between groups

- net divergence between groups

- XY-like variants per site

- ZW-like variants per site


## installation

### clone this repo
```
git clone https://github.com/mscharmann/delim_supergene_poolseq
```
### setup a conda environment

## create stepwise
```
conda create --name delimit_sexregions
conda activate delimit_sexregions
conda install snakemake=5.4 scipy bwa samtools bedtools seqtk vcftools bcftools tabix plink parallel freebayes -y
conda install -c conda-forge -c bioconda -c defaults vcflib -y
conda install -c conda-forge r-ggplot2 r-cowplot -y
```
## OR use this YAML:

### modify prefix of installation path in last line of this file, then
```
conda env create --file delimit_sexregions.2021-06-30.yml
```

### get poolSNP
This is a heuristic SNP caller for poolseq data:

```
cd delim_supergene_poolseq

git clone https://github.com/capoony/PoolSNP
```

### modify poolSNP
poolSNP will report only the SNPs, however, we also want the invariant, fixed sites reported.
Modify the script PoolSNP/scripts/PoolSnp.py by replacing lines 202 - 216 with the following code:

```
    # test if site is polymorphic
    is_polym = True
    if len(totalalleles) < 2:
 #       print (CHR,POS,"non-poly",totalalleles)
        is_polym = False

    # create output for VCF
    if is_polym:
        ALT = []
        # set alternative allele order:
        for i in ["A", "T", "C", "G"]:
            if i == REF:
                continue
            if i not in totalalleles:
                continue
            ALT.append(i)
    else:
        ALT = ["."]
```
Furthermore, we need BGZIP compressed output. Replace line 284 of the script PoolSNP/PoolSNP.sh with the following:
```
cat $out/temp/header.txt $out/temp/SNPs.txt | bgzip > $out.vcf.gz
```
We also do not need the "bad_sites" output from poolSNP. To save time, delete lines 288 ff from the script PoolSNP/PoolSNP.sh


## run

### config.yaml: provide meta-data and set parameters
This file should be mostly self-explanatory and also contains some comments. Further important points are:
- samples_units_fqs_map = a table listing the pools, samples and "mapping" to their corresponding fastq files. Each sample or pool may have multiple (pairs) of fastq files, called a 'unit'. Both absolute paths and paths relative to the working directory should be tolerated. Both single and paired-end data can be handled, also mixed.
- genome fasta file path
- regions_for_plot_bed: a BED file with a region of special interest, e.g. the sex chromosome. The pipeline produces statistics for all genomic windows and then extracts a subset (region). Plots will be made for all genomic windows, and for the subset region. This can be useful if the assembly contains many contigs that would make a messy plot.  
- windowsize = length of intervals (windows) in which to split all statistics and plots. Most useful values are between 1000 and 100000.
- VCF_MAX_DEPTH_PROPORTION = a poolSNP parameter. Will discard sites that exceed the Xth percentile of coverage. Must use to exclude highly repetitive regions.
- VCF_MIN_DEPTH = a poolSNP parameter. Minimuzm read depth across all pools to consider a site ; e.g. min-cov=10 will only consider positions with a minimum coverage >10
- VCF_MIN_COUNT = a poolSNP parameter. The minimum count of an allele over both pools to consider this allele.
- avg_samples_per_pool = a parameter required for assessment of group-specific kmers by KMC (-ci parameter). Searches for group-specific kmers will be run with multiples of this value, and the best will be selected as the one with the greatest difference in group-specific kmers between groups. This logic applies when only one of the two groups is expected to contain specific kmers.


### local:
```
snakemake -j24 --keep-going --rerun-incomplete
```
### on SLURM cluster:
```
snakemake -j 500 --cluster-config cluster.axiom.json --cluster "sbatch -p {cluster.partition} -t {cluster.time} -c {cluster.CPUs} --mem={cluster.RAM_memory}" --restart-times 3 --keep-going --rerun-incomplete
```
### on LSF cluster:
```
snakemake -j 500 --cluster-config cluster.EULER.json --cluster "bsub -W {cluster.time} -n {cluster.CPUs} -R {cluster.mem_and_span}" --restart-times 3 --keep-going --rerun-incomplete
```
## post-run
- Inspect the PDF files and .txt tables with statistics in the directory results_processed.
- results_processed are specific for the given windowsize. If a different windowsize is desired, intermediate results in the directory results_raw can be used to generate these relatively quickly, without going back to the read data and variant-calling. Just move the directory results_raw, change the windowsize in config.yaml, and start the pipeline again exactly as before. I usually do this several times for window sizes 1 kb, 10 kb, 25 kb, 50 kb, 100 kb.
- If a different subset region for plotting is desired, there is no need to run the pipeline again, just subset the file calls/allstats.txt and call scripts/plot_allstats.R . For example, to plot only assembly sequences larger than 5 Gbp  
```
seqtk comp data/genome_assembly.fa | awk '{{if($2>5000000) print $1"\t0\t"$2}}' > bigchroms.bed

bedtools intersect -wa -a calls/allstats.txt -b bigchroms.bed > subs.txt

Rscript scripts/plot_allstats.R subs.txt bigchroms.pdf
```
- cleaning up / archiving: Once the final stats for all desired windowsizes are produced, it makes sense to clean up large intermediate files. I suggest to keep only the config.yaml, the popmap and samples_reads_map, the results_processed, or if you have a bit of space, also results_raw. I would delete the large .BAM files and intermediate VCF files; in the worst case these can be re-calculated. To clean up in this way, run 
```
rm -rf mapped_reads FB_chunks FB_chunk_VCFs FB_chunk_VCFs_filtered FB_regionsALL.bed normalization_coefficients.txt
```
- I would also get rid of the hidden .snakemake directory, unless you really need it. It can be quite large/numerous tiny files.




# Simulate test data: an XY system with the sexchroms and one autosome

```
import numpy as np
import gzip

n_autosomes = 1

autosomes = []
for i in range(n_autosomes):
	autosomes.append( "".join( np.random.choice(["A","C","T","G"], size = 600000, replace = True) ) )

```
the sexchroms are diverged by a 2 MB region between 2Mbp to 4 Mbp

- Y-hemizygous region 1 Mbp
- X-hemizygous region 1 Mbp => Y chrom is SHORTER
- X-Y gametolog region: divergence = 10%
- PARs on both sides.
```
Xchrom = "".join( np.random.choice(["A","C","T","G"], size = 600000, replace = True) )

```
construct the Ychrom: PAR-X_Y_gametolog_region-Y_hemizygous

1	100000	PAR1

100001	200000	X_Y_gametolog_region

200001	300000	Y-hemizygous

300001	500000	PAR2


and therefore the X chrom is: with the PAR2 being 4-6 Mbp on the X, but 3-5 Mb on the Y (because the X-hemiz is larger (2x) than the Y-hemiz region)

1	100000	PAR1

100001	200000	X_Y_gametolog_region

200001	400000	X-hemizygous

400001	600000	PAR2

```
X_Y_gametolog_region = list( Xchrom[100000:200000] )
snpsites = np.random.uniform(0, len(X_Y_gametolog_region), size = int(0.1*len(X_Y_gametolog_region)) )
print len(snpsites)
for s in snpsites:
	s = int(s)
	isnuc = X_Y_gametolog_region[s]
	newnuc = np.random.choice([ x for x in ["A","C","T","G"] if not x == isnuc])
	X_Y_gametolog_region[s] = newnuc 


X_Y_gametolog_region = "".join(X_Y_gametolog_region)
Y_hemiz_region = "".join( np.random.choice(["A","C","T","G"], size = 100000, replace = True) ) 

Ychrom = Xchrom[:100000] + X_Y_gametolog_region + Y_hemiz_region + Xchrom[400000:]


len(Ychrom)
len(Xchrom)


with open("fakegenome.FEMALE.fa", "w") as O:
	cnt = 0
	for chr in autosomes:
		cnt += 1
		O.write(">chrom" + str(cnt) + "\n")
		O.write(chr + "\n")
	O.write(">chrom_X" + "\n")
	O.write(Xchrom + "\n")


with open("fakegenome.FEMALE.diploid.fa", "w") as O:
	cnt = 0
	for chr in autosomes:
		cnt += 1
		O.write(">chrom" + str(cnt) + "_1" + "\n")
		O.write(chr + "\n")
		O.write(">chrom" + str(cnt) + "_2" + "\n")
		O.write(chr + "\n")
	O.write(">chrom_X" + "_1" + "\n")
	O.write(Xchrom + "\n")
	O.write(">chrom_X" + "_2" + "\n")
	O.write(Xchrom + "\n")


with open("fakegenome.MALE.fa", "w") as O:
	cnt = 0
	for chr in autosomes:
		cnt += 1
		O.write(">chrom" + str(cnt) + "\n")
		O.write(chr + "\n")
	O.write(">chrom_Y" + "\n")
	O.write(Ychrom + "\n")


with open("fakegenome.MALE.diploid.fa", "w") as O:
	cnt = 0
	for chr in autosomes:
		cnt += 1
		O.write(">chrom" + str(cnt) + "_1" + "\n")
		O.write(chr + "\n")
		O.write(">chrom" + str(cnt) + "_2" + "\n")
		O.write(chr + "\n")
	O.write(">chrom_Y" + "\n")
	O.write(Ychrom + "\n")
	O.write(">chrom_X" + "\n")
	O.write(Xchrom + "\n")

```

now simulate some WGS read data!

use wgsim from samtools package

```
for i in {1..3}; do
wgsim -N 160000 -1 150 -2 150 fakegenome.MALE.diploid.fa sample_M_${i}.1.fastq sample_M_${i}.2.fastq
done


for i in {1..3}; do
wgsim -N 160000 -1 150 -2 150 fakegenome.FEMALE.diploid.fa sample_F_${i}.1.fastq sample_F_${i}.2.fastq
done

gzip *.fastq
```

