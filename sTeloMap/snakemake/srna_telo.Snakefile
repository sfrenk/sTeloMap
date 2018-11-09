import glob
import re
import sys
import os

### NEW AND IMPROVED PIPELINE ###
#	Genome mapping with butter


############################    PARAMETERS    ##############################

# Input parameters
BASEDIR = "raw_compiled"
EXTENSION = ".txt.gz"

# Filtering parameters
FILTER_BASE = "A,T,C,G"
SIZE = "18,30"
TRIM = False
MIN_TRIM_LENGTH = 0

# Mapping parameters
TELO_INDEX = "/nas/longleaf/home/sfrenk/proj/seq/telomere/bowtie/telomere"

# Other parameters
UTILS_DIR = "/nas/longleaf/home/sfrenk/scripts/util"

###############################################################################

# Get bowtie index for species
refs = {"elegans" : "/nas/longleaf/home/sfrenk/proj/seq/WS251/genome/bowtie/genome.fa", "remanei" : "/nas/longleaf/home/sfrenk/proj/seq/remanei/bowtie/genome.fa", "briggsae" : "/nas/longleaf/home/sfrenk/proj/seq/briggsae/WS263/bowtie/genome.fa"}

def get_ref(sample_name):

	species=re.search("elegans|briggsae|remanei", sample_name)
	if species is None:
		species="elegans"
	else:
		species = species.group(0)
	
	species_ref = refs[species]
	
	return(species_ref)

#if SPECIES not in refs:
#	print("ERROR: Unknown reference!")
#else:
#	REF = refs[SPECIES]

DATASET = re.search("([^/]*)/?$", BASEDIR).group(1)
SAMPLES = glob.glob(BASEDIR + "/*" + EXTENSION)
SAMPLES = [ re.search(BASEDIR + "/?([^/]+)" + EXTENSION, x).group(1) for x in SAMPLES ]


if len(SAMPLES) == 0:
	sys.exit("ERROR: no samples in base directory!")

rule all:
    input:
        expand("genome/{sample}_genome.bam", sample = SAMPLES),
        expand("results/{sample}.db", sample = SAMPLES),
        "results/" + DATASET + ".db"

rule filter_srna:
	input:
		BASEDIR + "/{sample}" + EXTENSION
	output:
		"filtered/{sample}.fa"
	params:
		filter_base = FILTER_BASE,
		size = SIZE,
		trim = TRIM,
		min_trim_length = MIN_TRIM_LENGTH,
		utils_dir = UTILS_DIR
	log:
		"logs/{sample}_filter.log"
	shell:
		"module add python; "
		"python3 {params.utils_dir}/small_rna_filter.py \
		-f {params.filter_base} \
		-s {params.size} \
		-t {params.trim} \
		-m {params.min_trim_length} \
		-o {output} \
		{input} > {log} 2>&1"

rule butter_mapping_genome:
	input:
		"filtered/{sample}.fa"
	output:
		bamfile = "genome/{sample}_genome.bam",
		bamidx = "genome/{sample}_genome.bam.bai"
	params:
		ref = lambda wildcards: get_ref('{sample}'.format(sample = wildcards.sample)),
		sample_name="{sample}"
	log:
		"logs/{sample}_map_genome.log"
	threads: 8
	shell:
		"module purge; " 
		"module add bowtie/1.1.2 perl samtools/0.1.19 && \
		butter --no_condense --aln_cores {threads} --max_rep 10000 --bam2wig none {input} {params.ref} > {log} 2>&1; "
		"mv filtered/{params.sample_name}.bam genome/{params.sample_name}_genome.bam; "
		"mv filtered/{params.sample_name}.bam.bai genome/{params.sample_name}_genome.bam.bai; "

rule map_telo:
	input:
		"filtered/{sample}.fa"
	output:
		"telo/{sample}.bam"
	params:
		telo_index = TELO_INDEX
	threads: 8
	log:
		"logs/{sample}_map_telo.log"
	shell:
		"module purge; " 
		"module add bowtie/1.1.2 samtools/1.8 && \
		bowtie -f --best --strata -M 1 -S -v 3 -p {threads} {params.telo_index} {input} | samtools view -bh -F 4 - | samtools sort -o {output} - > {log} 2>&1"

rule index_telo_bam:
	input:
		"telo/{sample}.bam"
	output:
		"telo/{sample}.bam.bai"
	log:
		"logs/{sample}_index_telo_bam.log"
	shell:
		"module add samtools/1.8; "
		"samtools index {input}"

rule extract_telo_reads:
	input:
		telo_bam = "telo/{sample}.bam",
		telo_idx = "telo/{sample}.bam.bai",
		genome_bam = "genome/{sample}_genome.bam",
		genome_idx = "genome/{sample}_genome.bam.bai"
	output:
		reads_file = "results/{sample}_reads.txt",
		alignments_file = "results/{sample}_alignments.txt",
		alignments_bam = "results/{sample}_telo.bam",
		alignments_bam_idx = "results/{sample}_telo.bam.bai"
	params:
		utils_dir = UTILS_DIR,
		dataset = DATASET,
		sample_name = "{sample}"
	log:
		"logs/{sample}_extract_telo_reads.log"
	shell:
		"python3 {params.utils_dir}/extract_telo_reads.py -d {params.dataset} -s {params.sample_name} -o results/{params.sample_name} -t {input.telo_bam} -g {input.genome_bam} &> {log}"

rule get_sample_info:
	input:
		telo_reads_file = "results/{sample}_reads.txt",
		genome_bam = "genome/{sample}_genome.bam",
		bam_index = "genome/{sample}_genome.bam.bai"
	output:
		"sample_info/{sample}_info.txt"
	params:
		dataset = DATASET,
		sample_name = "{sample}"
	log:
		"logs/{sample}_get_sample_info.log"
	run:
		"module add samtools; "
		shell('''telo_reads=$(wc -l < {input.telo_reads_file}); module add samtools/1.8 && samtools view -c {input.genome_bam} | awk -v var="$telo_reads" '{{ print "{params.dataset}_{params.sample_name}\t{params.dataset}\t{params.sample_name}\t"$0"\t"var }}' > {output} 2> {log}''')

rule store_results:
	input:
		alignment_file = "results/{sample}_alignments.txt",
		reads_file = "results/{sample}_reads.txt",
		sample_info = "sample_info/{sample}_info.txt"
	output:
		"results/{sample}.db"
	params:
		utils_dir = UTILS_DIR
	log:
		"logs/{sample}_store_results.log"
	shell:
		"module add python; "
		"python3 {params.utils_dir}/ngs_store.py -d {output} -s {input.sample_info} -r {input.reads_file} -a {input.alignment_file} &> {log}"

rule merge_db:
	input:
		expand("results/{sample}.db", sample = SAMPLES)
	output:
		"results/" + DATASET + ".db"
	params:
		utils_dir = UTILS_DIR
	log:
		"logs/merge_db.log"
	shell:
		"module add python; "
		"python3 {params.utils_dir}/merge_db.py -o {output} {input} &> {log}"
