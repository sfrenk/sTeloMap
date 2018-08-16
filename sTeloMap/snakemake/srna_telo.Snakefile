import glob
import re
import sys
import os

############################    DESCRIPTION    ##############################

# Pipeline for analyzing telomere-related small RNAs (telo-sRNAs).
# Reads mapping to telomeres with 0-3 mismatches are identified and stored in a SQLite3 database along with library metadata.

# Usage:
#	1. Use setup_working_dir.sh to set up the working directory for analysis:
		#	setup_working_dir.sh -p srna_telo -d <directory containing reads>
#	2. Run the snakemake script
		# bash run_snakemake.sh

############################    PARAMETERS    ##############################

# Input parameters
BASEDIR = ""
EXTENSION = ".fa.gz"

# Filtering parameters
FILTER_BASE = "A,T,C,G"
SIZE = "18,30"
TRIM = False
MIN_TRIM_LENGTH = 0

# Mapping parameters
TELO_INDEX = "/nas/longleaf/home/sfrenk/proj/seq/telomere/bowtie/telomere"
GENOME_INDEX = "/nas/longleaf/home/sfrenk/proj/seq/WS251/genome/bowtie/genome"
MULTI_FLAG = "-M 1"
MISMATCH = "0"

# DB parameters
UTILDIR = ""

###############################################################################

# By default, the dataset name is taken form the read directory name.
DATASET = re.search("/([^/]*)/?$", BASEDIR).group(1)

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
		util_dir = UTILDIR
	log:
		"logs/{sample}_filter.log"
	shell:
		"python3 {params.util_dir}/small_rna_filter.py \
		-f {params.filter_base} \
		-s {params.size} \
		-t {params.trim} \
		-m {params.min_trim_length} \
		-o {output} \
		{input} > {log} 2>&1"

rule bowtie_mapping_genome:
	input:
		reads_file = "filtered/{sample}.fa"
	output: "bowtie_out/{sample}.sam"
	params:
		idx = GENOME_INDEX,
		multi_flag = MULTI_FLAG,
		mismatch = MISMATCH
	log:
		"logs/{sample}_map_genome.log"
	threads: 8
	shell:
		"bowtie -f --best --strata -S \
		-p {threads} \
		{params.multi_flag} \
		-v {params.mismatch} \
		{params.idx} \
		{input.reads_file} {output} > {log} 2>&1"

rule convert_to_bam:
	input: "bowtie_out/{sample}.sam"
	output: 
		tempfile = "genome/{sample}/temp.bam",
		bamfile = "genome/{sample}_genome.bam"
	run:
		shell("samtools view -bh -F 4 {input} > {output.tempfile}")
		shell("samtools sort -o {output.bamfile} {output.tempfile}")

rule index_bam:
	input: "genome/{sample}_genome.bam"
	output:
		"genome/{sample}_genome.bam.bai"
	shell:
		"samtools index {input}"

rule map_telo_0_mismatch:
	input: "filtered/{sample}.fa"
	output:
		un = "fasta/{sample}_unmapped_0.fa",
		bam = "telo/{sample}_0.bam"
	params:
		telo_index = TELO_INDEX,
		al = "fasta/{sample}_0_al.fa",
		m = "fasta/{sample}_0_max.fa"
	threads: 8
	log:
		"logs/{sample}_map_0.log"
	shell:
		"bowtie -f --best --strata -M 1 -S -v 0 -p {threads} --al {params.al} --un {output.un} --max {params.m} {params.telo_index} {input} | samtools view -bh -F 4 - | samtools sort -o {output.bam} - > {log} 2>&1"

rule map_telo_1_mismatch:
	input: "fasta/{sample}_unmapped_0.fa"
	output:
		un = "fasta/{sample}_unmapped_1.fa",
		bam = "telo/{sample}_1.bam"
	params:
		telo_index = TELO_INDEX,
		al = "fasta/{sample}_1_al.fa",
		m = "fasta/{sample}_1_max.fa"
	threads: 8
	log:
		"logs/{sample}_map_1.log"
	shell:
		"bowtie -f --best --strata -M 1 -S -v 1 -p {threads} --al {params.al} --un {output.un} --max {params.m} {params.telo_index} {input} | samtools view -bh -F 4 - | samtools sort -o {output.bam} - > {log} 2>&1"

rule map_telo_2_mismatch:
	input: "fasta/{sample}_unmapped_1.fa"
	output:
		un = "fasta/{sample}_unmapped_2.fa",
		bam = "telo/{sample}_2.bam"
	params:
		telo_index = TELO_INDEX,
		al = "fasta/{sample}_2_al.fa",
		m = "fasta/{sample}_2_max.fa"
	threads: 8
	log:
		"logs/{sample}_map_2.log"
	shell:
		"bowtie -f --best --strata -M 1 -S -v 2 -p {threads} --al {params.al} --un {output.un} --max {params.m} {params.telo_index} {input} | samtools view -bh -F 4 - | samtools sort -o {output.bam} - > {log} 2>&1"

rule map_telo_3_mismatch:
	input: "fasta/{sample}_unmapped_2.fa"
	output:
		un = "fasta/{sample}_unmapped_3.fa",
		bam = "telo/{sample}_2.bam"
	params:
		telo_index = TELO_INDEX,
		al = "fasta/{sample}_3_al.fa",
		m = "fasta/{sample}_3_max.fa"
	threads: 8
	log:
		"logs/{sample}_map_3.log"
	shell:
		"bowtie -f --best --strata -M 1 -S -v 3 -p {threads} --al {params.al} --un {output.un} --max {params.m} {params.telo_index} {input} | samtools view -bh -F 4 - | samtools sort -o {output.bam} - > {log} 2>&1"

rule merge_telo_reads:
	# Combine all telo-sRNA reads together
	input: "fasta/{sample}_unmapped_3.fa"
	output: "fasta/{sample}_all.fa"
	params:
		sample_name = "{sample}"
	run:
		shell('''if compgen -G "fasta/{params.sample_name}_[0-3]_*.fa" > /dev/null; then cat fasta/{params.sample_name}_[0-3]_*.fa > {output}; else touch {output}; fi''')

rule get_sample_info:
	# Get metadata:
	#	1. total number of telo-sRNAs
	#	2. total mapped reads in library
	input:
		telo_reads_file = "fasta/{sample}_all.fa",
		genome_bam = "genome/{sample}_genome.bam",
		bam_index = "genome/{sample}_genome.bam.bai"
	output:
		"sample_info/{sample}_info.txt"
	params:
		dataset = DATASET,
		sample_name = "{sample}"
	run:
		shell('''telo_reads=$(grep -v ">" {input.telo_reads_file} | wc -l); samtools view -c {input.genome_bam} | awk -v var="$telo_reads" '{{ print "{params.dataset}_{params.sample_name}\t{params.dataset}\t{params.sample_name}\t"$0"\t"var }}' > {output}''')

rule map_telo_reads_to_genome:
	input: 
		genome_bam = "genome/{sample}_genome.bam",
		fasta = "fasta/{sample}_all.fa"
	output: "telo/{sample}_genome.bam"
	params:
		idx = GENOME_INDEX
	threads: 8
	shell:
		"if [[ $(wc -l < {input.fasta}) -eq 0 ]]; then samtools view -bH {input.genome_bam} > {output}; else bowtie -f --best --strata -a -S -v 0 -p {threads} {params.idx} {input.fasta} | samtools view -bh -F 4 - | samtools sort -o {output} -; fi"

rule convert_data:
	input: "telo/{sample}_genome.bam"
	output: 
		alignment_file = "results/{sample}_alignments.txt",
		reads_file = "results/{sample}_reads.txt"
	params:
		dataset = DATASET,
		sample_name = "{sample}",
		util_dir = UTILDIR
	run:
		shell('''if [[ $(samtools view {input} | wc -l) -eq 0 ]]; then touch {output.alignment_file}; else samtools view {input} | awk -v OFS="\t" -v d={params.dataset} -v s={params.sample_name} '{{print d"_"s"_"$1,$3,$4,1-($2/16),$6}}' > {output.alignment_file}; fi''')
		shell("python3 {params.util_dir}/db/compile_data.py -d {params.dataset} -b fasta/{params.sample_name} -o {output.reads_file}")

rule store_results:
	# Load library data and metadata into SQLite3 database
	input:
		alignment_file = "results/{sample}_alignments.txt",
		reads_file = "results/{sample}_reads.txt",
		sample_info = "sample_info/{sample}_info.txt"
	output:
		"results/{sample}.db"
	params:
		util_dir = UTILDIR
	shell:
		"python3 {params.util_dir}/db/ngs_store.py -d {output} -s {input.sample_info} -r {input.reads_file} -a {input.alignment_file}"

rule merge_db:
	# Add all library databases into one database for the dataset
	input:
		expand("results/{sample}.db", sample = SAMPLES)
	output: "results/" + DATASET + ".db"
	params:
		util_dir = UTILDIR
	shell:
		"python3 {params.util_dir}/db/merge_db.py -o {output} {input}"
