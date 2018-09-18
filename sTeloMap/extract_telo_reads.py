#!/usr/bin/env python3

import pysam
import argparse
import os

parser = argparse.ArgumentParser(description = "Process telomere reads")

parser.add_argument("-t", "--telo_bam", help = "bam file containing telomere alignments", required = True)
parser.add_argument("-g", "--genome_bam", help = "bam file containing genome alignments (output form Butter)", required = True)
parser.add_argument("-o", "--output", help = "output basename", required = True, type = str)
parser.add_argument("-d", "--dataset", help = "dataset name", required = True, type = str)
parser.add_argument("-s", "--sample", help = "sample name (by default, this is extracted from the bam file name)", default = "default")

args = parser.parse_args()

print("Extracting telomeric reads...")

if args.sample == "default":
	sample_name = os.path.basename(args.input_file).replace(".bam", "")
else:
	sample_name = args.sample


## Read info ##
read_ids = []
# Open input file
telo_bamfile = pysam.AlignmentFile(args.telo_bam, "rb")

with open(args.output + "_reads.txt", "w") as f:
	for read in telo_bamfile.fetch():

		# Get number of mismatches
		mm = int(read.get_tag("XM"))
		if mm <= 3:
			read_ids.append(read.query_name)
			f.write("\t".join([args.dataset + "_" + sample_name + "_" + read.query_name, args.dataset + "_" + sample_name, read.seq, str(mm)]) + "\n")

telo_bamfile.close()

## Alignment info ##
genome_bamfile = pysam.AlignmentFile(args.genome_bam, "rb")

with open(args.output + "_alignments.txt", "w") as f:
	for read in genome_bamfile.fetch():

		if read.query_name in read_ids:
			
			# Get chromosome
			chrom = genome_bamfile.get_reference_name(read.reference_id)

			# Get number of hits
			nh = read.get_tag("XX")
			# Get strand
			if read.is_reverse:
				strand = 0
			else:
				strand = 1

			f.write("\t".join([args.dataset + "_" + sample_name + "_" + read.query_name, chrom, str(read.reference_start), str(strand), str(nh)]) + "\n")

telo_bamfile.close()
