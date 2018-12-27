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
parser.add_argument("-b", "--bowtie", help = "reads were mapped with bowtie rather than butter, therefore XX, XY and XZ flags are not present and will be set to NA", action = "store_true", default = False)

args = parser.parse_args()

print("Extracting telomeric reads...")

# Extract sample name from bam file name unless specifically given using the --sample argument
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
		mm = int(read.get_tag("NM"))
		if mm <= 3:
			read_ids.append(read.query_name)
			f.write("\t".join([args.dataset + "_" + sample_name + "_" + read.query_name, args.dataset + "_" + sample_name, read.get_forward_sequence(), str(mm)]) + "\n")

telo_bamfile.close()

## Alignment info ##
genome_bamfile = pysam.AlignmentFile(args.genome_bam, "rb")

# Alignments are written to both a text file and a bam file
output_bamfile = pysam.AlignmentFile(args.output + "_telo.bam","wb", template = genome_bamfile)

with open(args.output + "_alignments.txt", "w") as f:
	for read in genome_bamfile.fetch():

		if read.query_name in read_ids:

			# Write read to bam file
			output_bamfile.write(read)
			
			# Get chromosome
			chrom = genome_bamfile.get_reference_name(read.reference_id)

			if args.bowtie:
				nhits_tag = "NA"
				placement_tag = "NA"
				prob_tag = "NA"
			else:
				# Get number of hits
				nhits_tag = read.get_tag("XX")

				# Get placement method:
				placement_tag = read.get_tag("XY")

				# Get probability of correct mapping
				prob_tag = read.get_tag("XZ")

			# Get strand
			if read.is_reverse:
				strand = 0
			else:
				strand = 1

			f.write("\t".join([args.dataset + "_" + sample_name + "_" + read.query_name, chrom, str(read.reference_start), str(strand), str(nhits_tag), str(placement_tag), str(prob_tag)]) + "\n")

genome_bamfile.close()
output_bamfile.close()
pysam.index(args.output + "_telo.bam", args.output + "_telo.bam.bai")
