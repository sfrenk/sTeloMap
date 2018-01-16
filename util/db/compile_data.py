#!/usr/bin/env python3

import os
from Bio import SeqIO
import subprocess
import argparse

def compile_reads_data(reads_base, reads_output_file, dataset):

	''' Combines fasta files of reads mapping to telomeres into a single table file '''
	
	with open(reads_output_file, "w") as f:
		
		for i in range(4):
			
			for j in ("al", "max"):
				reads_input_file = reads_base + "_" + str(i) + "_" + j + ".fa"
				
				if os.path.exists(reads_input_file):
					
					for record in SeqIO.parse(reads_input_file, "fasta"):
						f.write(dataset + "_" + os.path.basename(reads_base) + "_" + str(record.id) + "\t" + dataset + "_" + os.path.basename(reads_base) + "\t" + str(record.seq) + "\t" + str(i) + "\n")

def main():

	parser = argparse.ArgumentParser()

	parser.add_argument("-d", "--dataset", help = "Name of dataset")
	parser.add_argument("-b", "--base", help = "base name for fasta files (eg. ./fasta/sample1)")
	parser.add_argument("-o", "--output", help = "output filename")

	args = parser.parse_args()

	compile_reads_data(args.base, args.output, args.dataset)

if __name__ == "__main__":
	main()
