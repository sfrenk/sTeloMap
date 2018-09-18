#!/usr/bin/env python3

import sqlite3
import pandas as pd
import setup_db
import os
import argparse

# Store telo-sRNA table file as sqlite3 database

def sql_store(srna_db, sample_file, reads_file, alignments_file):

	''' Puts data into sqlite3 db'''
	
	# Open db connection
	conn = sqlite3.connect(srna_db)
	c = conn.cursor()

	# Store sample info
	sample_info = pd.read_csv(sample_file, sep = "\t", header = None, names = ["sample_id", "dataset", "name", "total_mapped", "telo_reads"])

	sample_info.to_sql("sample_info", conn, if_exists = 'append', index = False)

	# Store reads
	reads = pd.read_csv(reads_file, sep = "\t", header = None, names = ["read_id", "sample_id", "sequence", "mismatches"])

	#reads.read_id = reads.sample_id + "_" + reads.read_id

	reads.to_sql("reads", conn, if_exists = 'append', index = False)

	# Store alignments
	alignments = pd.read_csv(alignments_file, sep = "\t", header = None, names = ["read_id", "chrom", "pos", "strand", "n_hits"])

	alignments.to_sql("alignments", conn, if_exists = 'append', index = False)

	# Save changes
	conn.commit()

	# Close db
	conn.close()

def main():

	parser = argparse.ArgumentParser()

	parser.add_argument("-d", "--db", help = "Name of sqlite database for storage")
	parser.add_argument("-s", "--sample_file", help = "sample_info file")
	parser.add_argument("-r", "--reads_file", help = "reads file")
	parser.add_argument("-a", "--alignments_file", help = "alignments file")

	args = parser.parse_args()

	if not os.path.exists(args.db):
		setup_db.setup_db(args.db)

	sql_store(args.db, args.sample_file, args.reads_file, args.alignments_file)

if __name__ == "__main__":
	main()
	