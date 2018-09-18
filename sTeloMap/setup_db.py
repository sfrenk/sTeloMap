#!/usr/bin/env python3

import os
import sys
import sqlite3
import argparse

# Setup a new sRNA-telo sqlite3 database

def setup_db(db_name):

	if os.path.exists(db_name):
		sys.exit('ERROR: ' + db_name + ' already exists!')

	print("Creating new database: " + db_name)

	conn = sqlite3.connect(db_name)

	c = conn.cursor()

	## Create tables ##

	# sample_info
	c.execute("CREATE TABLE sample_info (sample_id text primary key, dataset text, name text, total_mapped integer, telo_reads integer)")
	c.execute("CREATE TABLE reads (read_id text primary key, sample_id text, sequence text, mismatches integer, foreign key(sample_id) references sample_info(sample_id))")
	c.execute("CREATE TABLE alignments (alignment_id integer primary key autoincrement, read_id text, chrom text, pos integer, strand integer, n_hits integer, foreign key(read_id) references reads(read_id))")

	# Save changes
	conn.commit()
	
	# Close connection
	conn.close()

def main():

	parser = argparse.ArgumentParser()

	parser.add_argument("db_name", help = "Name of database file (eg. sample1.db)")

	args = parser.parse_args()

	setup_db(args.db_name)

if __name__ == "__main__":
	main()
