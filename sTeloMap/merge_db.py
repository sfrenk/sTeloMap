#!/usr/bin/env python3

import sqlite3
import argparse
import setup_db
import os
import sys
import pandas as pd

# Merge telo-sRNA databases

def merge_db(db_output_name, databases):

	'''
	Merges multiple dbs into a single output db
	'''

	# Setup master db (this will be the output)
	if not os.path.exists(db_output_name):
		setup_db.setup_db(db_output_name)
	else:
		sys.exit("ERROR: " + db_output_name + " already exists!")

	db_master = sqlite3.connect(db_output_name)
	db_master_c = db_master.cursor()

	for database in databases:
	
		# Connect to input db
		db_input = sqlite3.connect(database)
		db_input_c = db_input.cursor()

		for table in ["sample_info", "reads", "alignments"]:

			# Convert table from input db into a pandas dataframe
			rows = pd.read_sql_query("SELECT * FROM " + table + ";", db_input)

			if table == "alignments":
				# Get rid of the index for alignments as this is an autoincramenting primary key
				rows = rows.drop("alignment_id", 1)		
			
			# Append input table to output table
			rows.to_sql(table, db_master, if_exists = "append", index = False)

	db_input.close()

	db_master.commit()
	db_master.close()

def main():

	parser = argparse.ArgumentParser()

	parser.add_argument("dbs", help = "List of databases to merge", nargs = "+")
	parser.add_argument("-o", "--output", help = "output database filename")

	args = parser.parse_args()

	merge_db(args.output, args.dbs)

if __name__ == "__main__":
	main()
