import sqlite3

conn = sqlite3.connect('/home/sfrenk/Documents/Projects/charlie_paper/srna_telo/pipeline/db/test/results/cde_1_input.db')
#de_1_input.db
c = conn.cursor()

c.execute("SELECT * FROM reads")
print([description[0] for description in c.description])

#print([description[0] for description in c.description])

conn.close()