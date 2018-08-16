#!/usr/bin/env python3

import pytest
import sqlite3
from sTeloMap.setup_db import setup_db
import os

# Setup database
if os.path.exists("test.db"):
	os.remove("test.db")

setup_db("test.db")

# Check the tables
con = sqlite3.connect("test.db")
cursor = con.cursor()
cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
tables = cursor.fetchall()
con.close()
os.remove("test.db")

def test_setup_db():
	assert tables == [('sample_info',), ('reads',), ('alignments',), ('sqlite_sequence',)]