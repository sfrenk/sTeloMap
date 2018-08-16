# Pipeline for analysis of telomere-related small RNAs (telo-sRNAs)

This pipeline searches for reads that map to telomeres with up to three mismatches. Results are stored in an sqlite3 database which contains sequence data and sample metadata.

## Installation

```
git clone https://github.com/sfrenk/sTeloMap.git

bash install.sh
```


## Usage

1. Use setup_working_dir.sh to set up the working directory for analysis:

```
setup_dir.sh -d <directory containing reads>
```

2. Run the snakemake script

```
bash run_snakemake.sh
```

3. Merge results databases

```
python3 merge_db.py -o <merged db name> <list of dbs to merge>
```