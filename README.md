# Pipeline for analysis of telomere-related small RNAs (telo-sRNAs)

## Usage

1. Use setup_working_dir.sh to set up the working directory for analysis:

```
setup_working_dir.sh -p srna_telo -d <directory containing reads>
```

2. Run the snakemake script

```
bash run_snakemake.sh
```

3. Merge results databases

```
python3 merge_db.py -o <merged db name> <list of dbs to merge>
```