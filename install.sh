#!/usr/bin/bash

# Find snakefile directory
script_dir="$(pwd)/sTeloMap"
snakefile_dir="${script_dir}/snakemake"

echo "Snakefile directory: $snakefile_dir"

# Edit setup_dir.sh to contain correct snakefile directory
#snakefile_dir_sed="\"${snakefile_dir}\""
#echo $snakefile_dir_sed
sed -i -e "s|^snakefile_dir.*|snakefile_dir=${snakefile_dir}|g" "${snakefile_dir}/setup_dir.sh"

echo "Done!"
