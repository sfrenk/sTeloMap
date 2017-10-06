#!/usr/bin/env bash

# Variables
script_dir="/nas/longleaf/home/sfrenk/scripts/charlie_paper/snakemake"

if [[ -z "$1" ]]; then
	echo "Please specify directory!"
	exit 1
fi

# Make subdirectory
srna_dir="$1"

base="$(basename ${srna_dir})"

mkdir "./${base}"

# Copy over the necessary files
cp "${script_dir}/Snakefile" ./${base}
cp "${script_dir}/cluster.json" ./${base}
cp "${script_dir}/config.json" ./${base}


# Edit base directory in config file
sed -i -r -e "s,\"basedir\": \".*\",\"basedir\": \"${srna_dir}\"," "./${base}/config.json"

# Determine file extension
extension="$(ls $srna_dir | grep -Eo "\.[^/]+" | sort | uniq)"

if [[ "${#extension[@]}" != 1 ]]; then
	echo "ERROR: All files must have the same extension"

else

	# Edit extension in config file
	sed -i -r -e "s,\"extension\": \".*\",\"extension\": \"${extension}\"," "./${base}/config.json"
fi
