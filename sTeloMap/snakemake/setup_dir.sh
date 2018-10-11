#!/usr/bin/env bash

# Hard variables

# Directory containing Snakemake and cluster.json files
# This should have been set by install.sh

snakefile_dir=""

usage="
Create directory with Snakemake files required for sTeloMap. Usage:

setup_dir.sh <working directory (default: current directory)>
"

help_args=("-h" "--help")

if [[ -z "$1" ]]; then
    dir="."
elif [[ ${help_args[@]} =~ "$1" ]]; then
	echo "$usage"
	exit 0
else
	dir="$1"
fi

if [[ $dir == $snakefile_dir ]]; then
	echo "ERROR: Working directory must be in different directory to snakefile"
	exit 1
fi

# Copy over the snakefile
snakefile="srna_telo.Snakefile"
cp ${snakefile_dir}/${snakefile} ./${snakefile}

# Edit base directory in Snakefile
# Remove trailing "/" from dir if it's there
dir_name="$(echo $dir |sed -r 's/\/$//')"
sed -r -i -e "s,^BASEDIR.*,BASEDIR = \"${dir_name}\"," "$snakefile"

# Determine file extension
extension="$(ls $dir | grep -Eo "\.[^/]+(\.gz)?$" | sort | uniq)"

# Check if there are multiple file extensions in the same directory
ext_count="$(echo $extension | wc -l)"

if [[ ext_count == 0 ]]; then
	echo "ERROR: Directory is empty!"
elif [[ ext_count != 1 ]]; then
	echo "WARNING: Multiple file extensions found: using .fastq.gz"
	extension=".fastq.gz"
fi

# Edit extension and utils_dir in Snakefile
extension="\"${extension}\""
sed -i -e "s|^EXTENSION.*|EXTENSION = ${extension}|g" "$snakefile"
utils_dir="${snakefile_dir%/snakemake}"
utils_dir="\"${utils_dir}\""
sed -i -e "s|^UTILS_DIR.*|UTILS_DIR = ${utils_dir}|g" "$snakefile"

# Create Snakmake command script
printf "#!/usr/bin/bash\n" > "run_snakemake.sh"
printf "#SBATCH -t 2-0\n\n" >> "run_snakemake.sh"
printf "module add python\n\n" >> "run_snakemake.sh"
printf "snakemake -s $snakefile --rerun-incomplete --cluster-config ${snakefile_dir}/cluster.json -j 100 --cluster \"sbatch -n {cluster.n} -N {cluster.N} -t {cluster.time}\"\n" >> "run_snakemake.sh"
