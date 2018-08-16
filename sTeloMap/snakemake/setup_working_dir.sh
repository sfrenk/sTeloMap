#!/usr/bin/env bash

# Hard variables

# Directories containing scripts
# Determine these by checking the location of this script
#snakedir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
snakedir="$(dirname "$(readlink -f "$0")")"
utildir="$(echo $snakedir | sed -r 's/snakemake$/util/')"

# Modules to be used in pipeline
modules="anaconda python bowtie/1.1.2 samtools subread"

##############################################################################

usage="Create directory with Snakemake files required for pipeline \n\n setup_dir -p <pipeline> -d <directory> \n\n pipelines: srna_telo, remap"

pipeline=""

if [ -z "$1" ]; then
    echo "$usage"
    exit
fi

while [[ $# > 0 ]]
do
    key="$1"
    case $key in
    	-p|--pipeline)
		pipeline="$2"
		shift
		;;
        -d|--dir)
        dir="$2"
        shift
        ;;
        -h|--help)
		printf "$usage"
		exit
		;;
    esac
    shift
done


if [[ ! -d $dir ]]; then
	echo "ERROR: Invalid directory"
	exit 1
fi

if [[ $pipeline == "" ]]; then
	echo "ERROR: Please select pipeline"
	exit 1
fi

# Determine pipeline file
case $pipeline in
	"srna_telo")
	snakefile="srna_telo.Snakefile"
	;;
	"remap")
	snakefile="remap.Snakefile"
	;;
	*)
	echo "ERROR: Invalid pipeline. Please select one of the following: srna_telo, remap"
	exit 1
	;;
esac

# Copy over the snakefile
cp ${snakedir}/${snakefile} .

# Edit base directory in Snakefile
base="$(basename ${dir})"
sed -r -i -e "s,^BASEDIR.*,BASEDIR = \"${dir}\"," "$snakefile"

# Determine file extension
extension="$(ls $dir | grep -Eo "\.[^/]+" | sort | uniq)"

# Check if there are multiple file extensions in the same directory
ext_count="$(echo $extension | wc -l)"

if [[ ext_count == 0 ]]; then
	echo "ERROR: Directory is empty!"
elif [[ ext_count != 1 ]]; then
	echo "WARNING: Multiple file extensions found: using .fastq.gz"
	extension=".fastq.gz"
fi

# Edit extension in Snakefile
extension="\"${extension}\""
sed -i -r -e "s@^EXTENSION.*@EXTENSION = ${extension}@g" "$snakefile"

# Edit util dir location in Snakefile
sed -i -r -e "s@^UTILDIR.*@UTILDIR = \"${utildir}\"@g" "$snakefile"

# Create Snakmake command script
printf "#!/usr/bin/bash\n\n" > "run_snakemake.sh"
printf "module add $modules\n\n" >> "run_snakemake.sh"
printf "snakemake -s $snakefile --cluster-config ${snakedir}/cluster.json -j 100 --cluster \"sbatch -n {cluster.n} -N {cluster.N} -t {cluster.time}\"\n" >> "run_snakemake.sh"
