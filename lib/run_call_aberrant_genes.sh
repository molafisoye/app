#!/bin/bash
usage="$(basename "$0") [-h] -i INPUT_GENE_COUNT -o OUTPUT_DIR -p PROJECT_NAME -t SUBSAMPLE_VALUE 

where:
    -h  show this help text and exit
    -i	input gene count table
    -o	output directory
    -p	project name
    -t	target to subsample to"
    
while getopts "hi:o:p:t:" OPTION
do
	case $OPTION in
		h)
    			echo "$usage"
			exit
			;;
		i)
			input=$OPTARG
			;;
		o)	
			output=$OPTARG
			;;
		p)
			project=$OPTARG
			;;
		t) 
			target=$OPTARG
			;;
		:) 	printf "missing argument for -%s\n" "$OPTARG" >&2
       			echo "$usage" >&2
       			exit 1
       			;;
   		\?) 	printf "illegal option: -%s\n" "$OPTARG" >&2
       			echo "$usage" >&2
       			exit 1
       			;;	

	esac
done

# Absolute path to this script. /home/user/bin/foo.sh
SCRIPT=$(readlink -f $0)
# Absolute path this script is in. /home/user/bin
SCRIPTPATH=`dirname $SCRIPT`

echo "Subsampling gene count table"

CMD_1="Rscript $SCRIPTPATH/src/subsample_V1.R -i $input -o "${output}/${project}_subsampled_${target}.csv" -t $target"

$CMD_1

echo "starting normalizer with input ${output}/${project}_subsampled_${target}.csv removing previous normalized.txt"

node --experimental-modules $SCRIPTPATH/src/normalizer.mjs "/app/${project}_subsampled_${target}.csv"


echo "Scaling and transforming data..."
CMD_2="Rscript $SCRIPTPATH/src/scale_and_transform_V1.R --input=/app/normalized.txt -o "${output}/${project}_normalized.csv""

$CMD_2

echo "Detecting aberrantly expressed genes..."
CMD_3="Rscript $SCRIPTPATH/src/find_abberants_V1_01.R -i /app/${project}_normalized.csv -o $output -p $project"

$CMD_3

## Move normalized.txt to the output directory (so it doesn't get overwritten by the next run)
#if [ -f "normalized.txt" ] && [ ! -f "${output}/${project}_normalized.txt" ]
#then
#	mv "normalized.txt" "${output}/${project}_normalized.txt"
#fi
#
#echo 'Analysis complete!'
#
