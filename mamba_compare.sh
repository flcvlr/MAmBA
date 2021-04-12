#!/bin/bash
# -t target_sample
# -c control_sample
# -m minimum_coverage
# -f minimum_frequency
# -s species
# -o output_directory
# -j cores
# -u short neme for target
# -d short name for control
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

show_help () { cat $SCRIPT_DIR/readme_compare.txt; }

while getopts "t:c:m:f:s:o:j:d:u:" opt; do
    case "$opt" in
    h|\?)
        show_help
        exit 0
        ;;
    t)  target=$OPTARG
        ;;
    c)  control=$OPTARG
        ;;
    f)  min_freq=$OPTARG
        ;;
    m)  min_cov=$OPTARG
        ;;
    s)  species=$OPTARG
        ;;
    o)  output_dir=$OPTARG
        ;;
    j)  cores=$OPTARG
        ;;
    d)  cont_name=$OPTARG
        ;;
    u)  target_name=$OPTARG
        ;;

    esac
done

shift $((OPTIND-1))

[ "${1:-}" = "--" ] && shift



if [ -d $output_dir ]; then
	echo -e "\n\n\t"
	read -p "Output directory exists already. Overwrite? " -n 1 -r
	if [[ $REPLY =~ ^[Yy]$ ]]
	then
		echo -e "\n\n\tOverwriting previous output"
		rm $output_dir/*.png
	else
		exit
	fi
else
	mkdir $output_dir $output_dir/logs 
fi	

intersectBed -a $target/hairpin_subbed.bed -b $control/hairpin_subbed.bed -wo -s | awk -v MC=$min_cov -v MF=$min_freq -v OFS='\t' '{if ( (($9 > MC) && ($18 > MC)) && (($7 > MF) || ($16 > MF)) ) print $1,$2,$3,$4,$8,$9,$17,$18}' > $output_dir/meth_data_to_test
Rscript --vanilla $SCRIPT_DIR/R_scripts/plot_differential.R $target $control $target_name $cont_name $SCRIPT_DIR/$species $output_dir
