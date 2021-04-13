#!/bin/bash

### command line options reminder
#	-a adaptor sequence (default set to "no", will cause adaptor trimming skip)
#	-s species (hsa|mmu)
#	-o output_directory (default: mamba_out.$infile)
#	-h|? show help 
#	-f fastq.gz input file
#	-i bowtie2 index for the genome of interest (assuming the corresponding fasta is in the same directory with the same basename. If missing, will be created)
#	-j number_of_cores
#	-n short_name to use in images
#	-l log/linear y axis in pictures
#	-C color to use in pictures (default: darkgreen)
#
#
#


# Initialize our own variables:

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
OPTIND=1         # Reset in case getopts has been used previously in the shell.
output_dir="mamba_out"
ADAPTER="no"
control_spike="no"
show_help ()  { cat $SCRIPT_DIR/readme.txt; }

### options parsing ###

while getopts "h?a:o:f:i:s:c:j:n:l:C:" opt; do
    case "$opt" in
    h|\?)
        show_help
        exit 0
        ;;
    a)  ADAPTER=$OPTARG
        ;;
    o)  output_dir=$OPTARG
        ;;
    f)  infile=$OPTARG
        ;;
    i)  index=$OPTARG
        ;;
    s)  species=$OPTARG
        ;;
    c)  control_spike=$OPTARG
        ;;
    j)  cores=$OPTARG
        ;;
    n)  short_name=$OPTARG
        ;;
    l)  log=$OPTARG
        ;;
    C)  color=$OPTARG
        ;;
    esac
done

shift $((OPTIND-1))

[ "${1:-}" = "--" ] && shift


#### Check for files and directories ####


if [ -z $color ]; then
	color="darkgreen";
fi


if [ -z $log ]; then
	log="log";
fi

if [ ! -s $infile ]; then
	1>&2 echo -e "\n\n\tError: Failed to launch MAmBA: input $infile does not exist or is empty. Please check filename and/or path.\n\n";
	exit;
fi	

if [ -d $output_dir ]; then
	echo -e "\n\n\t"
	read -p "Output directory exists already. Overwrite? " -n 1 -r
	if [[ $REPLY =~ ^[Yy]$ ]]
	then
		echo -e "\n\n\tOverwriting previous output"
		rm $output_dir/data_for_R/*
		rm $output_dir/images/*
		rm $output_dir/logs/*
		rm $output_dir/total_transcriptome_out/logs/*
		rm $output_dir/total_transcriptome_out/*
		

	else
		exit
	fi
else
	mkdir $output_dir $output_dir/logs
fi	

#### check for required software
exists()
{
	command -v $1 >/dev/null 2>&1 || { echo -e >&2 "\n\tI require $1 but it's not installed.  Aborting."; show_help; exit 0;}
}

exists tophat2
exists bowtie2
exists bedtools
exists samtools
exists Rscript
exists parallel
exists pigz

####  Check conversion efficiency using spike in
if [ ${#control_spike} -gt 10 ]
then
	perl $SCRIPT_DIR/perl_scripts/spike_in_count.pl $infile $control_spike $output_dir
	#read -p "Based on this C --> U conversion efficiency, do you wish to continue the analysis?(y/n)" -n 1 -r
	#if [[ $REPLY =~ ^[Nn]$ ]]
	#then
	#	echo -e "\n\n\texiting...\n\n"
	#	exit
	#fi

else
	1>&2 echo -e "\n\tskipping spike-in control, no control or too short control spike in small RNA was provided with -c option.\n\n"
fi


#### Adapter trimming
if [ ${#ADAPTER} -gt 10 ]
then
	trimmed=`basename $infile .gz`


	1>&2 echo -e "\n\n\tTrimming $infile using adapter $ADAPTER. Trimmed file will be saved as ${output_dir}/$trimmed.gz.\n\tNo quality trimming will be applied. Reads shorter than 18 nt after trimming will be discarded.\n\tPlease check trimming log in ${output_dir}/logs/$trimmed.trim.log\n"

	zcat $infile | cutadapt -j $cores -a $ADAPTER -m 18 - --trim-n 2> ${output_dir}/logs/$trimmed.trim.log | pigz >  ${output_dir}/$trimmed.gz
	infile="${output_dir}/$trimmed.gz"

else
	trimmed=`basename $infile .gz`


	1>&2 echo -e "\tTrimming terminal Ns in $infile reads. Trimmed file will be saved as ${output_dir}/$trimmed.gz.\n\tNo quality trimming will be applied. Reads shorter than 18 nt after trimming will be discarded.\n\tPlease check trimming log in ${output_dir}/logs/$trimmed.trim.log\n"

	zcat $infile | cutadapt -j $cores -m 18 - --trim-n 2> ${output_dir}/logs/$trimmed.trim.log | pigz >  ${output_dir}/$trimmed.gz
	infile="${output_dir}/$trimmed.gz"


fi

zcat $infile | perl -ne 'if (($a % 4) == 0) {$n=$n+1; $_="@".$n.$_;} $a +=1; print $_;' > ${output_dir}/trimmed.fastq
sed '2~4s/C/T/g;' ${output_dir}/trimmed.fastq > ${output_dir}/trimmed_CT_converted.fq

#### check C --> U conversion efficiency on 18S and 28S rRNAs

mkdir $output_dir/C_T_efficiency_18S
echo -e "\n18S alignement stats:\n" > $output_dir/logs/alignment.log
bowtie2 -u 300000 -q -N 0 -L 18 -p 8 --score-min L,0,-0.2 --norc -x $SCRIPT_DIR/$species/18S_C_T_converted -U $output_dir/trimmed_CT_converted.fq 2>> $output_dir/logs/alignment.log| samtools view -Sb -  > $output_dir/C_T_efficiency_18S/18s_C_T.bam 

samtools view -H $output_dir/C_T_efficiency_18S/18s_C_T.bam > $output_dir/C_T_efficiency_18S/header.txt

samtools view -F 256 -F 4 $output_dir/C_T_efficiency_18S/18s_C_T.bam | awk '{if (($6 !~ /[I,D]/) &&  ($6 !~ /N.+N/)) print}' - | LC_ALL=C sort -S 2G -k 1,1n -t '@' - | cat $output_dir/C_T_efficiency_18S/header.txt - | samtools view -Sb - > $output_dir/C_T_efficiency_18S/18s_C_T_clean.bam
bedtools bamtobed -i $output_dir/C_T_efficiency_18S/18s_C_T_clean.bam | bedtools getfasta -name -split -s -fi $SCRIPT_DIR/$species/18S.fa -bed - -fo $output_dir/C_T_efficiency_18S/18s_C_T.fa

perl $SCRIPT_DIR/perl_scripts/methylation_call_and_coverage.pl $output_dir/C_T_efficiency_18S/18s_C_T_clean.bam $output_dir/trimmed.fastq $output_dir/C_T_efficiency_18S/18s_C_T.fa $output_dir/C_T_efficiency_18S/
bedtools merge -i $output_dir/C_T_efficiency_18S/methylated_C_sorted.bed.gz -d -1 -s -c 6,4,5 -o first,count_distinct,sum | awk -v OFS='\t' '{print $1,$2,$3,$(NF-1),$(NF-0),$(NF-2)}' | pigz > $output_dir/C_T_efficiency_18S/methylated_C_summary_sorted.bed.gz
bedtools merge -i $output_dir/C_T_efficiency_18S/total_C_sorted.bed.gz  -d -1 -s -c 6,4,5 -o first,count_distinct,sum | awk -v OFS='\t' '{print $1,$2,$3,$(NF-1),$(NF-0),$(NF-2)}' | pigz > $output_dir/C_T_efficiency_18S/total_C_summary_sorted.bed.gz
bedtools intersect -sorted -a $output_dir/C_T_efficiency_18S/total_C_summary_sorted.bed.gz -b $output_dir/C_T_efficiency_18S/methylated_C_summary_sorted.bed.gz -s -wao | awk -v OFS='\t' '{if ($13 =="0") print $1,$2,$3,".",".",$6,"0","0",$5; else print $1,$2,$3,".",".",$6,$11/$5,$11,$5}' > $output_dir/C_T_efficiency_18S/total_transcriptome_sorted.bed
echo -e -n "\n18S alignement conversion efficiency:" >> $output_dir/logs/alignment.log
echo -e -n "\n18S alignement conversion efficiency:"
cat $output_dir/C_T_efficiency_18S/total_transcriptome_sorted.bed | Rscript --vanilla -e 'aa <- read.table("stdin"); 1-(sum(aa[,8])/sum(aa[,9]));' | cut -f 2 -d ' '
cat $output_dir/C_T_efficiency_18S/total_transcriptome_sorted.bed | Rscript --vanilla -e 'aa <- read.table("stdin"); 1-(sum(aa[,8])/sum(aa[,9]));' >> $output_dir/logs/alignment.log


mkdir $output_dir/C_T_efficiency_28S
echo -e "\n28S alignement stats:\n" >> $output_dir/logs/alignment.log
bowtie2 -u 300000 -q -N 0 -L 18 -p 8 --score-min L,0,-0.2 --norc -x $SCRIPT_DIR/$species/28S_C_T_converted -U $output_dir/trimmed_CT_converted.fq 2>> $output_dir/logs/alignment.log | samtools view -Sb - > $output_dir/C_T_efficiency_28S/28s_C_T.bam 

samtools view -H $output_dir/C_T_efficiency_28S/28s_C_T.bam > $output_dir/C_T_efficiency_28S/header.txt

samtools view -F 256 -F 4 $output_dir/C_T_efficiency_28S/28s_C_T.bam | awk '{if (($6 !~ /[I,D]/) &&  ($6 !~ /N.+N/)) print}' - | LC_ALL=C sort -S 2G -k 1,1n -t '@' - | cat $output_dir/C_T_efficiency_28S/header.txt - | samtools view -Sb - > $output_dir/C_T_efficiency_28S/28s_C_T_clean.bam
bedtools bamtobed -i $output_dir/C_T_efficiency_28S/28s_C_T_clean.bam | bedtools getfasta -name -split -s -fi $SCRIPT_DIR/$species/28S.fa -bed - -fo $output_dir/C_T_efficiency_28S/28s_C_T.fa

perl $SCRIPT_DIR/perl_scripts/methylation_call_and_coverage.pl $output_dir/C_T_efficiency_28S/28s_C_T_clean.bam $output_dir/trimmed.fastq $output_dir/C_T_efficiency_28S/28s_C_T.fa $output_dir/C_T_efficiency_28S/
bedtools merge -i $output_dir/C_T_efficiency_28S/methylated_C_sorted.bed.gz  -d -1 -s -c 6,4,5 -o first,count_distinct,sum | awk -v OFS='\t' '{print $1,$2,$3,$(NF-1),$(NF-0),$(NF-2) }' | pigz > $output_dir/C_T_efficiency_28S/methylated_C_summary_sorted.bed.gz
bedtools merge -i $output_dir/C_T_efficiency_28S/total_C_sorted.bed.gz  -d -1 -s -c 6,4,5 -o first,count_distinct,sum | awk -v OFS='\t' '{print $1,$2,$3,$(NF-1),$(NF-0),$(NF-2) }' | pigz > $output_dir/C_T_efficiency_28S/total_C_summary_sorted.bed.gz
bedtools intersect -sorted -a $output_dir/C_T_efficiency_28S/total_C_summary_sorted.bed.gz -b $output_dir/C_T_efficiency_28S/methylated_C_summary_sorted.bed.gz -s -wao | awk -v OFS='\t' '{if ($13 =="0") print $1,$2,$3,".",".",$6,"0","0",$5; else print $1,$2,$3,".",".",$6,$11/$5,$11,$5}' > $output_dir/C_T_efficiency_28S/total_transcriptome_sorted.bed
echo -e -n "\n28S alignement conversion efficiency:" >> $output_dir/logs/alignment.log
echo -e -n "\n28S alignement conversion efficiency:"
cat $output_dir/C_T_efficiency_28S/total_transcriptome_sorted.bed | Rscript --vanilla -e 'aa <- read.table("stdin"); 1-(sum(aa[,8])/sum(aa[,9]));' | cut -f2 -d ' '
cat $output_dir/C_T_efficiency_28S/total_transcriptome_sorted.bed | Rscript --vanilla -e 'aa <- read.table("stdin"); 1-(sum(aa[,8])/sum(aa[,9]));' >> $output_dir/logs/alignment.log



#### remove tRNA and rRNA fragments with bowtie2



1>&2 echo -e "\nrRNA/tRNA alignment (bowtie2)" >> ${output_dir}/logs/alignment.log
bowtie2 -q -N 0 -L 18 -p $cores --score-min L,0,-0.2 --norc --un-gz ${output_dir}/unmapped_on_rRNA_tRNA_CT_converted.fq.gz -x $SCRIPT_DIR/$species/rRNA_tRNA_converted/rRNAs_tRNAs_converted -U $output_dir/trimmed_CT_converted.fq > /dev/null 2>> ${output_dir}/logs/alignment.log

rm ${output_dir}/trimmed_CT_converted.fq

#### Alignment to transcriptome with tophat2

tophat2 -o ${output_dir}/total_transcriptome_out --b2-N 0 --b2-L 18 -p $cores --library-type fr-secondstrand --no-sort-bam --transcriptome-index=$SCRIPT_DIR/$species/total_transcriptome_converted/annotation -T $index ${output_dir}/unmapped_on_rRNA_tRNA_CT_converted.fq.gz

echo -e "\ntotal transcriptome alignment (tophat2)" >> ${output_dir}/logs/alignment.log
cat ${output_dir}/total_transcriptome_out/align_summary.txt >> ${output_dir}/logs/alignment.log


#### Identify methylated cytosines

##	Preparation of files for methylation_call.pl script: 
##
##	1) extraction of sorted primary aligned reads in .sam file 

##  save the header for reuse
samtools view -H ${output_dir}/total_transcriptome_out/accepted_hits.bam > ${output_dir}/header.txt 

##  extracts only primary alignments (-F 256) without indels and containing no more than 1 splicing junction (use awk to filter on CIGAR) sort by readname, put back header in place and save in bam format
samtools view -F 256 ${output_dir}/total_transcriptome_out/accepted_hits.bam | awk '{if (($6 !~ /[I,D]/) && ($21=="NH:i:1") && ($6 !~ /N.+N/)) print}' - | LC_ALL=C sort -S 2G -k 1,1n -t '@' - | cat ${output_dir}/header.txt - | samtools view -Sb - > ${output_dir}/total_transcriptome_primary_sorted.bam 

##	
##	2) extraction of references in sorted .fa file
##  convert bam file to bed (splitting on CIGAR N operation), get fasta sequence from reference genome. This is done preserving sorting of alignements 
if [ ! -s $index.fa ]
then
	echo "creating a fasta file corresponding to the bowtie2 index";
	bowtie2-inspect $index > $output_dir/genome.fa;
	bedtools bamtobed -splitD -bed12 -i ${output_dir}/total_transcriptome_primary_sorted.bam | bedtools getfasta -name -split -s -fi $output_dir/genome.fa -bed - -fo ${output_dir}/read_reference_sorted.fa
	rm $output_dir/genome.fa
else
	bedtools bamtobed -splitD -bed12 -i ${output_dir}/total_transcriptome_primary_sorted.bam | bedtools getfasta -name -split -s -fi $index.fa -bed - -fo ${output_dir}/read_reference_sorted.fa
	
fi



##	3) Running of methylation_call.pl script
## the script will compare, for each alignement in bam file, C/T residues in original fastq read with C residues in corresponding genome reference to identify unconverted C
## output is saved to:
##	I) methylated_C_sorted.bed.gz (1 nt bed intervals for each unconverted C in original fastq read)
##	II)total_C_sorted.bed.gz (1 nt intervals for each C in genome that was covered by a read)
##	in both cases read name is retained. These files may be exploited to verify association of methylated C on the same read in future development.

perl $SCRIPT_DIR/perl_scripts/methylation_call_and_coverage.pl ${output_dir}/total_transcriptome_primary_sorted.bam ${output_dir}/trimmed.fastq ${output_dir}/read_reference_sorted.fa ${output_dir}
#rm ${output_dir}/trimmed.fastq ${output_dir}/header.txt ${output_dir}/read_reference_sorted.fa 

## Summarizing coverage of single (methylated) cytosines in bed format 

bedtools merge -i ${output_dir}/methylated_C_sorted.bed.gz  -d -1 -s -c 6,4,5 -o first,count_distinct,sum | awk -v OFS='\t' '{print $1,$2,$3,$(NF-1),$(NF-0),$(NF-2)}' | pigz > ${output_dir}/methylated_C_summary_sorted.bed.gz
bedtools merge -i ${output_dir}/total_C_sorted.bed.gz  -d -1 -s -c 6,4,5 -o first,count_distinct,sum | awk -v OFS='\t' '{print $1,$2,$3,$(NF-1),$(NF-0),$(NF-2)}' | pigz > ${output_dir}/total_C_summary_sorted.bed.gz

## putting together info on coverage and methylation. fields legend:
## 1) chr
## 2) start
## 3) end
## 4) . (intentionally left empty to comply with bed format)
## 5) . (intentionally left empty to comply with bed format)
## 6) strand
## 7) ratio of methylated C (field 8) / coverage (field 9)
## 8) number of reads in which the C was methylated (not converted)
## 9) number of reads overlapping that position (converted + not converted) 

bedtools intersect -sorted -a ${output_dir}/total_C_summary_sorted.bed.gz -b ${output_dir}/methylated_C_summary_sorted.bed.gz -s -wao | LC_ALL=C awk -v OFS='\t' '{if ($13 =="0") print $1,$2,$3,".",".",$6,"0","0",$5; else print $1,$2,$3,".",".",$6,$11/$5,$11,$5}' > ${output_dir}/total_transcriptome_sorted.bed

#### Sequencing coverage data

## sort bam file by genome position

samtools sort  -@ $cores ${output_dir}/total_transcriptome_primary_sorted.bam -o ${output_dir}/total_transcriptome_primary_pos_sorted.bam

## compute strand-specific coverage on + and - strand at 1 nt resolution (not only C) and merge into a single file.

bedtools genomecov -ibam ${output_dir}/total_transcriptome_primary_pos_sorted.bam -dz -strand - -split | awk -v OFS='\t' '{print $1,$2,$2+1,".",$3,"-"}' > ${output_dir}/coverage.bed
bedtools genomecov -ibam ${output_dir}/total_transcriptome_primary_pos_sorted.bam -dz -strand + -split | awk -v OFS='\t' '{print $1,$2,$2+1,".",$3,"+"}' >> ${output_dir}/coverage.bed
LC_ALL=C sort -S 2G -k1,1 -k2,2n ${output_dir}/coverage.bed > ${output_dir}/coverage_sorted.bed

## extract coverage for each microRNA hairpin at 1 nt resolution.

bedtools intersect -wo -s -sorted -a  ${output_dir}/coverage_sorted.bed -b $SCRIPT_DIR/$species/hairpin.bed |  awk -v OFS='\t' '{if ($6 == "+") print $10,$2-$8,$5; if ($6 == "-") print $10,$9-$3,$5}' | sort -k1,1 -k2,2n > ${output_dir}/coverage_hairpin_sorted.bed

### clean-up

pigz -f ${output_dir}/coverage_sorted.bed
#rm ${output_dir}/coverage.bed ${output_dir}/total_transcriptome_primary_sorted.bam

## Final bed output/Star generation/Additional statistics generation

bedtools intersect -wo -s -sorted -a ${output_dir}/total_transcriptome_sorted.bed -b $SCRIPT_DIR/$species/hairpin.bed > ${output_dir}/hairpin_intermediate.bed
bedtools intersect -wo -s -v -sorted -a ${output_dir}/hairpin_intermediate.bed -b $SCRIPT_DIR/$species/mature.bed | awk -v OFS='\t' '{print $1,$2,$3,$13,$5,$6,$7,$8,$9}' - > ${output_dir}/unmature_subbed.bed

bedtools intersect -wo -s -sorted -a ${output_dir}/total_transcriptome_sorted.bed -b $SCRIPT_DIR/$species/mature.bed | awk -v OFS='\t' '{print $1,$2,$3,$13,"*",$6,$7,$8,$9}' - > ${output_dir}/mature_subbed.bed

cat ${output_dir}/mature_subbed.bed ${output_dir}/unmature_subbed.bed | sort -k1,1 -k2,2n - | bedtools intersect -wo -s -sorted -a - -b $SCRIPT_DIR/$species/hairpin.bed | awk -v OFS='\t' '{if ($6 == "+"){print $13,$2-$11,$3-$11,$4,$5,$6,$7,$8,$9} else if ($6 == "-"){print $13,$12-$3,$12-$2,$4,$5,$6,$7,$8,$9}}' - | sort -k1,1 -k2,2n - > ${output_dir}/hairpin_subbed.bed
awk -v OFS='\t' '{print $1,$3,$9,$7}'  ${output_dir}/hairpin_subbed.bed > ${output_dir}/bisRNA_input

Rscript --vanilla $SCRIPT_DIR/R_scripts/bisRNA_analysis.R ${output_dir}
 
awk -v OFS='\t' '{print $1,$2,$3,$7,$8,$9,$5,$4,$6}' ${output_dir}/hairpin_subbed.bed | awk '($6 >= 10) {print}' - | awk '($4 >= 0.05) {print}' > ${output_dir}/5mC_star.txt
cut -f 1 ${output_dir}/5mC_star.txt | uniq | fgrep -w -f - ${output_dir}/hairpin_subbed.bed | awk -v OFS='\t' '{print $1,$2,$3,$7,$8,$9,$5,$4,$6}' - > ${output_dir}/5mC_cov_star.txt

bedtools intersect -wb -f 1 -s -abam ${output_dir}/total_transcriptome_primary_pos_sorted.bam -b $SCRIPT_DIR/$species/hairpin.bed -bed |  awk -v OFS='\t' '{print $16}' - | sort | uniq -c | sed 's/^[ ]*//;s/ /\t/' | sort -n -r > ${output_dir}/hairpin_aligned_statistics.txt
echo -e "\nAligned to miRNAs hairpin:\n" >> $output_dir/logs/alignment.log
cut -f 1 ${output_dir}/hairpin_aligned_statistics.txt | paste -sd + - | bc >> $output_dir/logs/alignment.log
pigz -f ${output_dir}/total_transcriptome_sorted.bed
#rm ${output_dir}/hairpin_intermediate.bed ${output_dir}/${output_dir}/mature_subbed.bed ${output_dir}/unmature_subbed.bed
sed 's/_0\+/\t/' ${output_dir}/bisRNA_output | awk -v OFS='\t' '{print $1,$2-1,$2,".",$4,"."}' | bedtools intersect -b - -a ${output_dir}/hairpin_subbed.bed -wao | awk -v OFS='\t' '{ print $1,$2,$3,$4,$14,$6,$7,$8,$9,$5}' | sed 's/\t-1\t/\t1\t/' > ${output_dir}/hairpin_subbed_sig.bed
prepare_data_for_R () {
	grep -P "$1\t"  $2/coverage_hairpin_sorted.bed >  $2/data_for_R/$1.cov
	grep -P "$1\t"  $2/hairpin_subbed_sig.bed >  $2/data_for_R/$1.bed
	
}

export -f prepare_data_for_R

if [ ! -d ${output_dir}/data_for_R ]; then
	mkdir ${output_dir}/data_for_R
fi

cut -f 1 ${output_dir}/hairpin_subbed_sig.bed | uniq | parallel -k prepare_data_for_R {} ${output_dir}

if [ ! -d ${output_dir}/images ]; then
	mkdir ${output_dir}/images
	mkdir ${output_dir}/images/methylated
	mkdir ${output_dir}/images/not_methylated
fi


cut -f 1 ${output_dir}/hairpin_subbed_sig.bed | uniq | parallel -k Rscript --vanilla $SCRIPT_DIR/R_scripts/plot_mir_meth_log_1_sample.R {} ${output_dir}/  $SCRIPT_DIR/$species $short_name $color $log
cp ${output_dir}/hairpin_subbed_sig.bed ${output_dir}/Summary_report.txt
sed '1,1s/^/pre-miRNA\tstart\tend\tmiRNA\tadj P-value\tstrand (genome)\tNon-converted frequency\tnon-converted reads\ttotal reads\tis in mature miRNA?\n/;s/\.$/No/;s/\*$/Yes/' -i ${output_dir}/Summary_report.txt
exit
