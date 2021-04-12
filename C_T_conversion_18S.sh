#!/bin/bash

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
output_dir=$1
species="hsa"
#sed '2~4s/C/T/g;' $output_dir/trimmed.fastq > $output_dir/trimmed_CT_converted.fastq

mkdir $output_dir/C_T_efficiency_18S
bowtie2 -q -N 0 -L 18 -p 8 --score-min L,0,-0.2 --norc -x $SCRIPT_DIR/$species/18S_C_T_converted -U $output_dir/trimmed_CT_converted.fastq | samtools view -Sb - > $output_dir/C_T_efficiency_18S/18s_C_T.bam

samtools view -H $output_dir/C_T_efficiency_18S/18s_C_T.bam > $output_dir/C_T_efficiency_18S/header.txt

samtools view -F 256 -F 4 $output_dir/C_T_efficiency_18S/18s_C_T.bam | awk '{if (($6 !~ /[I,D]/) &&  ($6 !~ /N.+N/)) print}' - | LC_ALL=C sort -S 2G -k 1,1n -t '@' - | cat $output_dir/C_T_efficiency_18S/header.txt - | samtools view -Sb - > $output_dir/C_T_efficiency_18S/18s_C_T_clean.bam
bedtools bamtobed -i $output_dir/C_T_efficiency_18S/18s_C_T_clean.bam | bedtools getfasta -name -split -s -fi $SCRIPT_DIR/$species/18S.fa -bed - -fo $output_dir/C_T_efficiency_18S/18s_C_T.fa

perl /media/DATA/mamba/perl_scripts/methylation_call_and_coverage.pl $output_dir/C_T_efficiency_18S/18s_C_T_clean.bam $output_dir/trimmed.fastq $output_dir/C_T_efficiency_18S/18s_C_T.fa $output_dir/C_T_efficiency_18S/
mergeBed -i $output_dir/C_T_efficiency_18S/methylated_C_sorted.bed.gz -d -1 -s -c 4,5 -o count_distinct,sum | awk -v OFS='\t' '{print $1,$2,$3,$5,$6,$4 }' | pigz > $output_dir/C_T_efficiency_18S/methylated_C_summary_sorted.bed.gz
mergeBed -i $output_dir/C_T_efficiency_18S/total_C_sorted.bed.gz -d -1 -s -c 4,5 -o count_distinct,sum | awk -v OFS='\t' '{print $1,$2,$3,$5,$6,$4 }' | pigz > $output_dir/C_T_efficiency_18S/total_C_summary_sorted.bed.gz
intersectBed -sorted -a $output_dir/C_T_efficiency_18S/total_C_summary_sorted.bed.gz -b $output_dir/C_T_efficiency_18S/methylated_C_summary_sorted.bed.gz -s -wao | awk -v OFS='\t' '{if ($13 =="0") print $1,$2,$3,".",".",$6,"0","0",$5; else print $1,$2,$3,".",".",$6,$11/$5,$11,$5}' > $output_dir/C_T_efficiency_18S/total_transcriptome_sorted.bed
 cat $output_dir/C_T_efficiency_18S/total_transcriptome_sorted.bed | Rscript --vanilla -e 'aa <- read.table("stdin"); 1-(sum(aa[,8])/sum(aa[,9]));'
cat $output_dir/C_T_efficiency_18S/total_transcriptome_sorted.bed | Rscript --vanilla -e 'aa <- read.table("stdin"); 1-(sum(aa[,8])/sum(aa[,9]));' > $output_dir/C_T_efficiency_18S/efficiency.log
exit

