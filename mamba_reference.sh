#!/bin/bash

## (I) Index generation
## (I.1) rRNA_tRNA index generation


SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

show_help () { cat $SCRIPT_DIR/readme_ref_prep.txt; }



while getopts "h?r:g:i:s:b:" opt; do
    case "$opt" in
    h|\?)
        show_help
        exit 0
        ;;
    r)  rRNA_tRNA_file=$OPTARG
        ;;
    g)  gtf=$OPTARG
        ;;
    i)  index_bt2=$OPTARG
        ;;
    s)  species=$OPTARG
        ;;
    b)  genome_build=$OPTARG
        ;;
    c)  rRNA_18S=$OPTARG
        ;;
    esac
done
shift $((OPTIND-1))



[ "${1:-}" = "--" ] && shift

if [ ! -d $SCRIPT_DIR/$species ]; then
	mkdir $SCRIPT_DIR/$species/ 
fi

if [ ! -d $SCRIPT_DIR/$species/rRNA_tRNA_converted ]; then
mkdir $SCRIPT_DIR/$species/rRNA_tRNA_converted 
fi

mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e "select chrom, size from $genome_build.chromInfo"  > $SCRIPT_DIR/$species/chrom.sizes
sed '/^>/!s/C/T/g;' $rRNA_tRNA_file > $SCRIPT_DIR/$species/rRNA_tRNA_converted/rRNAs_tRNAs_converted.fa
sed '/^>/!s/C/T/g;' $rRNA_18S > $SCRIPT_DIR/$species/18S_C_T_converted
sed '/^>/!s/C/T/g;' $rRNA_28S > $SCRIPT_DIR/$species/28S_C_T_converted


echo -e "\n\n\tBuilding bowtie2 index for fasta file $rRNA_tRNA_file\n\n"
bowtie2-build $SCRIPT_DIR/$species/rRNA_tRNA_converted/rRNAs_tRNAs_converted.fa $SCRIPT_DIR/$species/rRNA_tRNA_converted/rRNAs_tRNAs_converted > $SCRIPT_DIR/$species/rRNA_tRNA_converted/index.log 2>&1 
bowtie2-build $SCRIPT_DIR/$species/18S_C_T_converted $SCRIPT_DIR/$species/18S_C_T_converted > $SCRIPT_DIR/$species/18S_index.log 2>&1 

## (I.2) Total transcriptome index generation

cp $gtf $SCRIPT_DIR/$species/annotation.gtf

if [ ! -d $SCRIPT_DIR/$species/total_transcriptome ]; then
	mkdir $SCRIPT_DIR/$species/total_transcriptome
fi

echo -e "\n\n\tRunning Tophat2 to build bowtie2 index for transcriptome in $gtf\n\n"

tophat2 -G $SCRIPT_DIR/$species/annotation.gtf --transcriptome-index $SCRIPT_DIR/$species/total_transcriptome $index_bt2 > $SCRIPT_DIR/$species/total_transcriptome/index.log 2> $SCRIPT_DIR/$species/total_transcriptome/indexerr.log

if [ ! -d $SCRIPT_DIR/$species/total_transcriptome_converted ]; then
	mkdir $SCRIPT_DIR/$species/total_transcriptome_converted/
fi

perl $SCRIPT_DIR/perl_scripts/CTGA_converter.pl $SCRIPT_DIR/$species/total_transcriptome/annotation.fa > $SCRIPT_DIR/$species/total_transcriptome_converted/annotation.fa


echo -e "\n\n\tBuild bowtie2 index for C --> U converted transcriptome (same transcripts as in $gtf)\n\n"

bowtie2-build $SCRIPT_DIR/$species/total_transcriptome_converted/annotation.fa $SCRIPT_DIR/$species/total_transcriptome_converted/annotation > $SCRIPT_DIR/$species/total_transcriptome_converted/index.log 2>&1 

mv $SCRIPT_DIR/$species/total_transcriptome/annotation.fa.tlst $SCRIPT_DIR/$species/total_transcriptome_converted/
mv $SCRIPT_DIR/$species/total_transcriptome/annotation.gff $SCRIPT_DIR/$species/total_transcriptome_converted/
mv $SCRIPT_DIR/$species/total_transcriptome/annotation.ver $SCRIPT_DIR/$species/total_transcriptome_converted/

rm -r $SCRIPT_DIR/$species/total_transcriptome/
rm $SCRIPT_DIR/$species/annotation.gtf

## (I.3) Sorting of hairpin and mature bed files
wget -O $SCRIPT_DIR/$species/microRNA.gff3 ftp://mirbase.org/pub/mirbase/CURRENT/genomes/$species.gff3
cut -f 1,4,5,7,9 $SCRIPT_DIR/$species/microRNA.gff3 | sed -n /^chr/p |  fgrep  MIMAT | sed 's/ID=.*Name=\(.*\);Derive.*/\1/' | awk -v OFS='\t' '{print $1,$2-1,$3,$5,".",$4}' | sort -k1,1 -k2,2n > $SCRIPT_DIR/$species/mature.bed
cut -f 1,4,5,7,9 $SCRIPT_DIR/$species/microRNA.gff3 | sed -n /^chr/p |  fgrep -v MIMAT | sed 's/ID=.*=//' | awk -v OFS='\t' '{print $1,$2-1,$3,$5,".",$4}' | sort -k1,1 -k2,2n > $SCRIPT_DIR/$species/hairpin.bed
 intersectBed -a $SCRIPT_DIR/$species/mature.bed -b $SCRIPT_DIR/$species/hairpin.bed -s -wo | awk -v OFS='\t' '{if ($6=="+") print $2-$8,$3-$8,$10;else print $9-$3,$9-$2,$10}' > $SCRIPT_DIR/$species/all_positions
if [ ! -d $SCRIPT_DIR/$species/anno_for_R ]; then
	mkdir $SCRIPT_DIR/$species/anno_for_R/
fi

generate_anno_for_R () {
	grep -w $1$ $2/all_positions > $2/anno_for_R/$1.pos
	grep -P "$1\t" $2/hairpin.bed | fastaFromBed -fi $3.fa -bed - -s -tab | sed -n 1,1p  | cut -f 2 | sed s/T/U/g > $2/anno_for_R/$1.fa
}
export -f generate_anno_for_R 
if [ -f $index_bt2.fa ]; then
	echo -e "\n\n genome fasta file exists: $index_bt2.fa; extracting micorRNA hairpins sequences\n\n"
else 
	echo -e "\n\n genome fasta file does not exist: generating $index_bt2 with bowtie2-inspect\n\n"
	bowtie2-inspect $index_bt2 > $index_bt2.fa
	echo -e "extracting micorRNA hairpins sequences\n\n"
fi

cut -f 4 $SCRIPT_DIR/$species/hairpin.bed | sort | uniq | parallel -k generate_anno_for_R {} $SCRIPT_DIR/$species $index_bt2
exit
