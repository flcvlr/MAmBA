	mamba_compare.sh 

This script will compare 5(h)mC rates in two different samples analyzed with mamba.sh.

Usage:	mamba_compare.sh -g annotation.gtf -r rRNA_tRNA.fa -s species -i bowtie2_index

Options:

-h/-?	print this help message

-t	target sample. A directory name corresponding to the output of a mamba
	analysis you already run.

-c	control sample. A directory name corresponding to the output of a mamba
	analysis you already run.

-f 	the minimum 5(h)mC frequency for testing. Only C residues attaining a 5(h)mC 
	higher than this in at least one of the two samples (-t and -c) will be tested.
	0.2 is usually a reasonable value.

-m	the minimum sequencing depth for testing (raw reads). Only C residues attaining 
	sequencing depth higher than this in at least one of the two samples (-t and -c) 
	will be tested. 30 is usually a reasonable value.

-s	the species annotation to use (hsa|mmu)

-o 	output directory. Figures in .png format for each miRNA will be saved here along with
	other files containing details on statistical testing, coverage and 5(h)mC frequency

-j	number of threads to be used

-d	short name to be used in figures for control sample

-u	short name to be used in figures for target sample
	

All options (apart -h/-?) are required.
