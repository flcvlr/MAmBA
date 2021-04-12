	mamba_reference.sh 

This script will generate a reference for mamba and put the relevant files in the appropriate directory.

Usage:	mamba_reference.sh -g annotation.gtf -r rRNA_tRNA.fa -s species -i bowtie2_index

Options:

-h/-?	print this help message

-g	a gtf file containing annotation of transcripts in gtf format
	the gtf file should be based on the same genome version as the
	bowtie2_index provided with the -i option and should include 
	pre-miRNA transcripts

-r	a fasta file containing the sequences of rRNA and tRNAs of the
	organism of interest. rRNA and tRNA fasta files for your species
	of interest can be downloaded from RNAcentral or directly from 
	other databases

-s	the species you want to build a reference for. You will have to 
	provide this exact name through the -s option when running mamba.
	This will be used to retrieve the miRNA annotation files from miRBase
	Therefore you must stick to the common 3 letters code (hsa = Homo sapiens,
	mmu = Mus musculus, dme = Drosophila elanogaster etc) and meake sure by 
	visiting miRBase website you are providing the correct code for your 
	favorite species according to miRBase.

-i	the basename of a valid bowtie2 index for the genome you want to use.
	Please double check that the gtf file provided with the -g option and
	the bowtie2 index are based on the same genome version (e.g hg38, mm10)

-b	The genome build of your bowtie index (e.g hg38, mm10). This will be used to retrieve 
	chromosome length information from UCSC genome browser, so you should make sure
	that the genome build you provide corresponds to a valid genome build 
	according to UCSC genome browser AND that it also corresponds to the 
	bowtie2 index provided with -i option and the gtf file profided with -g option

All options (apart -h/-?) are required.
