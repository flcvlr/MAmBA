

Usage: 

	mamba.sh [options] <-s species> <-f inputfile.fastq.gz> <-i bowtie2_index_basename> 


Mandatory Arguments:

  -f fastq.gz input file
	a gzip compressed fastq file containing reads from a small RNA-seq experiment
	reads should be "stranded" (second strand), i.e. the sequence of the read should
	be the same as the original small RNA. This is the usual configuration of 
	current small RNA-seq protocols

  -i bowtie2 index 
	basename of a bowtie2 index for the genome of interest. If your bowtie2 index is
	/home/user/data/hg38.1.bt2 you should use /home/user/data/hg38
	MAmBA assumes that the corresponding fasta is in the same directory 
	with the same basename. If missing, a corresponding fasta file will be created by MAmBA
	in the output directory and erased at the end of the analysis. If you wish to create
	a fasta file corresponding to the bowtie2 index in the same folder (so that MAmBA
	will not spend time to create it) you can run:
	
	bowtie2-inspect /path/to/your/bowtie2index/basename > /path/to/your/bowtie2index/basename.fa

  -s species (hsa|mmu)
	Currently there are only references for human (hg38) and mouse (mm10) which are 
	distributed here: https://sites.google.com/a/uniroma1.it/valeriofulci-eng/software 
	However, MAmBA includes a script (mamba_reference.sh) to generate custom annotation 
	provided that you have annotation files for your favorite species. 
	For details please run:
	
	mamba_reference.sh -h 
	
  -n short_name to use in images
	MAmBA will output .png images for all pre-miRNAs with sufficient coverage in the 
	output directory. You should provide a human readable name to be used in figures 
	(e.g patient_1 or HeLa_cells). Please avoid spaces in the short name, use underscore.

Command line options:

  -o output_directory (default: mamba_out.$infile)
	MAmBA will create this directory (default: mamba_out.fastq.gz) or, in case it exists 
	already ask for confirmation before overwriting previous output.

  -j number_of_threads (default: 1)
	The number of threads to be used by MAmBA 

  -a adaptor sequence 
	The sequence of the 3' adapter used in NGS. should be at least 10 nucleotides long 
	only A,C,G,T are allowed. By default is set to "no", will cause adaptor trimming skip, 
	assuming you have already trimmed adapters from your fastq file

  -c control spike in
	The sequence of a control spike-in 5' phosphorylated RNA oligo that you have added to the 
	RNA before bisulfite treatment. If provided, C to U conversion will be estimated 
	(assuming that your RNA oligo only contains unmethylated C). Alternatively, the 
	sequence of an endogenous non-methylated small RNA can be provided 
	(e.g. tRNA(Gly) 5' fragment)

  -l log/linear (default: log)
	If set to any value other than "log" a linear y axis (non-conversion frequency) 
	will be used in figures. 
	
  -C color (default: darkgreen)
	If set the argument will be passed to R as the color to use for bars in figures
	(useful to plot multiple samples in different colors)

  -h|? show this help 


