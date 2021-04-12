#!/usr/bin/perl
use strict;
use warnings;


#usage perl methylation_call.pl <BAM> <READS> <REFERENCE>

#opening of output sorted bed file
open(BEDMET, "| LC_ALL=C sort -k1,1 -k2,2n -S 40%  | pigz -c > $ARGV[3]/methylated_C_sorted.bed.gz") or die "cannot open output BED: $!";

#opening of output total coverage bed file
open(COV_C, "| LC_ALL=C  sort -k1,1 -k2,2n -S 40%  | pigz -c > $ARGV[3]/total_C_sorted.bed.gz") or die "cannot open output COV: $!";

#opening of alignment sam file
open(SAM, "samtools view -h $ARGV[0] |") or die "Can't open $ARGV[0]: $!";

#opening of fastq reads file

open(READS, "<"."$ARGV[1]") or die "Can't open $ARGV[1]: $!";

#opening of fasta reference file

open(REF, "<"."$ARGV[2]") or die "Can't open $ARGV[2]: $!";

my @n;
my @m;
my $strand;
my $start;
my $end;
my $i;
my @readline = ();
my @reference = ();
my $header_warn = 0;
my @sam;
my @refseq;
my @readseq;

while (my $samline = <SAM>) {
	chomp $samline;
	if ($samline =~ /^\@/) {
		if ($header_warn == 0) {warn "skipping SAM header line\n"; $header_warn=1;}
	next;
	}
	@sam = split(/\t/, $samline);
	while ($readline[0] = <READS>) {
		chomp $readline[0];
		my @cleanreadline = split(/ /, $readline[0]);
		if ($cleanreadline[0] eq "@".$sam[0]) {
			last;
		}
	}

	$readline[1] = <READS>;
	chomp $readline[1];
	$readline[2] = <READS>;
	chomp $readline[2];
	$readline[3] = <READS>;
	chomp $readline[3];

	$reference[0]= <REF>;
	chomp $reference[0];
	$reference[1]= <REF>;
	chomp $reference[1];

	@refseq=split(//,$reference[1]);
	@readseq=split(//,$readline[1]);


###### print bedgraph_output ########

if ($sam[5] =~ /N/) {
	@m = split(/\d+N/, $sam[5]);
	@n = split(/\d+M/, $sam[5]);
	chop $m[0];
	chop $m[1];
	chop $n[1];
	for ($i=0; $i < length($reference[1]); $i=$i + 1) {
		if ($refseq[$i] eq 'C') {
				if ($sam[1] == 0)	{
				if($i < $m[0]) {$start=$sam[3]-1+$i;}
				if($i >= $m[0]) {$start=$sam[3]-1+$i+$n[1];}
				$strand="+";
				}
			else {
				if($i < $m[1]) {$start=$sam[3]+length($reference[1])-$i-2+$n[1];}
				if($i >= $m[1]) {$start=$sam[3]+length($reference[1])-$i-2;}
				$strand="-";
				}
			$end = $start+1;
			print COV_C "$sam[2]\t$start\t$end\t$sam[0]\t1\t$strand\n";
			if ($readseq[$i] eq 'C') {
			print BEDMET "$sam[2]\t$start\t$end\t$sam[0]\t1\t$strand\n";
			}
		}	

	}

	next;
}

for ($i=0; $i < length($reference[1]); $i=$i + 1) {
		if ($refseq[$i] eq 'C') {
				if ($sam[1] == 0)	{
				$start=$sam[3]-1+$i;
				$strand="+";
				}
			else {
				$start=$sam[3]+length($reference[1])-$i-2;
				$strand="-";
				}
			$end = $start+1;
			print COV_C "$sam[2]\t$start\t$end\t$sam[0]\t1\t$strand\n";
			if ($readseq[$i] eq 'C') {
			print BEDMET "$sam[2]\t$start\t$end\t$sam[0]\t1\t$strand\n";
			}
		}
	}
}


close SAM;
close READS;
close REF;
close BEDMET;
close COV_C;

exit;
