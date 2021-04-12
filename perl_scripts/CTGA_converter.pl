#!/usr/bin/perl
use strict;
use warnings;

#opening of fasta file

open(REF, "<"."$ARGV[0]") or die "Can't open $ARGV[0]: $!";

#match and substitution
my $strand = "";
while (my $line = <REF>) {
	chomp $line;
	if ($line =~ />/) {
		print $line."\n";
		if ($line =~ /\+/) {$strand= "+";} else {$strand="-";}
	}
	else { if ($strand eq "+" ) {$line =~ s/C/T/g;} else {$line =~ s/G/A/g;}
	print $line."\n";
	}	
}
close REF;

exit;
