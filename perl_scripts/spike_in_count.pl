use warnings;
use strict;

open (IN, "zcat $ARGV[0] | sed -n 2~4p |") || die "cannnot open $ARGV[0]";
#$| = 10;
my $spike_seq = $ARGV[1];
$spike_seq =~ s/U/T/g;
my $out_folder=$ARGV[2];
my $sequence;
my $spike_in_number = 0;
my $count;
my $string; 
my $unconverted_C;
my $counter;
my $millions = 0;


print "\n\n\tSearching for spike-in seq $spike_seq in $ARGV[0] fastq file\n\n";
#print "\tWarning - This program only looks for perfect matches to spike-in \n\t(taking into account C --> U conversion may or may not happen)\n\tThis implies that those (hopefully few) spike-ins with sequencing errors will be missed\n\n";

my $num_C_per_spike = ($spike_seq  =~ s/C/[CT]/g);
#print "$spike_seq\n\n";
while ($sequence = <IN>) {
	chomp $sequence;
	$counter +=1;
	if ($sequence =~ m/^($spike_seq)/) {
		$spike_in_number += 1;
		
		
		$string = $1;
		$count = $string =~ tr/C/C/;
		$unconverted_C += $count;
		}
	$millions =sprintf("%.2f",$counter/1000000);
	print "\r\tFound $spike_in_number spike-in sequences in ${millions} Millions reads analysed thus far";
}
close IN;
#print "\n\n\n";

if ($spike_in_number < 5) {
	print "\tToo few spike-in reads found (please check carefully the control spike-in seq you provided)\n\tThis could be due to \n\t1) wrong spike in sequence provided (restart analysis providing the correct one)\n\t2)Concentration of spike-in RNA oligo is too low\n\t3)Unphsphorylated spike-in RNA oligo: did you use a 5' phsphorylated RNA oligo as spike-in?\n\n\n";
	open (SPIKE_LOG, ">$ARGV[2]/logs/spike_in_analysis.log") || die "cannot open logfile";
	print SPIKE_LOG "Analysis of file $ARGV[0]\nSpike-in sequence: $ARGV[1]\n";
	print SPIKE_LOG "Found only $spike_in_number reads matching the spike in sequence you provided.\n This could be due to \n1) wrong spike in sequence provided (restart analysis providing the correct one)\n2)Concentration of spike-in RNA oligo is too low\n3)Unphsphorylated spike-in RNA oligo: did you use a 5' phsphorylated RNA oligo as spike-in?\n\n"; 
	exit;
	}

my $spike_in_freq = sprintf("%.3f", $spike_in_number/$counter*1000000);
my $total_C_in_spikes = $spike_in_number * $num_C_per_spike;
my $conversion_rate=($total_C_in_spikes-$unconverted_C)/$total_C_in_spikes;
my $percent_conversion_rate = sprintf ("%.1f", $conversion_rate*100);
#print "\tSpike-in accounts for $spike_in_freq Reads per Million (RPM) in your library.\n\tVery low (< 10 RPM) values could be due to \n\t1) wrong spike in sequence provided (restart analysis providing the correct one)\n\t2)Concentration of spike-in RNA oligo is too low\n\t3)Unphsphorylated spike-in RNA oligo: did you use a 5' phsphorylated RNA oligo as spike-in?\n\n\n"; 
print "\n\tFound $spike_in_number spike-in sequences, accounting for a total of $total_C_in_spikes C before bisulfite treatment.\n";

print "\tOf those $total_C_in_spikes C residues, $unconverted_C were not converted into U, with a conversion efficiency of $percent_conversion_rate%.\n\n\n";
open (SPIKE_LOG, ">$ARGV[2]/logs/spike_in_analysis.log") || die "cannot open logfile";
print SPIKE_LOG "Analysis of file $ARGV[0]\nSpike-in sequence: $ARGV[1]\n";
print SPIKE_LOG "Spike-in accounts for $spike_in_freq Reads per Million (RPM) in your library.\nVery low (< 10 RPM) values could be due to \n1) wrong spike in sequence provided (restart analysis providing the correct one)\n2)Concentration of spike-in RNA oligo is too low\n3)Unphsphorylated spike-in RNA oligo: did you use a 5' phsphorylated RNA oligo as spike-in?\n\n\n\n"; 
print SPIKE_LOG "Found $spike_in_number spike-in sequences, accounting for a total of $total_C_in_spikes C before bisulfite treatment.\n";
print SPIKE_LOG "Of those $total_C_in_spikes C residues, $unconverted_C were not converted into U, with a conversion efficiency of $percent_conversion_rate%.\n\n\n";
close SPIKE_LOG;

exit;

