#!/usr/bin/perl -w

use Cwd 'abs_path';
use File::Basename;
use lib dirname( abs_path $0 );

use strict;
use warnings;
use CGP_lib qw/f_read_phy/;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use threads;


# command
# perl scrip_CGP2_phylip2sspDistribution.pl $segmentLength $inputPHYfile
my ($segmentLength,$inputPHYfile) = @ARGV;
my $Nthreads = 10; # maximum no. of threads for a multi-thread script



# read the phylip file
my ($Nstrains,$seqLength,$strain2seq) = @{f_read_phy($inputPHYfile)};
my @strainName = keys %$strain2seq;
my %strain2seqSegmentized;
while (my ($strain,$seq) = each %$strain2seq) {
	my @seq = unpack("(A$segmentLength)*",$seq);
	$strain2seqSegmentized{$strain} = \@seq;
}



# print out the strain pairs in STDOUT
# print out the SSP distribution in STDERR
for my $q (0..$#strainName-1) {
	my $strainQ = $strainName[$q];
	for (my $r0 = $q+1; $r0 < @strainName; $r0 += $Nthreads) {
		my @thr;
		for my $r ($r0..min($#strainName,$r0+$Nthreads-1)) {
			my $strainR = $strainName[$r];
			$thr[$r] = threads->new(sub{				
				my @NSNP = (0) x scalar(@{$strain2seqSegmentized{$strainQ}});
				for my $s (0..$#NSNP) {
					for my $t (0..length($strain2seqSegmentized{$strainQ}->[$s])-1) {
						if (uc(substr($strain2seqSegmentized{$strainQ}->[$s],$t,1)) ne uc(substr($strain2seqSegmentized{$strainR}->[$s],$t,1))) {
							$NSNP[$s]++;
						}
					}
				}				
				return [@NSNP];
			});
		}
		for my $r ($r0..min($#strainName,$r0+$Nthreads-1)) {
			my $strainR = $strainName[$r];
			my $tmpout = $thr[$r]->join();
			my (@NSNP) = @$tmpout;
			print "$strainQ\t$strainR\n";
			print STDERR join("\t",@NSNP)."\n";
		}
	}
}
