#!/usr/bin/perl -w

use Cwd 'abs_path';
use File::Basename;
use lib dirname( abs_path $0 );

use List::Util qw/min max/;
use strict;
use warnings;
use CGP_lib qw/f_print_fasta f_parse_fasta f_print_phy/;


# input command
# 	perl scrip_CGP_main.pl numberOfSitesPerSegment listOfGenomeID.txt alignmentCoreGene1.fasta alignmentCoreGene2.fasta ...
# the genes in each fasta file must of the same of genomeID
my ($NsitePerSegment,$outputPHYfile,$file_genome_names,@coreGeneFasta) = @ARGV;



# parse the genome names
my @genomes;
open FIN,"<$file_genome_names" or die "can't read $file_genome_names\n";
while (<FIN>) {
	chomp;
	push @genomes,$_ if $_;
}
close FIN;



# process the core gene fasta files
my @superGene_seq = ('') x @genomes;
for my $fasta_file (@coreGeneFasta) {
	my ($hash_genename2seq,$arr_genenames) = f_parse_fasta($fasta_file);
	
	# check and see if the sequences have equal length
	for my $n (1..$#$arr_genenames) {
		unless (length($hash_genename2seq->{$arr_genenames->[$n-1]}) == length($hash_genename2seq->{$arr_genenames->[$n]})) {
			die "genes in file $fasta_file have unequal lengths\n\n\n\n\n\n\n\n\n\n\n";
		}
	}
	
	# find out the positions without '-'
	my $range_genome_with_atgc_positions;
	for my $n (0..$#$arr_genenames) {
		my @cur_range = ([-2,-2]);
		my $cur_seq = $hash_genename2seq->{$arr_genenames->[$n]};
		for my $q (0..length($cur_seq)-1) {
			if (substr($cur_seq,$q,1)!~/\-/) {
				if ($cur_range[$#cur_range]->[1]+1==$q) {
					$cur_range[$#cur_range]->[1]++
				} else {
					push @cur_range,[$q,$q];
				}
			}
		}
		shift @cur_range;
		if ($range_genome_with_atgc_positions) {
			$range_genome_with_atgc_positions = f_update_list($range_genome_with_atgc_positions,\@cur_range);
		} else {
			$range_genome_with_atgc_positions = \@cur_range;
		}
	}
	my $seqs = [];
	for my $n (0..$#$arr_genenames) {
		my $cur_seq = $hash_genename2seq->{$arr_genenames->[$n]};
		for my $q (0..$#$range_genome_with_atgc_positions) {
			my ($start,$end) = ($range_genome_with_atgc_positions->[$q][0],$range_genome_with_atgc_positions->[$q][1]);
			$seqs->[$n] .= substr($cur_seq,$start,$end-$start+1);
		}
	}
	
	next unless length($seqs->[0]) >= $NsitePerSegment;
	
	# divide genes into segments, and save the concatenated segments in @superGene_seq
	for (my $q=0; $q<=length($seqs->[0])-1; $q+=$NsitePerSegment) {
		if ($q+$NsitePerSegment<=length($seqs->[0])) {
			for my $n (0..$#$seqs) {
				my $tmpseq = substr $seqs->[$n],$q,$NsitePerSegment;
				$superGene_seq[$n] .= $tmpseq;		
			}
		}
	}	
}



# print the concatenated genomes in a phylip file
f_print_phy($outputPHYfile,\@genomes,\@superGene_seq);
my $outputFASTAfile;
if ($outputPHYfile=~/^(.*)\.phy/) {
	$outputFASTAfile = $1.".fasta";
} else {
	$outputFASTAfile = $outputPHYfile.".fasta";
}
f_print_fasta($outputFASTAfile,\@genomes,\@superGene_seq);









sub f_update_list {
	my ($arr1,$arr2) = @_;
	
	my ($arrA,$arrB) = ($arr1->[0][0]<$arr2->[0][0]) ? ($arr1,$arr2) : ($arr2,$arr1);
	my @arr;
	my $ii_start_A = 0;
	for my $qB (0..$#$arrB) {
		while (1) {
			if (($arrA->[$ii_start_A][0]<=$arrB->[$qB][0]) && ($ii_start_A<$#$arrA) && ($arrA->[$ii_start_A+1][0]<=$arrB->[$qB][0])) {
				$ii_start_A++;
			} else {
				last;
			}
		}
		my $ii_end_A = $ii_start_A;
		for my $qA ($ii_start_A..$#$arrA) {
			if ($arrA->[$qA][1]>=$arrB->[$qB][1]) {
				$ii_end_A = $qA;
				last;
			}
		}
		for my $qA ($ii_start_A..$ii_end_A) {
			if (($arrA->[$qA][0]<=$arrB->[$qB][0]) && ($arrB->[$qB][0]<=$arrA->[$qA][1]) && ($arrA->[$qA][1]<=$arrB->[$qB][1])) {
				push @arr,[$arrB->[$qB][0],$arrA->[$qA][1]];
			} elsif (($arrB->[$qB][0]<=$arrA->[$qA][0]) && ($arrA->[$qA][0]<=$arrB->[$qB][1]) && ($arrB->[$qB][1]<=$arrA->[$qA][1])) {
				push @arr,[$arrA->[$qA][0],$arrB->[$qB][1]];
			} elsif (($arrA->[$qA][0]<=$arrB->[$qB][0]) && ($arrB->[$qB][1]<=$arrA->[$qA][1])) {
				push @arr,[$arrB->[$qB][0],$arrB->[$qB][1]];
			} elsif (($arrB->[$qB][0]<=$arrA->[$qA][0]) && ($arrA->[$qA][1]<=$arrB->[$qB][1])) {
				push @arr,[$arrA->[$qA][0],$arrA->[$qA][1]];
			}
		}		
	}
	@arr = sort {$a->[0] <=> $b->[0]} @arr;
	return \@arr;
}
