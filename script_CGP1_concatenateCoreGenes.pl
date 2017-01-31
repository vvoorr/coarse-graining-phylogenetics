use Cwd 'abs_path';
use File::Basename;
use lib dirname( abs_path $0 );

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
	
	# push the alignment in the array $seqs
	my $seqs = [];
	for my $genename (@$arr_genenames) {
		my @chars = split //,$hash_genename2seq->{$genename};
		push @$seqs,\@chars;
	}
	next unless scalar(@{$seqs->[0]}) >= $NsitePerSegment;
	
	# remove the positions with '-' in $seqs
	for my $q (reverse 0..$#{$seqs->[0]}) {
		my $hasDash = 0;
		for my $n (0..$#$seqs) {
			if ($seqs->[$n][$q] eq '-') {
				$hasDash = 1;
				last;
			}
		}
		if ($hasDash) {
			for my $n (0..$#$seqs) {
				splice @{$seqs->[$n]},$q,1;
			}
		}
	}
	next unless scalar(@{$seqs->[0]}) >= $NsitePerSegment;
	
	# divide genes into segments, and save the concatenated segments in @superGene_seq
	for (my $q=0; $q<=$#{$seqs->[0]}; $q+=$NsitePerSegment) {
		if ($q+$NsitePerSegment-1<=$#{$seqs->[0]}) {
			for my $n (0..$#$seqs) {
				my $tmpseq = join('',@{$seqs->[$n]}[$q..$q+$NsitePerSegment-1]);
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
