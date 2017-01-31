use Cwd 'abs_path';
use File::Basename;
use lib dirname( abs_path $0 );

use strict;
use warnings;
use CGP_lib qw/f_parse_fasta/;

my $dir;
if (abs_path($0)=~/^(.+)\/script_CGP_main\.pl$/) {
	$dir = $1;
}

# command example:
# perl script_CGP_main.pl --segment-length=20 --output-file=tmpoutput --name-file=listOfGenomes.txt test_b*fasta



my @coreGeneFastaFiles;
my ($segmentLength,$outputFileName,$file_genome_names,$file_input_tree) = (20,'output','','');
for my $option (@ARGV) {
	if ($option=~/^--segment-length=(.+)$/) {
		$segmentLength = $1;
	} elsif ($option=~/^--output-file=(.+)$/) {
		$outputFileName = $1;
	} elsif ($option=~/^--name-file=(.+)$/) {
		$file_genome_names = $1;
		if (not -e $file_genome_names) {
			print STDERR "file $file_genome_names does not exists, the genomes will be labelled as 1, 2, ...\n";
			$file_genome_names = '';
		}
	} elsif ($option=~/^--input-tree=(.+)$/) {
		$file_input_tree = $1;
		if (not -e $file_input_tree) {
			print STDERR "newick tree file $file_input_tree does not exists\n";
			$file_input_tree = '';
		}
	} else {
		if (-e $option) {
			push @coreGeneFastaFiles,$option;
		} else {
			print STDERR "file $option does not exists";
		}
	}
}



my $filename_supergene_phylip = $outputFileName.'_trimmedSupergene.phy';

# trim core gene fasta files, and concatenate the core gene sequences of a genome into a 'supergene', and write the supergene of each strain into a phylip file
print STDERR "trim the sequences the save them in a file\n\n";
if ($file_genome_names) {	
	my $core_gene_file_txt = join(' ',@coreGeneFastaFiles);
	system("perl $dir\/script_CGP1_concatenateCoreGenes.pl $segmentLength $filename_supergene_phylip $file_genome_names $core_gene_file_txt");
} else {
	# when genomes names are not given, then simply call them 1, 2, ...
	my ($hash_genename2seq,$arr_genenames) = f_parse_fasta($coreGeneFastaFiles[0]);
	my $Ngenomes = scalar @$arr_genenames;
	my $tmpfile = 'tmp.txt';
	open FOUT,">$tmpfile" or die "can't write file\n";
	for my $q (1..$Ngenomes) {
		print FOUT "$q\n";
	}
	close FOUT;
		
	my $core_gene_file_txt = join(' ',@coreGeneFastaFiles);
	system("perl $dir\/script_CGP1_concatenateCoreGenes.pl $segmentLength $filename_supergene_phylip $tmpfile $core_gene_file_txt");
	
	unlink $tmpfile;
}

# convert the phylip file into SSP distribution file for CGP algorithm
print STDERR "extract the SSP distribution of all the sequence pairs\n\n";
my $filename_strain_pairs = $outputFileName.'_pairName.dat';
my $filename_SSP_each_segment = $outputFileName.'_pairSSP.dat';
system("perl $dir\/script_CGP2_phylip2sspDistribution.pl $segmentLength $filename_supergene_phylip > $filename_strain_pairs 2> $filename_SSP_each_segment");
unlink $filename_supergene_phylip;

# perform the MCMC simlation
print STDERR "perform MCMC\n";
system("perl $dir\/script_CGP3_MCMC.pl $filename_strain_pairs $filename_SSP_each_segment $segmentLength 1000 $outputFileName $file_input_tree 2> /dev/null");
