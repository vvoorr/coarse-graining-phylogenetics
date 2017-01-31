package CGP_lib;
use Storable qw(nstore store_fd nstore_fd freeze thaw dclone);
use threads;
use strict;
use warnings;

use Exporter qw(import);
our @EXPORT_OK = qw/f_print_fasta f_parse_fasta f_print_phy f_read_phy/;




sub f_print_fasta {
	my ($filename,$arr_genenames,$arr_geneseq) = @_;
	die "unequal arry size $#$arr_genenames != $#$arr_geneseq" if $#$arr_genenames != $#$arr_geneseq;
	open FOUT,">$filename" or die "can't write file";
	for my $q (0..$#$arr_genenames) {	
		my @chunks = unpack("(A60)*",$arr_geneseq->[$q]);
		print FOUT join("\n",">".$arr_genenames->[$q],@chunks)."\n";	
	}
	close FOUT;
}



sub f_parse_fasta {
	my ($filename) = @_;
	local $/ = undef;
	
	my $fin;
	if ($filename=~/\.xz$/) {
		open $fin,"xzcat $filename |" or die "can't read $filename\n";
	} else {	
		open $fin,"<$filename" or die "can't read $filename";
	}
	my $content = <$fin>;	
	close $fin;
	chomp $content;
	
	my @geneNames;
	my %gene2content;	
	my @aagene = split /\>/,$content;
	shift @aagene unless $aagene[0];
	foreach my $gene (@aagene) {
		my @bbgene = split /\n/,$gene;
		my $genename = shift @bbgene;
		my $content = join("",@bbgene);
		$gene2content{$genename} = $content;
		push @geneNames,$genename;
	}	
	return (\%gene2content,\@geneNames);
}



sub f_print_phy {
	my ($filename,$arr_ncbis,$arr_geneseq) = @_;
	die "unequal arry size $#$arr_ncbis != $#$arr_geneseq" if $#$arr_ncbis != $#$arr_geneseq;
	my $segLength = length($arr_geneseq->[0]);
	my $Nstrain = scalar @$arr_ncbis;
	open FOUT,">$filename" or die "can't write file";
	print FOUT "$Nstrain $segLength\n";
	for my $q (0..$#$arr_ncbis) {	
		my @chunks = unpack("(A10)*",$arr_geneseq->[$q]);
		my $space = ' ' x (10 - length($arr_ncbis->[$q]));
		print FOUT "$arr_ncbis->[$q]$space".join(' ',@chunks)."\n";
	}
	close FOUT;
}



sub f_read_phy {
	my ($filename) = @_;
	my $fin;
	if ($filename=~/\.xz/) {
		open $fin,"xzcat $filename |" or die "can't read $filename\n";
	} else {
		open $fin,"<$filename" or die "can't read $filename\n";
	}
	my $firstLine = <$fin>;
	chomp $firstLine;
	my ($Nchromosome,$length) = split /\s+/,$firstLine;
	my $ncbi2seq = {};
	while (<$fin>) {
		chomp;
		my ($ncbi,@seqs) = split /\s+/,$_;
		my $seq = join('',@seqs);
		$ncbi2seq->{$ncbi} = $seq;
	}
	close $fin;
	return [$Nchromosome,$length,$ncbi2seq];
}
