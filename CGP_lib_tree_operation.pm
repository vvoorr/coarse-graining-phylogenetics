package CGP_lib_tree_operation;
use strict;
use warnings;
use Storable qw(nstore store_fd nstore_fd freeze thaw dclone);
use Scalar::Util qw(looks_like_number);
use List::Util qw(shuffle max min sum);
use POSIX;

use Exporter qw(import);
our @EXPORT_OK = qw/f_read_newick_file f_cut_and_graft f_read_file_pairName_pairSNP f_extract_clade_distMatrix_hash/;








sub f_combine_zero_height_clade {
	my ($input_o2y,$input_y2o) = @_;
	my $o2y = dclone $input_o2y;
	my $y2o = dclone $input_y2o;
	while (1) {
		my $hasChange = 0;
		foreach my $onode (keys %$o2y) {
#while (my ($y,$height) = each %{$o2y->{$onode}}) {
	#print STDERR "\t$onode\t$y\t$height\n";
	#if (exists $o2y->{$y}) {
		#print STDERR "\texists $y\n"
	#} else {
		#print STDERR "\tnot exists $y\n"
	#}
#}
			my @dsClades = grep {$o2y->{$onode}{$_}==0} grep {looks_like_number($o2y->{$onode}{$_})} grep {exists $o2y->{$onode}{$_}} grep {exists($o2y->{$_})} keys %{$o2y->{$onode}};
			if (@dsClades>0) {
				foreach my $dsClade (@dsClades) {
					my @dsdsClades = keys %{$o2y->{$dsClade}};
					foreach (@dsdsClades) {
						$o2y->{$onode}{$_} = $o2y->{$dsClade}{$_};
						$y2o->{$_} = $onode;
					}
					delete $o2y->{$onode}{$dsClade};
					delete $o2y->{$dsClade};
					delete $y2o->{$dsClade};
				}
				$hasChange = 1;
				last;
			}
		}
		last if not $hasChange;
	}
	my $node_to_tip = {};
	my $ff;
	$ff = sub {
		my ($cur_node) = @_;
		if (exists $node_to_tip->{$cur_node}) {
			return [keys %{$node_to_tip->{$cur_node}}];
		} elsif (exists($y2o->{$cur_node}) && not(exists $o2y->{$cur_node})) {
			$node_to_tip->{$cur_node}{$cur_node} = 1;
#print STDERR "tip\t$cur_node\n";
			return [$cur_node];
		} else {
			my @tipnodes;
			my @dsnode = keys %{$o2y->{$cur_node}};
			foreach (@dsnode) {
				my @tmmp = @{$ff->($_)};
				push @tipnodes,@tmmp;
			}
			foreach (@tipnodes) {
				$node_to_tip->{$cur_node}{$_} = 1;
			}
			return \@tipnodes;
		}
	};
	my @root = grep {not exists $y2o->{$_}} keys %$o2y;
	die join(' ',@root)."werwererwr" if @root!=1;
	my $ttttmp = $ff->($root[0]);
	my $node_to_name = {};
	foreach my $node (keys %$node_to_tip) {
		my $name = join('+',sort {$a cmp $b} keys %{$node_to_tip->{$node}});
#print STDERR "$node\:$name\n";
		$node_to_name->{$node} = $name;
	}
	my $out_o2y = {};
	my $out_y2o = {};
	while (my ($o,$v) = each %$o2y) {
#print STDERR "\t\t$o\n";
		my $oname = $node_to_name->{$o};
		for my $y (keys %$v) {
#print STDERR "\t\t\t$y\n";
			my $yname = $node_to_name->{$y};
			$out_o2y->{$oname}{$yname} = $o2y->{$o}{$y};
			$out_y2o->{$yname} = $oname;
		}
	}
	return ($out_o2y,$out_y2o);
}



sub get_tip_node {
	my ($tmpo2y,$node) = @_;
#print "$tmpo2y $node\n";
#print $tmpo2y->{$node}."\n";
	my @tipnodes;
	if (not exists $tmpo2y->{$node}) {
		return [$node];
	} else {
		for my $yNode (keys %{$tmpo2y->{$node}}) {
#print "$node\n";
			push @tipnodes,@{get_tip_node($tmpo2y,$yNode)};
		}
		return \@tipnodes;
	}
}



sub f_combine_single_branch_internal_node {
	my ($input_o2y,$input_y2o) = @_;
	my $o2y = dclone $input_o2y;
	my $y2o = dclone $input_y2o;
	while (1) {
		my $hasCollapsed = 0;
		($o2y,$y2o,$hasCollapsed) = f_combine_single_branch_internal_node_A($o2y,$y2o);
		last unless $hasCollapsed;
	}
	return ($o2y,$y2o);
}
sub f_combine_single_branch_internal_node_A {
	my ($input_o2y,$input_y2o) = @_;
	my $o2y = dclone $input_o2y;
	my $y2o = dclone $input_y2o;
	#my $hasCollapsed = 1;
	while (my ($o,$v) = each %$o2y) {
		if (scalar(keys %$v)==1) {
			my @ys = keys %$v;
			my $y = $ys[0];
			if (exists $y2o->{$o}) {
				my $oo = $y2o->{$o};				
				my $height = $o2y->{$o}{$y} + $o2y->{$oo}{$o};
#print STDERR "\tcollapse\t$y\t$o\t$oo\t$o2y->{$o}{$y} + $o2y->{$oo}{$o} = $height\n";
				delete $y2o->{$y};
				delete $o2y->{$o};
				delete $y2o->{$o};
				delete $o2y->{$oo}{$o};
				$y2o->{$y} = $oo;
				$o2y->{$oo}{$y} = $height;
				#$hasCollapsed = 1;
				return ($o2y,$y2o,1);
			} else {
				delete $y2o->{$y};
				delete $o2y->{$o};
				return ($o2y,$y2o,1);
			}
		}
	}
	return ($o2y,$y2o,0);
}


sub get_newick_lowest_level3 {
	my ($nw) = @_;
	my $toBeReplaced;
	my $bracketArea;
	my $cladeName;
	my $branchLength;
	my $annotation;
#print STDERR "$nw\n";
	my @bracketArea;	# $bracketArea[$q] = [nodename,branch_length,annotations]
	if ($nw!~/\:/) {
		if ($nw =~ m/(\([^\(\)]+\))/) {
			$bracketArea = $1;
			$toBeReplaced = $bracketArea;			
		} else {
			$cladeName = $1;
			#die "wererer $nw\n";
		}
		if ($bracketArea) {
			chop $bracketArea if $bracketArea=~/\)$/;
			$bracketArea = substr $bracketArea,1 if $bracketArea=~/^\(/;

			my @a = split /\,/,$bracketArea;
			for my $tmp (@a) {
				push @bracketArea,[$tmp,'',''];
			}
		}
	} else {
		if ($nw =~ m/(\([^\(\)]+\))([^\:]*)\:([^\[\]\,]*?)(\[.+?\])?([^\[\]\,\)]*)[\,\)]/) {
			$bracketArea = $1;
			$cladeName = $2;			
			$annotation = $4;
			$branchLength = $5;
			$toBeReplaced = $1.$2;
#print STDERR "aa\n";
		} elsif ($nw =~ m/(\([^\(\)]+\))([^\:]*)\:([^\[\]]*?)(\[.+?\])?([^\[\]\,\)]*)$/) {
			$bracketArea = $1;
			$cladeName = $2;			
			$annotation = $4;
			$branchLength = $5;
			$toBeReplaced = $1.$2;
#print STDERR "bb\n";
		} elsif ($nw =~ m/(\([^\(\)]+\))([^\:]*)\:([^\,]+?)[\,\)]/) {
			$bracketArea = $1;
			$cladeName = $2;
			$branchLength = $3;
			$toBeReplaced = $bracketArea.$cladeName;
#print STDERR "cc\n";
		} elsif ($nw =~ m/(\([^\(\)]+\))([^\:]*)\:([^\,]+?)$/) {
			$bracketArea = $1;
			$cladeName = $2;
			$branchLength = $3;
			$toBeReplaced = $bracketArea.$cladeName;
#print STDERR "dd\n";
		} elsif ($nw =~ m/(\([^\(\)]+\))([^\(\)\:]*)$/) {
			$bracketArea = $1;
			$cladeName = $2;
			$toBeReplaced = $bracketArea.$cladeName;
#print STDERR "ee\n";
		} elsif ($nw =~ m/^([^\(\)]+)(.*?)\:([^\[\]]+?)(\[.+?\])$/) {
			$cladeName = $1;
			$branchLength = $2;
			$annotation = $3;
#print STDERR "ff\n";
		} elsif ($nw =~ m/^([^\(\)]+)\:([^\,]*?)$/) {
			$cladeName = $1;
			$branchLength = $2;
#print STDERR "gg\n";
		} elsif ($nw =~ m/^([^\(\)\:]+)$/) {
			$cladeName = $1;
		} else {
			die "eerereor $nw\n";
		}
		
		if ($bracketArea) {
			chop $bracketArea if $bracketArea=~/\)$/;
			$bracketArea = substr $bracketArea,1 if $bracketArea=~/^\(/;

			if ($bracketArea=~/\[/) {
				$bracketArea=~s/\[[^\[\]\:]+?\]//g;
				my @a = split /\,/,$bracketArea;
				for my $tmp (@a) {
					if ($tmp=~/^(.+?)\:([^\[\]]+?)$/) {
						push @bracketArea,[$1,$2,''];
#print STDERR "\tuu\t$1\t$2\n";
					} elsif ($tmp=~/^([^\:]+?)$/) {
						push @bracketArea,[$1,'',''];
					} else {
						die "erewsdfgdrerer\t$tmp\n";
					}
				}
				#$bracketArea=~s/\]\,/\]\/\/\//g;
				#my @a = split /\/\/\//,$bracketArea;
				#my @a = split /\,/,$bracketArea;
				#for my $tmp (@a) {
					#if ($tmp=~/^([^\[\]\:]+?)\:([^\[\]]+?)(\[.+\])$/) {
						#push @bracketArea,[$1,$2,$3];
					#} elsif ($tmp=~/^([^\[\]\:]+?)\:(\[[^\[\]]+?\])([^\[\]]+?)$/) {
						#push @bracketArea,[$1,$3,$2];
#print STDERR "$1\n";
					#} else {
						#die "erewre3243rer\t$tmp\n";
					#}
				#}
			} else {
				my @a = split /\,/,$bracketArea;
				for my $tmp (@a) {
					if ($tmp=~/^(.+?)\:([^\[\]]+?)$/) {
						push @bracketArea,[$1,$2,''];
#print STDERR "\tvv\t$1\t$2\n";
					} elsif ($tmp=~/^([^\:]+?)$/) {
						push @bracketArea,[$1,'',''];
					} else {
						die "erewrerer\t$tmp\n";
					}
				}
			}
		}
	}
#print STDERR "$toBeReplaced\n";
#print STDERR "\tbracketArea:$bracketArea\n";
#print STDERR "\tcladeName:$cladeName\n";
#print STDERR "\tannotation:$annotation\n";
#print STDERR "\tbranchLength:$branchLength\n\n\n";
	return [$toBeReplaced,$bracketArea,$cladeName,$branchLength,$annotation,\@bracketArea];
}

sub f_node2downstreamNode_distance {
	my ($y2o,$o2y) = @_;
	my @root = grep {not exists $y2o->{$_}} keys %$o2y;
	if (scalar(@root)!=1) {
		die "number of root not equal to 1!!";
	}
	my $root = $root[0];
	my $o2allY_distance = {};
	my ($dummy1,$dummy2) = @{f_node2downstreamNode_distance_tmp($o2y,$o2allY_distance,$root)};
	return $o2allY_distance;
}
sub f_node2downstreamNode_distance_tmp {
	my ($o2y,$o2allY_distance,$curNode) = @_;
	my @tips;
	my @dw_nodes;
	if (exists $o2y->{$curNode}) {
		my @ys = keys %{$o2y->{$curNode}};
		for my $y (@ys) {
			$o2allY_distance->{$curNode}{$y} = $o2y->{$curNode}{$y};
			my ($tmp_tips,$tmp_downwardNodes) = @{f_node2downstreamNode_distance_tmp($o2y,$o2allY_distance,$y)};
			push @tips,@$tmp_tips;
			for my $dn (@$tmp_tips,@$tmp_downwardNodes) {
				$o2allY_distance->{$curNode}{$dn} = $o2allY_distance->{$y}{$dn} + $o2y->{$curNode}{$y};
			}
			push @dw_nodes,($y,@$tmp_downwardNodes);
		}
		return [\@tips,\@dw_nodes];
	} else {
		$o2allY_distance->{$curNode}{$curNode} = 0;
		return [[$curNode],[]];
	}
}



sub f_extract_clade_distMatrix_hash {
	my ($curnewick) = @_;
	my ($y2o,$o2y,$all_clade,$rootnode) = @{f_extract_clade($curnewick)};
	my $node2downstreamNode_distance = f_node2downstreamNode_distance($y2o,$o2y);
	my @tipnodes = grep {not exists $o2y->{$_}} keys %$y2o;
	my %tipnodes;
	$tipnodes{$_}=1 foreach @tipnodes;
	my $pairwise_distance = {};
	foreach my $o (grep {not exists $tipnodes{$_}} keys %$o2y) {
		my @ys = keys %{$o2y->{$o}};
		for my $q (0..$#ys-1) {
			my @qtips = grep {exists $node2downstreamNode_distance->{$ys[$q]}{$_}} @tipnodes;
			for my $r ($q+1..$#ys) {
				my @rtips = grep {exists $node2downstreamNode_distance->{$ys[$r]}{$_}} @tipnodes;
				for my $qtip (@qtips) {
					for my $rtip (@rtips) {
						my $tmpdist = ($node2downstreamNode_distance->{$o}{$qtip} + $node2downstreamNode_distance->{$o}{$rtip})/2;
						$pairwise_distance->{$qtip}{$rtip} = $tmpdist;
						$pairwise_distance->{$rtip}{$qtip} = $tmpdist;
#print STDERR "$o\t$qtip\t$rtip\t$tmpdist\n";
					}
				}
			}
		}
	}
	return $pairwise_distance;	
}


sub f_extract_clade {
	my ($curnewick,$hash_node_tipid2name) = @_;
	$curnewick =~ s/\;//g;
	my $id = -1;
	#my %node2tipgroup;
	
	# remove bracket area []
	while (1) {
		my $hasChanged = 0;
		my @txts = ($curnewick=~/(\[[^\]]+\])/g);
		for my $txt (@txts) {
			if ($txt=~/\=/) {
				$hasChanged=1;
				#my $toreplace = qr/$txt/;
#print "\t\t$txt\n";
				$curnewick=~s/\Q$txt\E//g;
#print "\t\t$curnewick\n";
			}
		}		
		last unless $hasChanged;
	}
#print STDERR "\t$curnewick\n";

	my %name2originalName;
	my $tmpy2o = {};
	my $tmpo2y = {};
	my $tmpy2height = {};
#print STDERR "aa$curnewick\n";
	while (1) {
		$id++;
		my ($toBeReplaced,$bracketArea,$tmpcladeName,$branchLength,$annotation,$arrsubnode) = @{get_newick_lowest_level3($curnewick)};
		if (not defined($bracketArea)) {
			last;
		}
#print STDERR "branchLength\t$branchLength\n" if $branchLength;
		my $nodeName;
		if ($tmpcladeName) {
			my $cladeName = $tmpcladeName;
			$cladeName =~ s/\[//g;
			$cladeName =~ s/\]//g;
			if ($cladeName=~/^\d+$/) {
				# this is not cladeName
				$nodeName = "node$id";
#print STDERR "zz1\n"
			} elsif ($cladeName=~/\=/) {
				# this is not cladeName
				$nodeName = "node$id";
#print STDERR "zz2\n"
			} elsif ($cladeName=~/\#/) {
				# this is not cladeName
				$nodeName = "node$id";
#print STDERR "zz2\n"
			} else {
				$nodeName = $cladeName;
#print STDERR "\tzz3cladeName $cladeName\n"
			}
		} else {
			$nodeName = "node$id";
#print STDERR "zz4\n"
		}
#print STDERR "nodeName\t$nodeName\n";
#print STDERR "$curnewick\n";
#print STDERR "\t$toBeReplaced\n" if $toBeReplaced;
#print STDERR "\t$nodeName\n";
		$toBeReplaced = quotemeta $toBeReplaced;
		$nodeName = quotemeta $nodeName;
		$curnewick =~ s/$toBeReplaced/$nodeName/;
		if (scalar(@$arrsubnode)>0) {
			for my $tmp (@$arrsubnode) {
				my ($tmp_ynode_name,$height,@aa) = @$tmp;
#print STDERR "$nodeName\t$tmp_ynode_name\t$height\n";
				my $ynode;
				if ($hash_node_tipid2name) {
					if (exists $tmpo2y->{$tmp_ynode_name}) {
						$ynode = $tmp_ynode_name;
					} elsif (not(exists $tmpo2y->{$tmp_ynode_name}) && exists($hash_node_tipid2name->{$tmp_ynode_name})) {
						#$ynode = $tmp_ynode_name;
						$ynode = $hash_node_tipid2name->{$tmp_ynode_name};
					} else {
#print STDERR "\tskipped\n";
						#next;
						$ynode = $tmp_ynode_name;
					}
				} else {
					$ynode = $tmp_ynode_name;
				}
				$tmpy2o->{$ynode} = $nodeName;
				$tmpo2y->{$nodeName}{$ynode} = $height;
				$tmpy2height->{$ynode} = $height;
#print STDERR "$nodeName\t$ynode\t$height\n";
			}			
		}
	}
#print STDERR "aaa\n\n\n";

	# remove tip nodes that are unnecessary
	($tmpo2y,$tmpy2o) = f_combine_single_branch_internal_node($tmpo2y,$tmpy2o);
	if ($hash_node_tipid2name) {
		my %tmp_tipname;
		for (values %$hash_node_tipid2name) {
			$tmp_tipname{$_} = 1;
		}
		my @unnecessarytip = grep {not exists $tmp_tipname{$_}} grep {not exists $tmpo2y->{$_}} keys %$tmpy2o;
		if (scalar(@unnecessarytip)>0) {
			for my $y (@unnecessarytip) {
				my $o = $tmpy2o->{$y};
#print STDERR "\t$o\t$y\n";
				delete $tmpo2y->{$o}{$y};
				delete $tmpy2o->{$y};
				delete $tmpy2height->{$y};
				($tmpo2y,$tmpy2o) = f_combine_single_branch_internal_node($tmpo2y,$tmpy2o);
			}
		}
	}
#print STDERR "bbb\n\n\n";
#print STDERR "\n";
	my %originalCladeName2name;
	my %name2originalCladeName;
	for my $nnode (keys %$tmpo2y,keys %$tmpy2o) {
#print STDERR "$nnode\n";
		next if exists $originalCladeName2name{$nnode};
		my @a = @{get_tip_node($tmpo2y,$nnode)};
		my %tmp;
		for (@a) {
			$tmp{$_} = 1;
		}
		my $name = join ('+',sort {$a cmp $b} keys %tmp);
		$originalCladeName2name{$nnode} = $name;
		$name2originalCladeName{$name} = $nnode;
#print "$nnode\t$originalCladeName2name{$nnode}\n";
	}
#print STDERR "ccc\n\n\n";
	my $y2o = {};
	my $o2y = {};
	while (my ($k,$v) = each %$tmpy2o) {
		my $ynode = $originalCladeName2name{$k};
		my $onode = $originalCladeName2name{$v};
		$y2o->{$ynode} = $onode;
	}
	while (my ($k,$v) = each %$tmpo2y) {		
		my $onode = $originalCladeName2name{$k};
		foreach my $y (keys %$v) {
			my $ynode = $originalCladeName2name{$y};
			$o2y->{$onode}{$ynode} = $tmpo2y->{$k}{$y}; #$tmpy2height->{$y};
#print STDERR "$onode\t$ynode\t$o2y->{$onode}{$ynode}\n";
		}
	}
	
	# remove clades with 0 length
#foreach my $internal_node (keys %$o2y) {
#print STDERR "int node 0 $internal_node\n";
#}
#exit 1;
	#my ($o2y_original,$y2o_original) = f_combine_zero_height_clade($tmpo2y,$tmpy2o);
#print STDERR "ddd\n\n\n";
	($o2y,$y2o) = f_combine_zero_height_clade($o2y,$y2o);
#print STDERR "eee\n\n\n";
	# find the complementary clades
	my @rootnode = grep {not exists $y2o->{$_}} keys %$o2y;
	my $rootnode = $rootnode[0];
	if (scalar(@rootnode) != 1) {
		die "bad root nodes";
	}
	my @alltips = split /\+/,$rootnode;
	my %all_clade;
	foreach my $internal_node (keys %$o2y) {
#print STDERR "$internal_node\n";
		my @tips = split /\+/,$internal_node;
		my %tips;
		$tips{$_} = 1 foreach @tips;
		my @outside_tips;
		foreach (@alltips) {
			if (not exists $tips{$_}) {
				push @outside_tips,$_;
			}
		}
		$all_clade{$internal_node} = \%tips;

		if (@outside_tips>0) {
			my $other_clade = join('+',sort {$a cmp $b} @outside_tips);
			my $tmp_hash = {};
			$tmp_hash->{$_}=1 foreach @outside_tips;
			$all_clade{$other_clade} = $tmp_hash;
		}
	}
	for my $q (0..$#alltips) {
		$all_clade{$alltips[$q]} = {$alltips[$q]=>1};
		
		my @tmptip = @alltips;
		splice @tmptip,$q,1;
		my $other_clade = join('+',sort {$a cmp $b} @tmptip);
		my $tmp_hash = {};
		$tmp_hash->{$_}=1 foreach @tmptip;
		$all_clade{$other_clade} = $tmp_hash;
	}
#print STDERR join("\t",keys %all_clade)."\n";
#print STDERR join("\t",values %all_clade)."\n";
#print STDERR "\n\n";
#foreach my $internal_node (keys %$o2y) {
#print STDERR "int node 0 $internal_node\n";
#}
#exit 1;
#print STDERR join("d\td",keys %$o2y)."\n\n";
	return [$y2o,$o2y,\%all_clade,$rootnode,\%originalCladeName2name,\%name2originalCladeName];
}

sub f_read_newick_file {
	my ($tree_file) = @_;
	open FIN,"$tree_file" or die "can't read $tree_file";
	my $line = <FIN>;
	close FIN;
	chomp $line;
	my @a = split /\s+/,$line;
	my $curnewick = $a[-1];
	$curnewick =~ s/\;//g;
	return $curnewick;
}



sub f_read_file_pairName_pairSNP {
	my ($fileWithName,$fileNSNP,$maxNSNP) = @_;
	
	my @line2strainPair;
	my $finName;
	my %strain;
	if ($fileWithName=~/\.xz/) {
		open $finName,"xzcat $fileWithName |" or die "can't read file $fileWithName\n";
	} else {
		open $finName,"$fileWithName" or die "can't read file $fileWithName\n";
	} 
	while (<$finName>) {
		chomp;
		my ($sn1,$sn2,@a) = split;
		push @line2strainPair,[$sn1,$sn2];
		$strain{$sn1} = 1;
		$strain{$sn2} = 1;
	}
	close $finName;
	my @strain = sort {$a cmp $b} keys %strain;
	for my $q (0..$#strain) {
		my $name = $strain[$q];
		$strain{$name} = $q;
	}
	
	my $strainPair2SNPdist = {};
	my $strainPair2lineNumber = {};
	my @data_SNP;
	my @data_SNP_dist;
	my $maxDelta = 0;
	my $max_ddelta = 0;
	my $finNSNP;
	if ($fileNSNP=~/\.xz$/) {
		open $finNSNP,"xzcat $fileNSNP |" or die "can't read $fileNSNP\n";	
	} else {
		open $finNSNP,"<$fileNSNP" or die "can't read $fileNSNP\n";
	}
	my $q = -1;
	while (<$finNSNP>) {
		$q++;
		chomp;
		my @a = split;
		my $Delta = (eval join '+',@a) / scalar(@a);
		$maxDelta = $Delta if $maxDelta<$Delta;
		$max_ddelta = max($max_ddelta,@a);
		my $dist = f_SNP2dist($maxNSNP,@a);
		push @data_SNP,\@a;
		push @data_SNP_dist,$dist;
		my ($strain1,$strain2) = @{$line2strainPair[$q]};
		$strainPair2SNPdist->{$strain1}{$strain2} = $dist;
		$strainPair2SNPdist->{$strain2}{$strain1} = $dist;
		$strainPair2lineNumber->{$strain1}{$strain2} = $q;
		$strainPair2lineNumber->{$strain2}{$strain1} = $q;
	}
	close $finNSNP;
	
	return [\@strain,\@line2strainPair,\@data_SNP,\@data_SNP_dist,$strainPair2SNPdist,$maxDelta,$strainPair2lineNumber,$max_ddelta];
}



sub f_SNP2dist {
	my ($maxNSNP,@NSNP) = @_;
	my @dist = (0) x (1+$maxNSNP);
	foreach (@NSNP) {
		my $aa = ($_>$maxNSNP) ? $maxNSNP : $_;
		$dist[$aa]++;
	}
	return \@dist;
}



sub f_cut_and_graft {
	my ($clusterData,$cluster2mv,$distmat) = @_;
	die "erereooasdfadvcxvrr\n" if not exists $clusterData->{$cluster2mv};
	return [] if not exists $clusterData->{$cluster2mv}{upstreamCluster};
	my $root2tip = {};
	my $tip2root = {};
	my @output_dist;
	foreach my $cluster (keys %$clusterData) {
		if (exists $clusterData->{$cluster}{downstreamCluster}) {
			my @dsnodes = keys %{$clusterData->{$cluster}{downstreamCluster}};
			foreach my $dsn (@dsnodes) {
				$root2tip->{$cluster}{$dsn} = 1;
				$tip2root->{$dsn} = $cluster;
			}
		}
	}
	
	my $f_node2tip;
		$f_node2tip = sub {
			my ($curnode,$node2tip,$tmp_root2tip) = @_;
			if (exists $node2tip->{$curnode}) {
				return [keys %{$node2tip->{$curnode}}];
			} elsif (exists $distmat->{$curnode}) {
				$node2tip->{$curnode}{$curnode} = 1;
				return [$curnode];
			} elsif (not(exists $tmp_root2tip->{$curnode}) || (scalar(keys %{$tmp_root2tip->{$curnode}})==0)) {
				$node2tip->{$curnode} = {};
				return [];
			} else {
				my @dsnodes = keys %{$tmp_root2tip->{$curnode}};
				my @tmptip;
				for my $ds (@dsnodes) {
					my $tip1 = $f_node2tip->($ds,$node2tip,$tmp_root2tip);
					push @tmptip,@$tip1;
				}
				for (@tmptip) {
					$node2tip->{$curnode}{$_} = 1;
				}
				return \@tmptip;
			}
		};
		
	my $cluster2mv_height = $clusterData->{$cluster2mv}{height};
	my $cluster2mv_upnode = $clusterData->{$cluster2mv}{upstreamCluster};
	delete $root2tip->{$cluster2mv_upnode}{$cluster2mv};
	delete $tip2root->{$cluster2mv};
	my @clusterMvTo0 = grep {exists($clusterData->{$_}{downstreamCluster}) && not(exists($clusterData->{$_}{downstreamCluster}{$cluster2mv}))} grep {$clusterData->{$_}{height}>$cluster2mv_height} keys %$clusterData;
	for my $targetCluster (@clusterMvTo0) {
		my $tmp_root2tip = dclone $root2tip;
		my $tmp_tip2root = dclone $tip2root;
		$tmp_root2tip->{$targetCluster}{$cluster2mv} = 1;
		$tmp_tip2root->{$cluster2mv} = $targetCluster;
		
		my $node2tip = {};
		my @rootnode = grep {not exists $tmp_tip2root->{$_}} keys %$tmp_root2tip;
		die "werwererer\t".scalar(@rootnode)."\n" if @rootnode != 1;
		my @tttmp = $f_node2tip->($rootnode[0],$node2tip,$tmp_root2tip);
		my $tmp_dist;
		my $Npair = 0;
		foreach my $usnode (keys %{$tmp_root2tip}) {
			my $tmp_height = $clusterData->{$usnode}{height};
			my @dsnode = keys %{$tmp_root2tip->{$usnode}};
			if (@dsnode<=1) {
				next;
			}
			for my $q (0..$#dsnode-1) {
				my @qtip = keys %{$node2tip->{$dsnode[$q]}};
				for my $r ($q+1..$#dsnode) {
					my @rtip = keys %{$node2tip->{$dsnode[$r]}};					
					foreach my $qqtip (@qtip) {
						foreach my $rrtip (@rtip) {
							$tmp_dist->{$qqtip}{$rrtip} = $tmp_height;
							$tmp_dist->{$rrtip}{$qqtip} = $tmp_height;
							$Npair++;
						}
					}
				}
			}
		}
		push @output_dist,$tmp_dist;
	}
	
	my $unodeHeight = $clusterData->{$cluster2mv_upnode}{height};
	my @clusterMvTo2 = grep {$_ ne $cluster2mv_upnode} grep {$clusterData->{$clusterData->{$_}{upstreamCluster}}{height}>$unodeHeight} grep {exists($clusterData->{$_}{upstreamCluster})} grep {$clusterData->{$_}{height}<$unodeHeight} keys %$clusterData;
	for my $targetLowCluster (@clusterMvTo2) {
		my $targetHighCluster = $clusterData->{$targetLowCluster}{upstreamCluster};
		my $tmp_root2tip = dclone $root2tip;
		my $tmp_tip2root = dclone $tip2root;
		delete $tmp_root2tip->{$targetHighCluster}{$targetLowCluster};
		
		my $tmpnode = 'tmpnodeaa';
		$tmp_root2tip->{$targetHighCluster}{$tmpnode} = 1;
		$tmp_tip2root->{$tmpnode} = $targetHighCluster;
		$tmp_root2tip->{$tmpnode}{$targetLowCluster} = 1;
		$tmp_tip2root->{$targetLowCluster} = $tmpnode;
		$tmp_root2tip->{$tmpnode}{$cluster2mv} = 1;
		$tmp_tip2root->{$cluster2mv} = $tmpnode;
		
		my @rootnode = grep {not exists $tmp_tip2root->{$_}} keys %$tmp_root2tip;
		die "werwerererwer\t".scalar(@rootnode)."\n" if @rootnode != 1;
		my $node2tip = {};
		my @tttmp = $f_node2tip->($rootnode[0],$node2tip,$tmp_root2tip);
		my $tmp_dist;
		my $Npair = 0;
		foreach my $usnode (keys %{$tmp_root2tip}) {
			my $tmp_height = (exists $clusterData->{$usnode}) ? $clusterData->{$usnode}{height} : $unodeHeight;
			my @dsnode = keys %{$tmp_root2tip->{$usnode}};
			if (@dsnode<=1) {
				next;
			}
			for my $q (0..$#dsnode-1) {
				my @qtip = keys %{$node2tip->{$dsnode[$q]}};
				for my $r ($q+1..$#dsnode) {
					my @rtip = keys %{$node2tip->{$dsnode[$r]}};					
					foreach my $qqtip (@qtip) {
						foreach my $rrtip (@rtip) {
							$tmp_dist->{$qqtip}{$rrtip} = $tmp_height;
							$tmp_dist->{$rrtip}{$qqtip} = $tmp_height;
							$Npair++;
						}
					}
				}
			}
		}
		push @output_dist,$tmp_dist;
	}
	
	
	my @clusterMvTo = grep {$_ ne $cluster2mv_upnode} grep {$clusterData->{$clusterData->{$_}{upstreamCluster}}{height}>$cluster2mv_height} grep {exists($clusterData->{$_}{upstreamCluster})} grep {$clusterData->{$_}{height}<$cluster2mv_height} keys %$clusterData;
	for my $targetLowCluster (@clusterMvTo) {
		my $targetHighCluster = $clusterData->{$targetLowCluster}{upstreamCluster};
		my $tmp_root2tip = dclone $root2tip;
		my $tmp_tip2root = dclone $tip2root;
		delete $tmp_root2tip->{$targetHighCluster}{$targetLowCluster};
		
		$tmp_root2tip->{$targetHighCluster}{$cluster2mv} = 1;
		$tmp_tip2root->{$cluster2mv} = $targetHighCluster;
		$tmp_root2tip->{$cluster2mv}{$targetLowCluster} = 1;
		$tmp_tip2root->{$targetLowCluster} = $cluster2mv;				
		
		my @rootnode = grep {not exists $tmp_tip2root->{$_}} keys %$tmp_root2tip;
		die "werwerererwer\t".scalar(@rootnode)."\n" if @rootnode != 1;
		my $node2tip = {};
		my @tttmp = $f_node2tip->($rootnode[0],$node2tip,$tmp_root2tip);
		my $tmp_dist;
		my $Npair = 0;
		foreach my $usnode (keys %{$tmp_root2tip}) {

			my $tmp_height = $clusterData->{$usnode}{height};
			my @dsnode = keys %{$tmp_root2tip->{$usnode}};
			if (@dsnode<=1) {
				next;
			}
			for my $q (0..$#dsnode-1) {
				my @qtip = keys %{$node2tip->{$dsnode[$q]}};
				for my $r ($q+1..$#dsnode) {
					my @rtip = keys %{$node2tip->{$dsnode[$r]}};					
					foreach my $qqtip (@qtip) {
						foreach my $rrtip (@rtip) {
							$tmp_dist->{$qqtip}{$rrtip} = $tmp_height;
							$tmp_dist->{$rrtip}{$qqtip} = $tmp_height;
							$Npair++;
						}
					}
				}
			}
		}
		push @output_dist,$tmp_dist;
	}
	
	return \@output_dist;
}





1;
