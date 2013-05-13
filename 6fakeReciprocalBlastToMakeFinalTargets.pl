#############################################################################################
# a script to find assembled contigs that match to a given target, accounting for the fact  #
# that one target might have multiple contigs that match and that one contig might have     # 
# multiple targets that match                                                               #
#                                                                                           #
# requires blastn & formatdb & cd-hit-est, all in path                                      #              
#                                                                                           #
# v 1.6, fixed error in how blast output file is defined                                    #               
# written by Sonal Singhal, sonal [dot] singhal 1 [at] gmail [dot] com, 23 Jan 2013         #
#############################################################################################

use warnings;
use strict;
use List::Util qw(max min);

my $seqfile = 'refloci.fa';
my $assembly = 'JP1428_assembly.cleaned.fa.final';
my $finalSeq = 'JP1428_finalref2.fa';
my $errorFile = 'JP1428_notcovered.out';

my $cluster = 0.99; #how much should you cluster the targets and assemblies at the get go
my $maxOverlap = 0.5; #how much overlap is allowed between adjoining assembled contigs mapping to the same target
my $perMatch = 80; #how similar does the assembled contig have to be to the target (note this is out of 100)
my $eval = 1e-10; #used in the initial BLAST step

open(SEQ, ">$finalSeq");
open(ERR, ">$errorFile");

my $call1 = system("cd-hit-est -i $seqfile -o tmp -c $cluster");
my $call2 = system("mv tmp $seqfile");
my $call3 = system("cd-hit-est -i $assembly -o tmp -c $cluster");
my $call4 = system("mv tmp $assembly");
my $call5 = system("rm tmp*");
my $call6 = system("formatdb -i $seqfile -p F");
my $out = $assembly . "_" . $1 . ".blast.out" if $seqfile =~ m/\/([^\/]+)$/;
my $call7 = system("blastn -query $assembly -db $seqfile -out $out -evalue $eval -outfmt 6");
my $call8 = system("rm formatdb.log");

my %tmp;
open(IN, "<$out");
while (<IN>) {
	chomp(my $line = $_);
	my @d = split(/\t/,$line);
	push(@{$tmp{$d[0]}},\@d) if $d[2] >= $perMatch;
	}
close(IN);

my %matches;
foreach my $id (keys %tmp) {	
    my $mArray = removeOverlap1($tmp{$id});
    my @mArray = @{$mArray};
    for (my $i = 0; $i < scalar(@mArray); $i++) {
	push(@{$matches{$mArray[$i][1]}}, \@{$mArray[$i]});
	}
    }	
undef %tmp;	
			
foreach my $id (keys %matches) {
	my $mArray = removeOverlap(\@{$matches{$id}});
	$matches{$id} = $mArray;
	}
		
my $seq = readSeq($seqfile);
my %seq = %{$seq};

my $contigs = readSeq($assembly);
my %contigs = %{$contigs};	
	
my %print;	
foreach my $id (keys %seq) {
	if ($matches{$id}) {
		my %length;
		for (my $i = 0; $i < scalar(@{$matches{$id}}); $i++) {
			my $start = $matches{$id}[$i][8];
			my $end = $matches{$id}[$i][9];
			for (my $n = min($start,$end); $n <= max($start,$end); $n++) {
				$length{$n}++;
				}			
			$print{$matches{$id}[$i][0]}{$id}++;
			}
		my $overlap = sprintf("%.3f",scalar(keys %length) / length($seq{$id}));	
		print ERR $id, "\t", $overlap, "\n";	
		}
	else {
		print SEQ ">", $id, "\n", $seq{$id}, "\n";
		print ERR $id, "\t", "NA\n";
		}
	}	

my %ids;	
foreach my $c (keys %print) {
	my $newid = join("_", keys %{$print{$c}}) . "_1";
	if ($ids{$newid}) {
		$ids{$newid}++;	
		if ($newid =~ m/(\S+)_(\d*)/) {
			my $core = $1;
			my $no = $ids{$newid};
			$newid  = $core . '_' . $no;
			}
		}
	else {
		$ids{$newid}++;
		}	
	print SEQ ">", $newid, "\n", $contigs{$c}, "\n"; 
	}
close ERR; close SEQ;		
	
#curse the recursive function!
sub removeOverlap {
	my ($array) = @_;
	
	for (my $i = 0; $i < scalar(@$array); $i++) {
		my $start1 = $array->[$i]->[8];
		my $end1 = $array->[$i]->[9];
		my %bp;
		for (my $n = min($start1,$end1); $n <= max($start1,$end1); $n++) {
			$bp{$n}++;
			}
		for (my $j = $i+1; $j < scalar(@$array); $j++) { 	
			my $start2 = $array->[$j]->[8];
			my $end2 = $array->[$j]->[9];
			my $overlap = 0;
			for (my $n = min($start2,$end2); $n <= max($start2,$end2); $n++) {
				$overlap++ if $bp{$n};
				}
			$overlap = $overlap / min(abs($start1 - $end1),abs($start2 - $end2));	
			if ($overlap > $maxOverlap) {
				if (abs($start1 - $end1) < abs($start2 - $end2)) {
					splice(@$array,$i,1);							
					}
				else {
					splice(@$array,$j,1);
					}	
				removeOverlap($array);	
				}
			}
		}	
	return($array);	
	}	
	
#curse the recursive function!
sub removeOverlap1 {
	my ($array) = @_;
	
	for (my $i = 0; $i < scalar(@$array); $i++) {
		my $start1 = $array->[$i]->[6];
		my $end1 = $array->[$i]->[7];
		my %bp;
		for (my $n = min($start1,$end1); $n <= max($start1,$end1); $n++) {
			$bp{$n}++;
			}
		for (my $j = $i+1; $j < scalar(@$array); $j++) { 	
			my $start2 = $array->[$j]->[6];
			my $end2 = $array->[$j]->[7];
			my $overlap = 0;
			for (my $n = min($start2,$end2); $n <= max($start2,$end2); $n++) {
				$overlap++ if $bp{$n};
				}
			$overlap = $overlap / min(abs($start1 - $end1),abs($start2 - $end2));	
			if ($overlap > $maxOverlap) {
				if (abs($start1 - $end1) < abs($start2 - $end2)) {
					splice(@$array,$i,1);							
					}
				else {
					splice(@$array,$j,1);
					}	
				removeOverlap1($array);	
				}
			}
		}	
	return($array);	
	}	
	
sub readFile {
	my ($file,$hash,$base) = @_;
	if (-f $file) {
		open(IN, "<$file");	
		my $id; my $tracker = 1;
		while(<IN>) {
			chomp(my $line = $_);
			if ($line =~ m/>(\S+)/) {
				$id = $base . $tracker;
				$tracker++;
				}
			else {
				$hash->{$id} .= $line;
				}
			}
		close(IN);	
		}	
	return($hash);
	}	
	
sub readSeq {
	my ($seqfile) = @_;
	my %seq; my $id;
	open(IN, "<$seqfile");
	while(<IN>) {
		chomp(my $line = $_);
		if ($line =~ m/>(\S+)/) {
			$id = $1;
			}
		else {
			$seq{$id} .= $line;
			}
		}
	close(IN);
	return(\%seq);
	}
