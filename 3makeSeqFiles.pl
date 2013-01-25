##################################################################################
# a script to prepare sequence files to use for assembly		                 #
# external dependencies: none													 #
# written by Sonal Singhal, sonal.singhal1 [at] gmail.com, 15 January 2013       #
##################################################################################

use warnings;
use strict;

#home directory with all my files; make sure to end with a backslash
my $dir = '/home/singhal/introgression2/';
#this file has all the info about my library, and creates a directory of all my files
my $lib = $dir . 'libraryInfo.txt';
#how many libraries do you have per lineage? will not make these files unless all the libraries are there.
my $numLib = 9; 

#defines the library
open(IN, "<$lib");
my %lib;
while (<IN>) {
	chomp(my $line = $_);
	my @d = split(/\t/, $line);
	my $lib = $d[2] . '/' . $d[0] . '/';
	$lib{$d[2]}{$d[0]} = $lib; 
	}
close(IN);

###########################
# run the subroutines     #
###########################

foreach my $lib (keys %lib) {
	my (@p1, @p2, @un);
	foreach my $name (keys %{$lib{$lib}}) {
		my $subdir = $dir . $lib{$lib}{$name};
		if (-d $subdir) {
			push(@p1, $subdir . $name . '_1_final.fastq.gz') if (-f $subdir . $name . '_1_final.fastq.gz');
			push(@p2, $subdir . $name . '_2_final.fastq.gz');
			push(@un, $subdir . $name . '_u_final.fastq.gz');
			}
		}	
	if (scalar(@p1) == $numLib) {
		print "Making files in $lib!\n";
		makeCombo(\@p1,\@p2,\@un,$lib);
		}
	}

###########################
# behold the subroutines! #
###########################

sub makeCombo {
	my ($p1,$p2,$un,$lib) = @_;
	
	my $sequ = $dir . $lib . '/' . $lib . '_u.fastq.gz';	
	my $files3 = join(" ", @{$un});
	my $call1 = system("gzip -dc -1 $files3 | gzip -c -1 > $sequ") unless (-f $sequ);	
	
	my $seq1 = $dir . $lib . '/' . $lib . '_1.fastq';
	my $seq2 = $dir . $lib . '/' . $lib . '_2.fastq';
	
	unless(-f $seq1 . '.gz') {
	open(OUT, ">$seq1");
	foreach my $file (@{$p1}) {
		my $call = system("gunzip $file"); $file =~ s/\.gz//;
		open(IN, "<$file");
		while(<IN>) {
			chomp(my $seqhead = $_);
			my $seq = <IN>; my $qualid = <IN>; my $qual = <IN>;
			print OUT $seqhead . '/1' . "\n" . $seq . $qualid . $qual;
			}
		close(IN);
		my $call2 = system("gzip -1 $file");
		 }
		 close(OUT);
}
	
	unless(-f $seq2 . '.gz') {
	open(OUT, ">$seq2");
	foreach my $file (@{$p2}) {
		my $call = system("gunzip $file"); $file =~ s/\.gz//;
		open(IN, "<$file");
		while(<IN>) {
			chomp(my $seqhead = $_);
			my $seq = <IN>; my $qualid = <IN>; my $qual = <IN>;
			print OUT $seqhead . '/2' . "\n" . $seq . $qualid . $qual;
			}
			close(IN);
			my $call2 = system("gzip -1 $file");
		}
	 close(OUT);	
}
	
	my $call2 = system("gzip $seq1"); my $call3 = system("gzip $seq2");
	}
