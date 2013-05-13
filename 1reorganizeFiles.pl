use warnings;
use strict;

my $dir = '/media/Elements/introgression/';
my $libInfo = $dir . 'libraryInfo.txt';

my %lib;
open(IN, "<$libInfo");
while(<IN>) {
	chomp(my $line = $_);
	my @d = split(/\t/,$line);
	#contact and then library
	$lib{$d[2]}{$d[0]}++;
	}
close(IN);

foreach my $contact (keys %lib) {
	my $subdir = $dir . $contact . '/';
	mkdir($subdir) unless(-d $subdir);
	foreach my $lib (keys %{$lib{$contact}}) {
		my $subsubdir = $subdir . $lib . '/';
		my $fwdSeq = $subsubdir . $lib . "_1.fastq.gz";
		my $revSeq = $subsubdir . $lib . "_2.fastq.gz";
		mkdir($subsubdir) unless(-d $subsubdir);
		my $call1 = system("mv $dir" . "$contact/*$lib*gz $subsubdir");
		unless(-f $fwdSeq) {
			my $call2 = system("gzip -1 -dc $subsubdir*R1* | gzip -1 -c > $fwdSeq");
			my $call3 = system("rm $subsubdir*R1*");
			}
		unless(-f $revSeq) {
			my $call4 = system("gzip -1 -dc $subsubdir*R2* | gzip -1 -c > $revSeq");
			my $call5 = system("rm $subsubdir*R2*");
			}
		print "Finished $lib for $contact!\n";
		}
	}

