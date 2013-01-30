use warnings;
use strict;

##################################################################################
# a script to clean up reads (adaptor/duplicate/contamination removal), trimming #
# external dependencies: flash, trimmomatic, bowtie2, cutadapt, cope             #         
# written by Sonal Singhal, sonal.singhal1 [at] gmail.com, 28 Dec 2012           #
# version 1.5 -- made a little faster, still                                     #
##################################################################################

#add gzipping and change file structure

#assumes a library naming convention of ONLY letters & numbers with a _[1|2].fastq.gz ending
#assumes HiSeq reads; will work with MiSeq with a small modification.

#directory with all read files
my $homedir = '/home/singhal/introgression/';
#this file has a list of the indexed adaptor OR primer sequences for my library, as fasta format, with each sequence named INDEXNUM_whatever.
my $adaptorFile = $homedir . 'genomes/adaptorSeqs.fa';
#this file has information about the libraries, here the first column has the library names and the second column has the adaptor name the third column has the species name
my $libInfo = $homedir . 'test.txt';
#contaminant genome in fasta format; here I am using human
my $contam = $homedir . 'genomes/contaminants.fa';
#the universal adaptor sequence
my $uniad = 'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT';

#paths to different executables
my $trimmomatic = '/home/singhal/bin/trimmomatic-0.22.jar';
my $cutadapt = '/home/singhal/cutadapt-1.2.1/bin/cutadapt';
my $flash = 'flash';
my $cope = 'cope';
my $bowtie = 'bowtie2';
my $flashslow =  'flash_slower';

my $readLength = 100;
my $nper = 0.6; #will get rid of reads for which more than $nper of bases are NNs
my $aper = 0.5; #will get rid of reads with any runs of bases longer than $aper*$readLength
my $minLength = 36;

my %lib;
open(IN, "<$libInfo");
while(<IN>) {
    chomp(my $line = $_);
    my @d = split(/\t/,$line);
    $lib{$d[0]} = {'species' => $d[2], 'ad'=> $d[1]};
}
close(IN);

foreach my $lib (keys %lib) {
    my $outdir = $homedir . $lib{$lib}{'species'} . '/' . $lib . '/';
    my $file1 = $outdir . $lib . '_1.fastq.gz';
    my $file2 = $file1;
    $file2 =~ s/_1/_2/;
    my $species = $lib{$lib}{'species'};    

    $file1 = unzip($file1); $file2 = unzip($file2);
 
	my $start1 = time;	
	my $dup = $outdir . $lib . '.duplicates.out';
	duplicates($file1, $file2, $dup);
	my $time1 = int((time - $start1)/60);
	print "Found duplicates in $lib in $time1 minutes! Now this is something...\n";

	my $start2 = time;
	my $low = $outdir . $lib . '.lowComplexity.out';
	removeLowComplexity($file1,$low);
	removeLowComplexity($file2,$low);
	my $time2 = int((time - $start2)/60);
	print "Found low complexity reads in $lib in $time2 minutes! Whew, almost there!\n";

	my $start3 = time;
	my %reads = ('1' => $file1, '2' => $file2);
	my $ad = getAdaptors($lib,$adaptorFile,$libInfo);
	my @clean1;
	for my $dir (keys %reads) {
		my $out1 = trimmomatic($lib,$ad,$reads{$dir},$reads{$dir},'trim1');
		my $out2 = cutadapt($ad,$out1,$reads{$dir},'trim2');		
		my $out3 = cutadapt($ad,$out2,$reads{$dir},'trim3');
		my $final = $reads{$dir} . '_cleaned1';
		my $call = system("mv $out3 $final");
		unlink($out1,$out2);
		push(@clean1,$final);
		}
	my $trim1 = fixMatePair($lib,\%reads,\@clean1,"trim1");
	my $reads2 = mergeReads($lib,\%reads,$trim1,"trim2");
	my %reads2 = %{$reads2};
	my @clean2;
	for my $dir (keys %reads2) {
	    my $out1 = trimmomatic($lib,$ad,$reads2{$dir},$reads2{$dir},'trim1');
	    my $out2 = cutadapt($ad,$out1,$reads2{$dir},'trim2');
	    my $out3 = cutadapt($ad,$out2,$reads2{$dir},'trim3');
	    my $final = $reads2{$dir} . '_cleaned2';
	    my $call = system("mv $out3 $final");
	    unlink($out1,$out2);
	    push(@clean2,$final);
	    }
	my $trim3 = fixMatePair($lib,\%reads,\@clean2, "trim2");
	my $reads3 = mergeReads($lib,\%reads,$trim3,"trim3");
	my $reads4 = reallyMergeReads($lib,\%reads,$reads3,"trim4","trim5");
	my $time3 = int((time - $start3)/60);
	print "Trimmed and merged in $lib in $time3 minutes! Whew, almost there!\n";

	my $start4 = time;
	my $contaminants =  $outdir . $lib . '.contam.out';
	removeContamination($reads4,$contam,$contaminants);
	my $time4 = int((time - $start4)/60);
	print "Removed contamination in $lib in $time4 minutes! It's going...\n";	

	makeFinal($outdir,$lib,$reads4,$dup,$low,$contaminants);
    	compress($file1); compress($file2);

	my $call_rm1 = system("rm $outdir*$lib*clean*");
	my $call_rm2 = system("rm $outdir*$lib*trim*");
	}

sub unzip  {
    my ($file) = @_;
    my $call = system("gunzip $file");
    my $new = $file;
    $new =~ s/\.gz//;
    return($new);
}

sub compress {
    my ($file) = @_;
    my $call = system("gzip $file");
}

sub makeFinal {
    my ($outdir,$lib,$reads,$dup,$low,$contam) = @_;
    
    my %junk;

    open(IN, "<$dup");
    while(<IN>) {
		chomp(my $line = $_);
		$junk{$line}++;
    	}
    close(IN);
    open(IN, "<$low");
    while(<IN>){
		chomp(my $line = $_);
		$junk{$line}++;
    	}
    close(IN);
    open(IN, "<$contam");
    while(<IN>) {
		chomp(my $line = $_);
		$junk{$line}++;
    	}
    close(IN);

    my $new1 = $outdir . $lib . "_1_final.fastq";
    my $new2 = $outdir . $lib . "_2_final.fastq";
    my $newu = $outdir . $lib . "_u_final.fastq";
    my %new = ('1' => $new1, '2' => $new2, 'u' => $newu);
    my %reads = %{$reads};
    foreach my $type (keys %reads) {
		open(OUT, ">$new{$type}");
		open(IN, "<$reads{$type}");
		while(<IN>) {
		    chomp(my $line = $_);
		    if ($line =~ m/^@(\S+)/) {
				my $id = $1;
				$id =~ s/_\S+//;
				my $seq = <IN>; my $qualid = <IN>; my $qual = <IN>;
				unless($junk{$id}){
		    		print OUT "@", $id, "\n", $seq, $qualid, $qual;
					}
			}
		}
		close(IN); close(OUT);
    	}
    compress($new1); compress($new2); compress($newu);
	}		

sub getAdaptors {
	my ($lib,$adaptorFile,$libInfo) = @_;
	
	my %lib1;
	open(IN, "<$libInfo");
	while(<IN>) {
		chomp(my $line = $_);
		my @d = split(/\t/,$line);
		$lib1{$d[0]} = $d[1];
		}
	close(IN);
	
	my %adaptors;
	open(IN, "<$adaptorFile");
	while(<IN>) {
		chomp(my $line = $_);
		if ($line =~ m/>(\d+)_/) {
			my $bc = $1;
			chomp(my $seq = <IN>);
			$adaptors{$bc} = $seq;
			}
		}
	close(IN);
	
	my %ad = ("uni" => $uniad, "uni_rc" => rc($uniad), "index" => $adaptors{$lib1{$lib}}, "index_rc" => rc($adaptors{$lib1{$lib}}));
	
	return(\%ad);
	}

sub rc {
	my ($seq) = @_;
	my $rc = $seq;
	$rc = reverse($rc);
	$rc =~ tr/ATGCatgc/TACGtacg/;
	return($rc);
	}

sub removeContamination {
	my ($reads,$contam,$contaminants) = @_;
	
	unless (-f $contam . ".3.bt2") {
		my $bw2build = $bowtie . "-build";
		my $call1 = system("$bw2build $contam $contam");
		}
		
	my %reads = %{$reads};
	
	my $contamout1 = $reads{'1'} . ".contam.sam1";
	my $call2 = system("$bowtie -x $contam -1 $reads{'1'} -2 $reads{'2'} --fast -S $contamout1 --sam-nohead --sam-nosq");
	my $contamout2 = $reads{'1'} . ".contam.sam2";
	my $call3 = system("$bowtie -x $contam -U $reads{'u'} --fast -S $contamout2 --sam-nohead --sam-nosq");
	
	parseSAM($contaminants,$contamout1);
	parseSAM($contaminants,$contamout2);
	}

sub parseSAM {
	my ($contaminants,$contam) = @_;
	open(OUT, ">>$contaminants");
	open(IN, "<$contam");

	while(<IN>) {
	    chomp(my $line = $_);
	    my @d = split(/\t/,$line);
	    if ($d[2] !~ m/\*/) {
			if ($d[5] =~ m/^\d+M$/) {
			    my $md = $1 if $line =~ m/(MD\:Z\S+)/;
			    my @a = ($md =~ m/([ATGC])/g);
			    if (scalar(@a) < 2) {
					$d[0] =~ s/\:[1|2]$//;
					print OUT $d[0], "\n";
	    			}
				}
    		}
		}
	close(IN); close(OUT);	
	unlink($contam);
	}	



sub duplicates {
	my ($file1,$file2,$dup) = @_;	
	my $sorted1 = sortFile($file1);
	my $sorted2 = sortFile($file2);	
	my $dub_ref = getDuplicates($sorted1);
	removeDuplicates($dub_ref, $sorted2, $dup);
	unlink($sorted1); unlink($sorted2);
	}
	
sub sortFile {
	my ($file) = @_;
	my $out = $file . "2";
	open(OUT, ">$out");
	my $sorted = $file . ".sorted";
	open(IN, "<$file");
	while(<IN>){
		chomp(my $line = $_);
		if ($line =~ m/^(\@HS\S+)/) {
		    my $id = $1 if $line =~ m/^(\@HS\S+)/;
			chomp(my $seq = <IN>);
			chomp(my $qualID = <IN>);
			chomp(my $qual = <IN>);
			print OUT $id, "\t", $seq, "\t", $qualID, "\t", $qual, "\n";
			}
		}
	close(IN); close(OUT);
	my $call = system("sort -k 2,2 $out > $sorted");
	unlink($out);
	return($sorted);
	}

sub getDuplicates {
	my ($sorted) = @_;
	open(SORT, "<$sorted");
	my ($old_id, $old_seq);
	my $tracker = 1;
	my %dup;
	while(<SORT>){
		chomp(my $line = $_);
		my @d = split(/\t/, $line);
		if ($old_id) {
			unless ($old_seq eq $d[1]) {
				delete $dup{$tracker} if scalar(@{$dup{$tracker}}) < 2;
				$tracker++;			
				}
			}
		push(@{$dup{$tracker}}, $d[0]);
		$old_id = $d[0]; $old_seq = $d[1];
		}
	close(SORT);
	return(\%dup);
	}	

sub removeDuplicates {
	my ($dupref, $sorted2, $dup) = @_;
	my %dup = %$dupref;
	my ($old_id, $old_seq);	
	my %dup_rev;

	foreach my $tracker (keys %dup) {
		foreach my $id (@{$dup{$tracker}}) {
		    $id =~ s/\:[1|2]$//;
			$dup_rev{$id} = $tracker;
			}
		}

	open(SORT, "<$sorted2");
	open(DUP, ">$dup");
	while(<SORT>){
		chomp(my $line = $_);
		my @d = split(/\t/, $line);
		if ($old_id) {
		    $d[0] =~ s/\:[1|2]$//; $old_id =~ s/\:[1|2]$//;
			if ($dup_rev{$old_id} && $dup_rev{$d[0]}) {	
				#these are duplicates in reverse direction
				if ($old_seq eq $d[1]) {
				    #these are duplicates in the forward direction too	
				    if ($dup_rev{$old_id} eq $dup_rev{$d[0]}) {
					$d[0] =~ s/^\@//;
						print DUP $d[0], "\n";
						}
					}
				}
			}
		$old_id = $d[0]; $old_seq = $d[1];
		}
	close(SORT); close(DUP); undef(%dup); undef(%dup_rev);
	}	

sub removeLowComplexity {
	my ($file,$low) = @_;
	open(IN, "<$file");
	open(OUT, ">>$low");
	while(<IN>) {
		chomp(my $line = $_);		
		if ($line =~ m/^@(\S+)/) {
			my $id = $1;
			chomp(my $seq = <IN>);
			chomp(my $qualid = <IN>); chomp(my $qual = <IN>);
			my $n = int($nper*length($seq));
			my $a = int($aper*length($seq));
			my $ncounter = ($seq =~ m/N/g);
			if ($seq =~ m/[A]{$a}/i || $seq =~ m/[T]{$a}/i || $seq =~ m/[G]{$a}/i || $seq =~ m/[C]{$a}/i || $ncounter >= $n) {
				print OUT $id, "\n";
				}
			}
		}	
	close(IN); close(OUT);	
	}

sub mergeReads {
    my ($lib,$orig,$reads,$base) = @_;
    my %reads = %{$reads};

    my $newread1 = $orig->{'1'} . '_' . $base . '_p1';
    my $newread2 = $orig->{'2'} . '_' .$base .'_p2';
    my $newreadu = $orig->{'1'} . '_' .$base .'_u';

    my $call1 = system("$flash $reads{'1'} $reads{'2'} -m 5 -x 0.01 -o $lib");
    my $call2 = system("cat $reads{'u'} $lib\.extendedFrags.fastq > $newreadu");
    my $call3 = system("mv $lib\.notCombined_1.fastq $newread1");
    my $call4 = system("mv $lib\.notCombined_2.fastq $newread2");
    my $call5 = system("rm $lib\.extendedFrags.fastq $lib\.hist*");
 
    my %newreads = ('1' => $newread1,'2' => $newread2, 'u' => $newreadu);
    return(\%newreads);
}

sub reallyMergeReads {
	my ($lib,$orig,$reads,$base1,$base2) = @_;	
	
	my %reads =  %{$reads};	
	my $newread1 = $orig->{'1'} . '_' . $base1 . '_p1';
    	my $newread2 = $orig->{'2'} . '_' .$base1 .'_p2';
    	my $newreadu = $orig->{'1'} . '_' .$base1 .'_u';
	my $call1 = system("cope -a $reads{'1'} -b $reads{'2'} -o $lib\.copemerged -2 $newread1 -3 $newread2 -m 0 -l 5 -c 0.9");

	my $newerread1 = $orig->{'1'} . '_' . $base2 . '_p1';
    	my $newerread2 = $orig->{'2'} . '_' .$base2 .'_p2';
    	my $newerreadu = $orig->{'1'} . '_' .$base2 .'_u';
    	my %newerreads = ('1' => $newerread1,'2' => $newerread2, 'u' => $newerreadu);

	my $call2 = system("$flashslow $newread1 $newread2 -m 5 -x 0.01 -o $lib");
		
	my $call3 = system("cat $lib\.copemerged $lib\.extendedFrags.fastq $reads{'u'} > $newerreadu");
	my $call4 = system("mv $lib\.notCombined_1.fastq $newerread1");
	my $call5 = system("mv $lib\.notCombined_2.fastq $newerread2");
	my $call6 = system("rm $lib\.copemerged $lib\.extendedFrags.fastq $lib\.hist*");
	return(\%newerreads);
   	}
			
sub fixMatePair {
	my ($lib,$read,$readarray,$base) = @_;
	my @trim = @{$readarray};
	my %pair;	
	foreach my $reads (@trim) {
		open(IN, "<$reads");
		while(<IN>) {
			chomp(my $line = $_);
			if ($line =~ m/^(@\S+)/) {
			    $pair{$1}++;
			    chomp(my $seq = <IN>); chomp(my $qualid = <IN>); chomp(my $qual = <IN>);
				}
			}
		close(IN);	
		}
	my %reads = %{$read};
	my $out1 = $reads{'1'} . '_' . $base . '_p1';
	my $out2 = $reads{'2'} . '_' . $base . '_p2';
	my $outu = $reads{'1'} . '_' . $base . '_u';
	open(OUT1, ">$out1"); 
	open(OUT2, ">$out2"); 
	open(OUTU, ">$outu"); 
	my %newpairs = ('1' => $out1, '2' => $out2, 'u' => $outu);
	foreach my $reads (@trim) {
		open(IN, "<$reads");
		my $file = $1 if $reads =~ m/$lib\_(\d)/;
		while(<IN>) {
			chomp(my $line = $_);	
			if ($line =~ m/^(@\S+)/) {
				my $id = $1;
				my $seq = <IN>;
				my $qualid = <IN>;
				my $qual = <IN>;
				if ($pair{$id} == 2) {
					if ($file == 1) {
						print OUT1 $id . "\n" . $seq . $qualid . $qual;
						}
					else {
						print OUT2 $id . "\n" . $seq . $qualid . $qual;
						}
					}
				else {
					print OUTU $id . "\n" . $seq . $qualid . $qual;
					}
				}	
			}
		close(IN);	
		}	
	close(OUT1); close(OUT2); close(OUTU);	
	return(\%newpairs);
	}

sub trimmomatic {
	my ($lib,$ad,$in,$base,$suffix) = @_;
	my $out  = $base . '_' . $suffix;	
	my $adfile = $lib . "_adfile.fa";
	open(OUT, ">$adfile");
	foreach my $name (keys %{$ad}) {
		print OUT ">", $name, "\n", $ad->{$name}, "\n";
		}
	close(OUT);	
	my $call1 = system("java -classpath $trimmomatic org.usadellab.trimmomatic.TrimmomaticSE -phred33 $in $out ILLUMINACLIP:$adfile:2:40:15 SLIDINGWINDOW:4:20 MINLEN:$minLength LEADING:3 TRAILING:3");
	unlink($adfile);
	return($out);
	}

sub cutadapt {
	my ($ad,$in,$base,$suffix) = @_;
	my $out  = $base . '_' . $suffix;	
	my $curRead = $in;
	my $tracker = 1;
	my %ad = %{$ad};
	foreach my $key (keys %ad) {
		my $out = $curRead . $tracker;
		my $call = system("$cutadapt -b $ad{$key} -O 4 -n 5 -f fastq $curRead -o $out -m $minLength");
		unlink($curRead) unless($curRead eq $in);
		$curRead = $out;
		$tracker++;
		}
	my $call2 = system("mv $curRead $out");
	return($out);
	}

sub bowtie {
	my ($lib,$ad,$in,$base,$suffix) = @_;
	my $out  = $base . '_' . $suffix;	
	my $file = $lib. "_out.sam";
	
	my $adfile = $lib . "_adfile_norev.fa";
	open(OUT, ">$adfile");
	foreach my $name (keys %{$ad}) {
		print OUT ">", $name, "\n", $ad->{$name}, "\n" unless $name =~ m/rc/;
		}
	close(OUT);	
	
	my $bw2build = $bowtie . "-build";
	my $call1 = system("$bw2build $adfile $adfile");
	my $call2 = system("$bowtie --local -D 15 -R 2 -N 1 -L 10 -i S,1,0.75 -k 1 -x $adfile -U $in -S $file");
	my $call3 = system("rm $adfile" . "*");
	
	open(IN, "<$file");
	open(OUT, ">$out");

	while(<IN>) {
		chomp(my $line1 = $_);
		my @d1 = split(/\t/,$line1);
		
		my $seq1 = $d1[9];
		my $qual1 = $d1[10];
		
		unless($line1 =~ m/^@/) {
			if ($line1 !~ m/\*/) {			
				if ($d1[5] =~ m/^(\d+)S\d+M$/) {
					my $l = $1;
					$seq1 = substr $seq1, 0, $l;
					$qual1 = substr $qual1, 0, $l;
					}
				elsif ($d1[5] =~ m/^(\d+)M\d+S$/) {
					my $start = $1;
					$seq1 = substr $seq1, $start;
					$qual1 = substr $qual1, $start;
					}
				else {
					my @s;
					while ($d1[5] =~ m/(\d+)S/g) {
						push(@s,$1);
						}
					@s = sort {$a <=> $b} @s;
					if ($s[$#s] >= $minLength) {	
						if ($d1[5] =~ m/^(\S*)$s[$#s]/) {
							my $match = $1;
							my $length = $s[$#s];
							my $start = 0;
							while ($match =~ m/(\d+)/g) {
								$start += $1;
								}
							$seq1 = substr $seq1, $length;
							$qual1 = substr $qual1, $length;	
							}													
						}
					else {
						$seq1 = 'N'; $qual1 = 'N';
						}
					}
				}
			
			if (length($seq1) >= $minLength) {	
				print OUT "@" . $d1[0] . "\n" . $seq1 . "\n" . '+' . "\n" . $qual1 . "\n";	
				}
			}	
		}	
	unlink($file);
	close(IN); close(OUT);	
	return($out);
	}
