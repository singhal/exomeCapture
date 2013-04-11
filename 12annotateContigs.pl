use strict;
use warnings;
use List::Util qw(min max);

#main directory with all the files
my $dir = '/Users/singhal/Desktop/activeWork/annotateExome/';
#transcriptome assembly
my $seqTrans = $dir . 'chipmunkTranscript.fa';
#exome assembly
my $contig = $dir . 'finalExomeAssembly.fa';
#reference protein assembly from model organism
my $prot = $dir . 'mouseProtein.fa';
#number of processors you can use for blasting
my $np = 4;
#name of annotated exome assembly -- you can change if you want
my $out = $contig . ".annotated";

my ($s,$prot_tmp) = parseSeq($seqTrans,$prot);	
my $blast = blast($contig,$prot_tmp);
my $c = parseContig($contig);	
$c = annotateContigs($s,$c,$blast);	
printContigs($c,$out);

sub blast {
	my ($exome,$prot) = @_;
	
	my $out = $dir . "blastResult.out";
	my $call1 = system("formatdb -i $prot -p T") unless (-f $prot . ".phr");
	my $call2 = system("blastall -p blastx -d $prot -i $exome -m 8 -o $out -v 1 -b 1 -e 1e-10 -a $np") unless (-f $out);
	
	my %blast;
	open(IN, "<$out");
	while(<IN>) {
		chomp(my $line = $_);
		my @d = split(/\t/,$line);
		$blast{$d[0]} = $d[1] unless $blast{$d[0]};
		}
	close(IN);	

	return(\%blast);
	}
	
sub printContigs {
	my ($c,$out) = @_;
	
	open(OUT, ">$out");
	
	foreach my $contig (keys %{$c}) {
		print OUT ">", $contig, "\t", $c->{$contig}->{'info'}, "\n", $c->{$contig}->{'seq'}, "\n";
		}
	
	close(OUT);
	}	
		
sub annotateContigs {
	my ($s,$c,$blast) = @_;
	
	foreach my $contig (keys %{$c}) {	
		my $info; my $seq = $c->{$contig}->{'seq'};
		if ($blast->{$contig}) {
		
			my $cout = "contig.fa";
			my $gene = $blast->{$contig};
			my $gout = "gene.fa";
			open(COUT, ">$cout"); open(GOUT, ">$gout");
			print COUT ">", $contig, "\n", $c->{$contig}->{'seq'}, "\n";
			print GOUT ">", $gene, "\n", $s->{$gene}, "\n";		
			close(COUT); close(GOUT);
				
			my @call = `exonerate -m protein2genome -t $cout -q $gout --showtargetgff TRUE --showalignment NO --showvulgar 0 -n 1 | grep 'cds' | sort -k 4,4 -n`;
						
			if (@call) {		
				my (@gs, @ge); my $reverse = 0;
				foreach my $line (@call) { 
					my @d = split(/\s+/,$line);
					
					my ($gs,$ge);
					
					if ($d[6] =~ m/-/) {
						my $length = length($seq);
						$gs = $length - $d[4] + 1;
						$ge = $length - $d[3] + 1;
						
						push(@gs,$gs); push(@ge,$ge);
						$reverse++;												
						}
					else {
						$gs = $d[3]; $ge = $d[4];
						push(@gs,$gs); push(@ge,$ge);
						}
					}	
				
				if ($reverse) {
					$seq = reverse($seq);
					$seq =~ tr/atgcATGC/tacgTACG/;
					@gs = reverse(@gs);
					@ge = reverse(@ge);
					}
					
				$info = $gene;	
				for (my $x = 0; $x < scalar(@gs); $x++) { 	
					$info .= '_gs' . $gs[$x] . '_ge' . $ge[$x];
					}				
				}
			#a blast match but no exonerate match?	
			else {
				$info = $gene;
				}	
			unlink($cout); unlink($gout);	
			}
		else {
			my $cout = "contig.fa";
			open(COUT, ">$cout");
			print COUT ">", $contig, "\n", $c->{$contig}->{'seq'}, "\n";
			close(COUT);
			
			my $call = system("blat $seqTrans $cout blatTmp -out=blast8");
			open(IN, "<blatTmp");
			chomp(my $match = <IN>);
			
			if ($match) {
				my @d = split(/\t/,$match);
				if ($d[10] <= 1e-10) {
					my $u5 = $1 if $d[1] =~ m/gs(\d+)/;
					my $u3 = $1 if $d[1] =~ m/ge(\d+)/;
			
					my $max = max($d[8],$d[9]);
					my $min = min($d[8],$d[9]);
					$info = $1 if $d[1] =~ m/([A-Z]*ENS[0-9|A-Z]+)/;
					
					my $overlap = 0;
					for (my $i = $min; $i <= $max; $i++) {
						if ($i > $u5 && $i < $u3) {
							$overlap++;
							}
						}	
				
					if ( ($overlap / ($max - $min) ) < 0.2 ) {
						$info .= '_utr';
						}
					}	
				else {
					$info = 'NAannotation';
					}
				}
			else {
				$info = 'NAannotation';
				}
			
			unlink($cout); my $call2 = system("rm blatTmp");
			}
		$c->{$contig}->{'info'} = $info;
		$c->{$contig}->{'seq'} = $seq;
		}
	return($c);	
	}
	
sub parseContig {
	my ($s) = @_;
	
	my %seq;
	open(IN, "<$s");
	while(<IN>) {
		chomp(my $line = $_);
		if ($line =~ m/>(\S+)/) {
			my $c = $1;
			chomp(my $seq = <IN>);			
			$seq{$c}{'seq'} = $seq;
			}
		}
	close(IN);
	
	return(\%seq);
	}	

sub parseSeq {
	my ($s,$g) = @_;
	
	my %seq;
	open(IN, "<$s");
	while(<IN>) {
		chomp(my $line = $_);
		if ($line =~ m/>(\S+)/) {
			my $c = $1;
			my $start = $1 if $line =~ m/gs(\d+)/;
			my $end = $1 if $line =~ m/ge(\d+)/;			
			my $match = $1 if $line =~ m/([A-Z]*ENS[0-9|A-Z]+)/;
			chomp(my $seq = <IN>);
			
			my $length = $end - $start + 1;
			unless( $length % 3) {
				$length = 3 * int($length / 3);
				}
			
			$seq = substr $seq, $start - 1, $length;
			my $aa = translate($seq);
			
			$seq{$match} = $aa;
			}
		}
	close(IN);
	
	my %prot; my $id;
	open(IN, "<$g");
	while(<IN>) {
		chomp(my $line = $_);
		if ($line =~ m/>(\S+)/) {
			$id = $1;
			}
		else {
			$prot{$id} .= $line;
			}
		}	
	close(IN);
	
	foreach my $c (keys %prot) {
		unless ($seq{$c}) {
			$seq{$c} = $prot{$c};
			}
		}	
	
	my $tmp = $dir . "tmpProteinFile.fa";
	
	open(TMP, ">$tmp");
	foreach my $c (keys %seq) {
		print TMP ">", $c, "\n", $seq{$c}, "\n";
		}
	close(TMP);	
	
	return(\%seq,$tmp);
	}		
	
sub translate {
	my $string = shift;
	$string = uc($string);
	my @codons = $string =~ m/(\S\S\S)/g;
	my %codons = (	'ATG'=>'M','ACG'=>'T','CTG'=>'L','CCG'=>'P','GTG'=>'V','GCG'=>'A','TTG'=>'L','TCG'=>'S',
					'ATA'=>'I','ACA'=>'T','CTA'=>'L','CCA'=>'P','GTA'=>'V','GCA'=>'A','TTA'=>'L','TCA'=>'S',
					'ATC'=>'I','ACC'=>'T','CTC'=>'L','CCC'=>'P','GTC'=>'V','GCC'=>'A','TTC'=>'F','TCC'=>'S',
					'ATT'=>'I','ACT'=>'T','CTT'=>'L','CCT'=>'P','GTT'=>'V','GCT'=>'A','TTT'=>'F','TCT'=>'S',
					'AGG'=>'R','AAG'=>'K','CGG'=>'R','CAG'=>'Q','GGG'=>'G','GAG'=>'E','TGG'=>'W','TAG'=>'*',
					'AGA'=>'R','AAA'=>'K','CGA'=>'R','CAA'=>'Q','GGA'=>'G','GAA'=>'E','TGA'=>'*','TAA'=>'*',
					'AGC'=>'S','AAC'=>'N','CGC'=>'R','CAC'=>'H','GGC'=>'G','GAC'=>'D','TGC'=>'C','TAC'=>'Y',
					'AGT'=>'S','AAT'=>'N','CGT'=>'R','CAT'=>'H','GGT'=>'G','GAT'=>'D','TGT'=>'C','TAT'=>'Y');
	my $translate;
	foreach(@codons) {
		if ($codons{$_}) {
			$translate = $translate . $codons{$_};
			}
		else {
#			print "ERROR: ILLEGAL PASS TO CODON TRANSLATION: $_ is not a codon!\n";
			$translate = $translate . 'X';
			}
		}
	return($translate);
	}
	