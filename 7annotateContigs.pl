use strict;
use warnings;

my $dir = '/Users/singhal/thesisWork/introgression/';
my @contacts = qw(carlia sjo nBMB gillies);
my @names = qw(Carlia Sapro Lampro1 Lampro2);

for (my $i = 0; $i < scalar(@contacts); $i++) {
	my $contact = $contacts[$i];
	my $name = $names[$i];
	my $seq = $dir . 'fullTranscripts/' . $contact . '.fa';
	my $contig = $dir . 'targetSequences/final/' . $contact . '_targets.fa.final'; 
	my $out = $contig;
	$out =~ s/\.final/\.annotated/;

	my $s = parseSeq($seq);	
	my $c = parseContig($contig);	
	$c = annotateContigs($s,$c,$name);	
	printContigs($c,$out);
	}
	
sub printContigs {
	my ($c,$out) = @_;
	
	open(OUT, ">$out");
	
	foreach my $contig (keys %{$c}) {
		if ($c->{$contig}->{'info'}) {
			print OUT ">", $contig, "\t", $c->{$contig}->{'info'}, "\n", $c->{$contig}->{'seq'}, "\n";
			}
		else {
			#do ORF finding!
			my $info = ORFfinder($c->{$contig}->{'seq'});
			print OUT ">", $contig, "\t", $info, "\n", $c->{$contig}->{'seq'}, "\n";
			}
		}
	
	close(OUT);
	}	
	
sub ORFfinder {
	my ($seq) = @_;
	my %match;
	for (my $i = 0; $i < 3; $i++) {
		my $subseq = substr $seq, $i;
		my $aa = translate($subseq);
		my $match = ''; my $location;
		while ($aa =~ m/([A-Z]+)/g) {
			$match = $1 if length($1) > length($match);
			}
		if ($aa =~ m/($match)/) {	
			$location = $-[0];
			}
		$match{$i}{'length'} = length($match) * 3;
		$match{$i}{'loc'} = ($location + 1) * 3 - 2 + $i;
		}			
	my $rseq = $seq;	
	$rseq = reverse($rseq);
	$rseq =~ tr/ATGC/TACG/;		
	for (my $i = 0; $i < 3; $i++) {
		my $subseq = substr $rseq, $i;
		my $aa = translate($subseq);
		my $match = ''; my $location;
		while ($aa =~ m/([A-Z]+)/g) {
			$match = $1 if length($1) > length($match);
			}
		if ($aa =~ m/($match)/) {	
			$location = $-[0];
			}
		$match{'r' . $i}{'length'} = length($match) * 3;
		my $tmp = ($location + 1) * 3 - 2 + $i;
		$tmp = length($seq) - $tmp;
		$match{'r' . $i}{'loc'} = $tmp - length($match) * 3 + 2;
		}		
	 #winner is $keys[0]; does not account for multiple good matches
	 my @keys = sort {$match{$b}{'length'} <=> $match{$a}{'length'}} keys %match;
	 my $info;
	 if ($keys[0] =~ m/r/) {
		my $gs = $match{$keys[0]}{'loc'};
		my $ge = $match{$keys[0]}{'length'} + $gs - 1;
		$info = 'R_gs' . $gs . '_ge' . $ge;
		}
	else {
		my $gs = $match{$keys[0]}{'loc'};
		my $ge = $match{$keys[0]}{'length'} + $gs - 1;
		$info = 'gs' . $gs . '_ge' . $ge;		
		}				
	return($info);
	}
	
sub annotateContigs {
	my ($s,$c,$name) = @_;
	
	foreach my $contig (keys %{$c}) {
		my $cout = $contig . "_contig.fa";
		my $gene = $1 if $contig =~ m/^(ENS[A-Z|0-9]+)/;
		my $gout = $gene . "_gene.fa";
		open(COUT, ">$cout"); open(GOUT, ">$gout");
		print COUT ">", $contig, "\n", $c->{$contig}->{'seq'}, "\n";
		print GOUT ">", $gene, "\n", $s->{$gene}, "\n";		
		close(COUT); close(GOUT);
				
		my @call = `exonerate -m protein2genome -t $cout -q $gout --showtargetgff TRUE --showalignment NO --showvulgar 0 -n 1 | grep 'cds' | sort -k 4,4 -n`;
						
		if (@call) {
			my $info;
			foreach my $line (@call) { 
				my @d = split(/\s+/,$line);
				
				#identify gene start and stop
				my $tmpinfo = 'gs' . $d[3] . '_ge' . $d[4];
				
				if ($info) {
					$info .= '_' . $tmpinfo;
					}
				else {
					$info = 'R_' if $d[6] =~ m/-/;
					$info .= $tmpinfo;
					}
				}				
			$c->{$contig}->{'info'} = $info;
			}
		#a blast match but no exonerate match?	
		else {
			my $tc = $contig; $tc =~ s/_\d+//;
			my $aln = $dir . 'probeDesign/results_' . $name . '/' . $tc . ".fa.aln";
			open(IN, "<$aln");
			my %tmp; my $tmpc;
			while(<IN>) {
				chomp (my $line = $_);
				if ($line =~ m/>(\S+)/) {
					$tmpc = $1;
					}
				else {
					$line =~ s/-//g;
					$tmp{$tmpc} .= $line;
					}
				}
			close(IN);
			
			my @tmp = sort{length($tmp{$b}) <=> length($tmp{$a})} keys %tmp;
			my $winner = $tmp{$tmp[0]};
			my %match;
			for (my $i = 0; $i < 3; $i++) {
				my $subseq = substr $winner, $i;
				my $aa = translate($subseq);
				if ($aa =~ m/([A-Z]+)/) {
					$match{$i} = length($1);
					}			
				}
			my @keys = sort {$match{$b} <=> $match{$a}} keys %match;
			$winner = substr $winner, $keys[0];
			$winner = translate($winner);
			
			open(GOUT, ">$gout");
			print GOUT ">", $gene, "\n", $winner, "\n";		
			close(GOUT);
			
			my @call = `exonerate -m protein2genome -t $cout -q $gout --showtargetgff TRUE --showalignment NO --showvulgar 0 -n 1 | grep 'cds' | sort -k 4,4 -n`;
						
			if (@call) {
				my $info;
				foreach my $line (@call) { 
					my @d = split(/\s+/,$line);
					
					#identify gene start and stop
					my $tmpinfo = 'gs' . $d[3] . '_ge' . $d[4];
				
					if ($info) {
						$info .= '_' . $tmpinfo;
						}
					else {
						$info = 'R_' if $d[6] =~ m/-/;
						$info .= $tmpinfo;
						}
					}				
				$c->{$contig}->{'info'} = $info;	
				}
			}	
		unlink($cout); unlink($gout);	
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
			$seq{$c}{'seq'} = $seq if $c =~ m/ENS/;
			}
		}
	close(IN);
	
	return(\%seq);
	}	

sub parseSeq {
	my ($s) = @_;
	
	my %seq;
	open(IN, "<$s");
	while(<IN>) {
		chomp(my $line = $_);
		if ($line =~ m/>(\S+)/) {
			my $c = $1;
			my $start = $1 if $line =~ m/gs(\d+)/;
			my $end = $1 if $line =~ m/ge(\d+)/;			
			my $match = $1 if $line =~ m/(ENS\S+)/;
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
	
	open(IN, "</Users/singhal/thesisWork/introgression/probeDesign/genomes/anolisProt.fa");
	while(<IN>) {
		chomp(my $line = $_);
		if ($line =~ m/>(\S+)/) {
			my $c = $1;
			chomp(my $seq = <IN>);
			unless($seq{$c}) {
				$seq{$c} = $seq;
				}
			}
		}
	close(IN);	
	
	return(\%seq);
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
	