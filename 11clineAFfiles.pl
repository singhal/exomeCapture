use warnings;
use strict;

my $mincov = 50; #only consider allele frequency for a population if coverage is this high
my $snpqual = 20; #only want to consider SNPs with at least this much quality or higher; 999 is the max
my @contacts = qw(sjo carlia gillies nBMB);
my @lib = qw(10kN 2kN 1kN nTail center sTail 1kS 2kS 10kS);
my $dir = '/media/Elements/introgression/';

foreach my $contact (@contacts) {
    my ($bed,$snp,$cline) = makeBed($dir,$contact);
    my $files = makeFiles($dir,$contact,\@lib);
    my $out = runPile($dir,$contact,$files,$bed);
    $cline = parsePile($dir,$contact,\@lib,$out,$snp,$cline);
    printCline($dir,$contact,$cline);
}

sub printCline {
    my ($dir,$contact,$cline) = @_;
    
    my $out = $dir . 'ancMapping/' . $contact . '/' . $contact . '.cline.out';
    open(OUT, ">$out");
   
    my @tmplib = @lib;
    push(@tmplib,'ancS'); unshift(@tmplib,'ancN');

    my %cline = %{$cline};
    foreach my $c (sort {$a cmp $b} keys %cline) {
	foreach my $pos (sort {$a <=> $b} keys %{$cline{$c}}) {
	    my $diff1 = 0; my $diff2 = 0;
	    if (exists $cline{$c}{$pos}{'ancN'}) {
		$diff1 = abs($cline{$c}{$pos}{'ancN'} - $cline{$c}{$pos}{'ancS'});
	    }
	    if (exists $cline{$c}{$pos}{'10kN'}) {
				
		my $inverse = 0;
		if (exists $cline{$c}{$pos}{'ancN'}) {
		    $inverse = 1 if $cline{$c}{$pos}{'ancN'} >= 0.5;
		    }
		else {
		    $inverse = 1 if $cline{$c}{$pos}{'10kN'} >= 0.5;
		    }

		if ($inverse) {
            		    for (my $i = 0; $i < scalar(@tmplib); $i++) {
		       if (exists $cline{$c}{$pos}{$tmplib[$i]}) {
			  if ($cline{$c}{$pos}{$tmplib[$i]} eq 'NA') {
			      print OUT $c, "\t", $pos, "\t",$tmplib[$i], "\t", "NA", "\n";
			      }
			   else {
			      my $x = sprintf("%.3f",1-$cline{$c}{$pos}{$tmplib[$i]});
			      print OUT $c, "\t", $pos, "\t", $tmplib[$i], "\t", $x, "\n";
			      }
			   }
		        else {
			    print OUT $c, "\t", $pos, "\t",$tmplib[$i], "\t", "NA", "\n";
			    }
			}
		    }
		    else {
			for (my$i = 0; $i < scalar(@tmplib); $i++) {
                            if (exists $cline{$c}{$pos}{$tmplib[$i]}) {
				if ($cline{$c}{$pos}{$tmplib[$i]} eq 'NA') {
                                    print OUT $c, "\t", $pos, "\t",$tmplib[$i], "\t", "NA", "\n";
				}
				else {
                                    my $x = sprintf("%.3f",$cline{$c}{$pos}{$tmplib[$i]});
                                    print OUT $c, "\t", $pos, "\t", $tmplib[$i], "\t", $x, "\n";
				}
			    }   
                            else {
				print OUT $c, "\t", $pos, "\t",$tmplib[$i], "\t", "NA", "\n";
			    }   
                        }
		    }
		}
	    }
	}
    }

sub makeBed {
    my ($dir,$contact) = @_;
    
    my %snp; my %bed; my %cline;
    
    my $outdir = $dir . 'bed/';
    mkdir($outdir) unless (-d $outdir);
    my $out = $outdir . $contact . ".bed";
    open(OUT, ">$out");

    my $vcf1 = $dir . 'ancMapping/' . $contact . ".vcf";
    open(IN, "<$vcf1");
    while(<IN>){
	chomp(my $line = $_);
	my @d = split(/\t/,$line);
	if (scalar(@d) > 10) {

	    $snp{$d[0]}{$d[1]}{'ref'} = $d[3];
            $snp{$d[0]}{$d[1]}{'alt'} = $d[4];

	    my $start = scalar(@d) - 10; my $af1 = 0; my $af2 = 0; 
	    for (my $i = $start; $i < $start + 5; $i++) {
		$af1++ if $d[$i] =~ m/1\//;
		$af1++ if $d[$i] =~ m/\/1/;
	    }
	    $af1 = $af1 / 10;
	    $start = scalar(@d) - 5;
	    for (my $i = $start; $i < $start + 5; $i++) {
		$af2++ if $d[$i] =~ m/1\//;
		$af2++ if $d[$i] =~ m/\/1/;
	    }
	    $af2 = $af2 / 10;
	    	    
	    $cline{$d[0]}{$d[1]}{'ancN'} = $af1;
	    $cline{$d[0]}{$d[1]}{'ancS'} = $af2;

	    unless ( ($af1 == 0 && $af2 == 0) || ($af1 == 1 && $af2 == 1) ) {
		$bed{$d[0]}{$d[1]}++;
	    }
	}
    }
    close(IN);

    my $vcf2 = $dir . 'ancMapping/' . $contact . "_HZ.vcf";
    open(IN, "<$vcf2");
    while(<IN>) {
	chomp(my $line = $_);
	my @d = split(/\t/,$line);
	if (scalar(@d) > 10) {
	    my $start = scalar(@d) - 9; my $af1 = 0; my $af2 = 0; 
 
	    $snp{$d[0]}{$d[1]}{'ref'} = $d[3];
            $snp{$d[0]}{$d[1]}{'alt'} = $d[4];
       
	    for (my $i = $start; $i < $start + 9; $i++) {
		$af1++ if $d[$i] =~ m/1\//;
		$af1++ if $d[$i] =~ m/\/1/;
		$af2++ if $d[$i] =~ m/0\//;
		$af2++ if $d[$i] =~ m/\/0/;
	    }
        
	    unless ($af1 == 18 || $af2 == 18) {
		$bed{$d[0]}{$d[1]}++;  
	    }
	}
    }
    close(IN);

    foreach my $c (sort {$a cmp $b} keys %bed) {
	foreach my $pos (sort {$a <=> $b} keys %{$bed{$c}}) {
	    print OUT $c, "\t", $pos - 1, "\t", $pos, "\n";
	}
    }
    close(OUT);
    return($out,\%snp,\%cline);
}

sub parsePile {
    my ($dir,$contact,$lib,$in,$snp,$cline) = @_;

    my @pop = @{$lib};
    open(IN, "<$in");
    while(<IN>) {
	chomp(my $line = $_);
	my @d = split(/\t/,$line);
    
	my $j = 0;
    
	for (my $i = 4; $i <= scalar(@d); $i += 3) {
	    my $var = $snp->{$d[0]}->{$d[1]}->{'alt'};

	    unless($var =~ m/,/) {

		while ($d[$i] =~ m/(\+\d+|\-\d+)/g) {
		    my $match = $1;
		    my $num = $1 if $match =~ m/(\d+)/;
		    $match = '\\' . $match;
		    $d[$i] =~ s/$match[atgcn]{$num}//i;
		}

		my $af = 0; my $other = 0;
		$af++ while $d[$i] =~ m/\./g;
		$af++ while $d[$i] =~ m/\,/g;
		
		$other++ while $d[$i] =~ m/$var/gi;

		my $tot = $af + $other;
		
		if ($tot <= $mincov) {
		    $af = 'NA';
		}
		else {
		    $af = $other / $tot;
		    $af = sprintf("%.03f",$af);
		}

		$cline->{$d[0]}->{$d[1]}->{$pop[$j]} =  $af;
		$j++;

	    }
	}
    }
    close(OUT); close(IN);
    return($cline);
}

sub runPile {
    my ($dir,$contact,$files,$bed) = @_;

    my $out = $dir . 'ancMapping/' . $contact . '/' . $contact . '_HZ.mpileup.out';
    my $target = $dir . 'targetSequences/final/' . $contact . "_targets.fa.final"; 
    
    my $filelist = join("\t",@{$files});

    unless(-f $out) {
	my $call = system("samtools mpileup -d 10000 -f $target -l $bed $filelist > $out");
    }
    return($out);
}

sub makeFiles {
    my ($dir, $contact, $lib) = @_;

    my @files;
    foreach my $l (@{$lib}) {
	my $file = $dir . 'ancMapping/' . $contact . '/' . $l . ".sorted.bam";
	push(@files,$file);
    }
    
    return(\@files);
}
