use warnings;
use strict;

my @contacts = qw(gillies nBMB sjo carlia);
my @lib = qw(10kN 2kN 1kN nTail center sTail 1kS 2kS 10kS);
my $dir = '/media/Elements/introgression/';

foreach my $contact (@contacts) {
    my $files = makeFiles($dir,$contact,\@lib);
    my $out = runPile($dir,$contact,$files);
    my $snp = SNPS($dir,$contact);
    parsePile($dir,$contact,\@lib,$out,$snp);
}

sub SNPS {
    my ($dir,$contact) = @_;

    my $vcf = $dir . 'ancMapping/' . $contact . ".vcf";
    my %snp;
    open(IN, "<$vcf");
    while(<IN>) {
	chomp(my $line = $_);
	my @d = split(/\t/,$line);
	if (scalar(@d) > 4) {
	    $snp{$d[0]}{$d[1]}{'ref'} = $d[3];
	    $snp{$d[0]}{$d[1]}{'alt'} = $d[4];
	}
    }
    close(IN);

    return(\%snp);
}

sub parsePile {
    my ($dir,$contact,$lib,$in,$snp) = @_;

    my $out = $dir . 'benchtruth/'. $contact . '.af.out'; 
    my @pop = @{$lib};
    open(IN, "<$in");
    open(OUT, ">$out");
    while(<IN>) {
	chomp(my $line = $_);
	my @d = split(/\t/,$line);
    
	my $j = 0;
    
	for (my $i = 4; $i <= scalar(@d); $i += 3) {

	    while ($d[$i] =~ m/(\+\d+|\-\d+)/g) {
		my $match = $1;
		my $num = $1 if $match =~ m/(\d+)/;
		$match = '\\' . $match;
		$d[$i] =~ s/$match[atgcn]{$num}//i;
	    }

	    my $af = 0; my $other = 0;
	    $af++ while $d[$i] =~ m/\./g;
	    $af++ while $d[$i] =~ m/\,/g;

	    my $var = $snp->{$d[0]}->{$d[1]}->{'alt'};
	    $other++ while $d[$i] =~ m/$var/gi;
	   
	    my $tot = $af + $other;

	    $af = $af / ($af + $other);
	    $af = sprintf("%.03f",$af);

	    print OUT $d[0], "\t", $pop[$j], "\t", $af, "\n";
	    $j++;
	}
    }
    close(OUT); close(IN);
}

sub runPile {
    my ($dir,$contact,$files) = @_;

    my $bed = $dir . 'benchtruth/SNPs/' . $contact . "_UTR.bed"; 
    my $out = $dir . 'benchtruth/' . $contact . '.mpileup.out';
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
