use warnings;
use strict;

my %contacts = ('carlia' => {'1' => 'Carlia_N', '2' => 'Carlia_S'},
                'sjo' => {'1' => 'Sapro_C', '2' => 'Sapro_S'},
		'gillies' => {'1' => 'Lampro_C', '2' => 'Lampro_S'},
		'nBMB' => {'1' => 'Lampro_N', '2' => 'Lampro_C'});
my $Idir = '/media/Elements/introgression/';
my $Sdir = '/media/DataDrive/sutureGenomics/';
my $Slib = $Sdir . 'library';
my $np = 4;

foreach my $contact (keys %contacts) {
    my $target = bw($Idir,$contact);
    my ($lib1,$lib2) = readFiles($Slib,$Sdir,$contact,\%contacts);
    my $map1 = mapFiles($contact,$Sdir,$lib1,$target,$np,$Idir);
    my $map2 = mapFiles($contact,$Sdir,$lib2,$target,$np,$Idir);
    callVariants($contact,$Idir,$map1,$map2,$target);
}

sub callVariants {
    my ($contact,$Idir,$map1,$map2,$target) = @_;
    
    my $outVCF = $Idir . 'ancMapping/' . $contact . '.vcf';

    my $call1 = system("samtools faidx $target");
    my $reads1 = join("\t",@{$map1});
    my $reads2 = join("\t",@{$map2});
    my $reads = $reads1 . "\t" . $reads2;
    my $call2 = system("samtools mpileup -uf $target $reads | bcftools view -bvcg - > var.raw.bcf");
    my $call3 = system("bcftools view var.raw.bcf | vcfutils.pl varFilter -w 0 > $outVCF");
}

sub mapFiles {
    my ($contact,$indir,$lib,$target,$np,$outdir) = @_;
    
    $outdir .= 'ancMapping/';
    mkdir($outdir) unless (-d $outdir);
    $outdir .= $contact . '/';
    mkdir($outdir) unless(-d $outdir);

    my @map;

    foreach my $lineage (keys %{$lib}) {
	foreach my $lib (sort {$a cmp $b} keys %{$lib->{$lineage}}) {
	    my $file1 = $indir . $lineage . '/' . $lib . '/' . $lib . '_1p_final.fastq.gz';
	    my $file2 = $file1; $file2 =~ s/1p_/2p_/;
	    my $fileu = $file1; $fileu =~ s/1p_/u_/;

	    my $out1 = $lineage . '_' . $lib . '1';
	    my $out2 = $lineage . '_' . $lib . '2';
	    my $finalout = $outdir . $contact . '_' . $lib . '.sorted';

	    my $call1 = system("bowtie2 -5 5 -3 5 -p $np -x $target -U $file1,$file2,$fileu -S $out1\.sam --un-gz $out1\.gz");
	    my $call2 = system("bowtie2 -5 5 -3 5 -p $np -x $target -U $out1\.gz -S $out2\.sam --local");
	    my $call3 = system("samtools view -bS $out1\.sam > $out1\.bam");
	    my $call4 = system("samtools view -bS $out2\.sam > $out2\.bam");
	    my $call5 = system("samtools merge tmp.bam $out1\.bam $out2\.bam");
	    my $call6 = system("samtools sort tmp.bam $finalout");
	    my $call7 = system("rm $out1\.sam $out2\.sam $out1\.gz $out1\.bam $out2\.bam tmp.bam");
	    $finalout .= '.bam';
	    push(@map,$finalout);
	}
    }
    return(\@map);
}

    

sub readFiles {
    my ($lib,$dir,$contact,$c) = @_;
    my %lib1; my %lib2;
    my $lib1 = $c->{$contact}->{'1'};
    my $lib2 = $c->{$contact}->{'2'};
    open(IN, "<$lib");
    while(<IN>) {
	chomp(my $line = $_);
	my @d = split(/\t/,$line);
	if ($d[2] eq $lib1) {
	    $lib1{$lib1}{$d[0]}++;
	}
	if ($d[2] eq $lib2) {
	    $lib2{$lib2}{$d[0]}++;
	}
    }
    close(IN);
    return(\%lib1,\%lib2);
}

sub bw {
    my ($dir,$contact) = @_;
    my $file = $dir . 'targets/' . $contact . '_targets.fa';
    unless(-f $file . ".1.bt2") {
	my $call = system("bowtie2-build $file $file");
    }
    return($file);
}
