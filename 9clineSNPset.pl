use warnings;
use strict;

my @contacts = qw(gillies nBMB);
my @lib = qw(10kN 2kN nTail 1kN center 1kS sTail 2kS 10kS);
my $Idir = '/media/Elements/introgression/';
my $np = 4;

foreach my $contact ( @contacts) {
    my $target = bw($Idir,$contact);
    my $map = mapFiles($contact,$Idir,\@lib,$target,$np);
    callVariants($contact,$Idir,$map,$target);
}

sub callVariants {
    my ($contact,$Idir,$map,$target) = @_;
    
    my $outVCF = $Idir . 'ancMapping/' . $contact . '_HZ.vcf';

    my $call1 = system("samtools faidx $target");
    my $reads = join("\t",@{$map});
    my $call2 = system("samtools mpileup -B -q 20 -uf $target $reads | bcftools view -bvcg - > $contact\.HZraw.bcf");
    my $call3 = system("bcftools view -I $contact\.HZraw.bcf > $outVCF");
    my $call4 = system("rm $contact\.HZraw.bcf");
}

sub mapFiles {
    my ($contact,$dir,$library,$target,$np) = @_;
    
    my $outdir = $dir . 'ancMapping/';
    mkdir($outdir) unless (-d $outdir);
    $outdir .= $contact . '/';
    mkdir($outdir) unless(-d $outdir);

    my @map; 

    foreach my $lib (@$library) {
	my $final = $outdir . $lib . '.sorted.bam';
	unless(-f $final) {
	    my $file1 = $dir . $contact . '/' . $lib . '/' . $lib . '_1_final.fastq.gz';
	    my $file2 = $file1; $file2 =~ s/1_/2_/;
	    my $fileu = $file1; $fileu =~ s/1_/u_/;

	    my $out1 = $contact . '_' . $lib . '1';
	    my $out2 = $contact . '_' . $lib . '2';
	    my $finalout = $outdir . $lib . '.sorted';

	    my $call1 = system("bowtie2 -5 5 -3 5 -p $np -x $target -U $file1,$file2,$fileu -S $out1\.sam --un-gz $out1\.gz");
	    my $call2 = system("bowtie2 -5 5 -3 5 -p $np -x $target -U $out1\.gz -S $out2\.sam --local");
	    my $call3 = system("samtools view -bS $out1\.sam > $out1\.bam");
	    my $call4 = system("samtools view -bS $out2\.sam > $out2\.bam");
	    my $call5 = system("samtools merge $lib\.bam $out1\.bam $out2\.bam");
	    my $call6 = system("samtools sort $lib\.bam $finalout");
	    my $call7 = system("rm $out1\.sam $out2\.sam $out1\.gz $out1\.bam $out2\.bam $lib\.bam");
	}
	push(@map,$final);
    }
    return(\@map);
}

sub bw {
    my ($dir,$contact) = @_;
    my $file = $dir . 'targetSequences/final/' . $contact . '_targets.fa.final';
    unless(-f $file . ".1.bt2") {
	my $call = system("bowtie2-build $file $file");
    }
    return($file);
}
