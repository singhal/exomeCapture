##############################################################################################################################
# a script to prepare shell scripts to use for assembly with multiple de novo assemblers on TACC                             #
# external dependencies: to run the script, none but ... to use the scripts: abyss, velvet, trinity, oases, soapdenovo-trans #
# written by Sonal Singhal, sonal.singhal1 [at] gmail.com, 16 Dec 2012                                                       #
##############################################################################################################################

use warnings;
use strict;

#home directory with all my files; make sure to end with a backslash
my $scratch = '/scratch/01302/ssinghal/';
my $home = '/work/01302/ssinghal/';
my $email = 'sonal.singhal1@gmail.com';
my $insert = 30; #how long is the expected distance between the two reads

#parameters for soap
my @soapkmer = qw(21 31 41 51 61 71 81 91); #max 96
my $soapdir = '/home/01302/ssinghal/bin/';

#parameters for abyss
my @abysskmer = qw(21 31 41 51 61 71 81 91); #max 96

#parameters for velvet
my @velvetkmer = qw(31 41 51 61 71 81 91); #max 96
my $velvetdir = '/home/01302/ssinghal/bin/'; #directory where velvet executables are

my @lib = qw(carlia sjo nBMB gillies);

###########################
# run the subroutines     #
###########################

foreach my $lib (@lib) {
	my $file_1 = $scratch . $lib . '_1.fastq';
	my $file_2 = $scratch . $lib . '_2.fastq';
	my $file_u = $scratch . $lib . '_u.fastq';
	my $file_1s = $scratch . $lib . '_1.fastq';
	my $file_2s = $scratch . $lib . '_2.fastq';
	my $file_us = $scratch . $lib . '_u.fastq';
	if (-f $file_1) {
		makeABYSS($file_1,$file_2,$file_u,$lib);
		}
	if (-f $file_us) {
		makeSOAP($file_1s,$file_2s,$file_us,$lib);
		}
	if (-f $file_us) {		
		makeVelvet($file_1s,$file_2s,$file_us,$lib);
		}
	}

###########################
# behold the subroutines! #
###########################

sub makeSOAP {
	my ($file_1,$file_2,$file_u,$lib) = @_;

	my $resultsDir =  $scratch . 'soapResults';
	my $runfileDir =  $scratch . 'soapScripts';
	mkdir($resultsDir) unless (-d $resultsDir);
	mkdir($runfileDir) unless (-d $runfileDir);
	foreach my $k (@soapkmer) {
		my $job = 'soap' . "_" . $lib . "_k" . $k;
		open(CFG, ">$runfileDir" . "/" . $job . ".config");
		open(OUT, ">$runfileDir" . "/" . $job . ".sh");

		print OUT "#!/bin/bash\n"; 
		print OUT "#\$ -N $job\n"; 
		print OUT "#\$ -j y\n"; 
		print OUT "#\$ -P hpc\n"; 
		print OUT "#\$ -o $resultsDir" , "/", $job, ".out\n"; 
		print OUT "#\$ -pe 1way 16\n"; 
		print OUT "#\$ -q largemem\n"; 
		print OUT "#\$ -l h_rt=08:00:00\n"; 
		print OUT "#\$ -M $email\n"; 
		print OUT "#\$ -m e\n";
		print OUT "#\$ -cwd\n"; 
		print OUT "#\$ -V\n";
		print OUT $soapdir . "SOAPdenovo-127mer all -s $runfileDir" . "/" . $job . ".config -K $k -o $resultsDir" . "/" . "$job -p 16\n";		
			
		print CFG "max_rd_len=200\n";	
		print CFG "\[LIB\]\n";
		print CFG "avg_ins=$insert\n";
		print CFG "asm_flags=3\n";
		print CFG "rank=1\n";
		print CFG "q1=$file_1\n";
		print CFG "q2=$file_2\n";
		print CFG "q=$file_u\n";

		close(CFG); close(OUT);
		}
    }

sub makeVelvet {
	my ($file_1,$file_2,$file_u,$lib) = @_;
	
	my $runfileDir = $scratch . 'velvetScripts/';
	mkdir($runfileDir) unless (-d $runfileDir);
	my $allresults = $scratch . 'velvetResults/';
	mkdir($allresults) unless (-d $allresults);

	foreach my $k (@velvetkmer) {
		my $jobH = 'velveth' . "_" . $lib . "_k" . $k; 
		my $jobG = 'velvetg' . "_" . $lib . "_k" . $k; 

		my $resultsDir = $allresults . 'velvet_' . $lib . "_k" . $k; 
		mkdir($resultsDir) unless (-d $resultsDir);
		
		open(OUTH, ">$runfileDir" . $jobH . ".sh");
		open(OUTG, ">$runfileDir" . $jobG . ".sh");

		print OUTH "#!/bin/bash\n"; 
		print OUTH "#\$ -N $jobH\n"; 
		print OUTH "#\$ -j y\n"; 
		print OUTH "#\$ -P hpc\n"; 
		print OUTH "#\$ -o $resultsDir" , "/", $jobH, ".out\n"; 
		print OUTH "#\$ -pe 1way 8\n"; 
		print OUTH "#\$ -q largemem\n"; 
		print OUTH "#\$ -l h_rt=08:00:00\n"; 
		print OUTH "#\$ -M $email\n"; 
		print OUTH "#\$ -m e\n";
		print OUTH "#\$ -cwd\n"; 
		print OUTH "#\$ -V\n";
		print OUTH "$velvetdir" . 'velveth ' . $resultsDir . " $k " . '-shortPaired -separate -fastq ' . $file_1 . ' ' . $file_2 . ' -short2 -fastq ' . $file_u . "\n";

		print OUTG "#!/bin/bash\n"; 
		print OUTG "#\$ -N $jobG\n"; 
		print OUTG "#\$ -j y\n"; 
		print OUTG "#\$ -P hpc\n"; 
		print OUTG "#\$ -o $resultsDir" , "/", $jobG, ".out\n"; 
		print OUTG "#\$ -pe 1way 8\n"; 
		print OUTG "#\$ -q largemem \n"; 
		print OUTG "#\$ -l h_rt=08:00:00\n"; 
		print OUTG "#\$ -M $email\n"; 
		print OUTG "#\$ -m e\n";
		print OUTG "#\$ -cwd\n"; 
		print OUTG "#\$ -V\n";
		print OUTG "$velvetdir" . 'velvetg ' . $resultsDir . ' -exp_cov auto -cov_cutoff auto -ins_length ' . $insert .  "\n";

		close(OUTG); close(OUTH); 
		}
	}
	
sub makeABYSS {
	my ($file_1,$file_2,$file_u,$lib) = @_;
	my $resultsDir =  $home . 'abyssResults';
	my $runfileDir =  $home . 'abyssScripts';
	mkdir($resultsDir) unless (-d $resultsDir);
	mkdir($runfileDir) unless (-d $runfileDir);
	foreach my $k (@abysskmer) {			
	    my $job = 'abyss' . "_" . $lib . "_k" . $k;
	    my $dir = $resultsDir . '/' . $job;
	    mkdir($dir) unless (-d $dir);
	    open(OUT, ">/$runfileDir" . "/" . $job . ".sh");
	    print OUT "#!/bin/bash\n"; 
	    print OUT "#\$ -N $job\n"; 
	    print OUT "#\$ -j y\n"; 
	    print OUT "#\$ -o $dir" , "/", $job, ".out\n"; 
	    print OUT "#\$ -pe 16way 128\n"; 
	    print OUT "#\$ -q normal\n"; 
	    print OUT "#\$ -l h_rt=24:00:00\n"; 
	    print OUT "#\$ -M $email\n"; 
	    print OUT "#\$ -m e\n";
	    print OUT "#\$ -cwd\n"; 
	    print OUT "#\$ -V\n";
	    print OUT "cd $dir\n";
	    print OUT "abyss-pe k=$k mpirun=/opt/apps/pgi7_2/openmpi/1.3/bin/mpirun np=32 n=5 s=200 in=\'$file_1 $file_2\' se=$file_u name=$dir/$job\n";
		}
	}

