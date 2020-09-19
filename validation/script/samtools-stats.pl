#!/usr/bin/env perl
#
### satmtools-stats.pl
##################################################################################
#  Simple/crude perl script to find BAMs and execute samtools stat for BAMs using either sbatch or qsub
#  Modify and update this script for more general usage or pipeline integration
#
### HISTORY
#######################################################################################
# Version               Date                Coder                     Comments
# 0.01                  2019-09-18         Takafumi Yamaguchi        Initial development.
#
#
###INCLUDES
######################################################################################
use warnings;
use strict;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use File::Path;
use File::Basename;
use Carp;
### COMMAND LINE DEFAULT ARGUMENTS ################################################################
# list of arguments and default values go here as hash key/value pairs
our %opts = (
	path_nf    => '/hot/pipelines/development/slurm/align-DNA/outputs/temp-caden',
	path_hpci  => '/data/users/shutao/sge_germline_variants/mapping/*/*/output_bwa/*_DNA_alignments/bwa/0.7.15/*/',
	cluster => 'Slurm'
	);

### MAIN CALLER ###################################################################################
my $result = main();
exit($result);

### FUNCTIONS #####################################################################################

### main ##########################################################################################
# Description:
#   Main subroutine for program
# Input Variables:
#   %opts = command line arguments
# Output Variables:
#   N/A

sub main {
	
	# get the command line arguments
	GetOptions(
		\%opts,
		"help|?",
		"man",
		"path_nf=s"            => \$opts{'path_nf'},
		"path_hpci=s"          => \$opts{'path_hpci'}
	) or pod2usage(64);
	
    if ($opts{'help'}) { pod2usage(1) };
    if ($opts{'man'}) { pod2usage(-exitstatus => 0, -verbose => 2) };
	while(my ($arg, $value) = each(%opts)){
		if (!($arg=~/\:/) and !defined $value) {
		print "ERROR: Missing argument $arg\n";
		pod2usage(2);
		}
	  }

	my @bams = split /\n/, `find $opts{path_nf} -name "*bam"`;

	foreach my $bam (sort @bams){
		#print $bam."\n";
		
		my $job='samtools-stats-'.basename($bam);
		my ($size, $pipeline) = ('', 'v3-update-merge');

		if($bam=~/full/){
			
			$size = 'full';
			
			} elsif($bam=~/partial|mini/){
			
			$size = 'partial';
			
			} else {
			
			$size = 'full';
			
			}
		
		(my $output = basename($bam)) =~s/bam/stats/;
		
		$output = join('-', $size, $pipeline, $output);
		
		my $cmd = join(' ', 'samtools', 'stats --threads 4', $bam, '>', '/hot/pipelines/development/slurm/align-DNA/validation/'.$output);
		
		# submit samtool stat job using sbatch or qsub (using 4 cores)
		if($opts{cluster} eq 'Slurm'){
			
			my $sbatch = join(' ', 'sbatch -p midmem -J', $job, '-o', $job.'.err', '--mem=5G -c 4 --wrap='.'"'.$cmd.'"');
			#print $sbatch."\n";
			#system($sbatch);

		} else {

			my $qsub = join(' ', 'qsub -cwd -b y -m aes -N', $job, '-l h_vmem=2g,slot_type=midmem -pe smpslots 4 ', '"'.$cmd.'"');
			#print $qsub."\n";
			#system($qsub);
			#exit;

		}

		# simple grep command to extract core bam stats 
		(my $output2 = basename($output)) =~s/stats/core-stats/;

		my $grep_cmd = join(' ', 'grep SN /hot/pipelines/development/slurm/align-DNA/validation/'.$output, '| cut -f 2-3 | grep -v -P "#"', '>', '/hot/pipelines/development/slurm/align-DNA/validation/'.$output2);
		#print $grep_cmd."\n";
		#system($grep_cmd);

		#exit;
	}
}
