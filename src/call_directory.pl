#!/usr/bin/perl
# Run all chips
use strict 'vars';
use POSIX;
use vars
  qw($dirname $tail %matches @data @files %map_id %count %mapper_name @filess $m $name @fields @line $i $j $k $step_size $snps);

$dirname = ".";

my $caller = "/home/dcutler/Software/bin/pecaller";

if ( @ARGV < 2 || @ARGV > 4 ) {
  print "\n Usage: ${0} directory genome.sdx [guide_file_bed_format] [ped file]\n";
  exit(1);
}

chdir("$ARGV[0]") || die "\nCan't Change Directory to $ARGV[0] \n";
opendir( DIR, $dirname ) || die "Cannot open directory";
@files = grep { /pileup.gz/ } readdir(DIR);
close(DIR);
my $tot_files = @files + 5;

open( FILE, ">$ARGV[0]_call.sh" );

if ( @ARGV == 2 ) {
  print FILE "\n\n$caller pileup $ARGV[1] $tot_files $ARGV[0] 0.95 0.001 n 24 n \n\n";
}
elsif ( @ARGV == 3 ) {
  print FILE
    "\n\n$caller pileup $ARGV[1] $tot_files $ARGV[0] 0.95 0.001 n 24 n $ARGV[2]\n\n";
}
else {
  print FILE
    "\n\n$caller pileup $ARGV[1] $tot_files $ARGV[0] 0.95 0.001 n 24 y $ARGV[3] 1e-8 $ARGV[2]\n\n";
}

close(FILE);
system("qsub  -v USER -v PATH -cwd -q pe.q -o /dev/null -j y $ARGV[0]_call.sh");

