#!/usr/bin/perl

# The code itself is Copyright (C) 2016, by David J. Cutler.
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

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

