#!/usr/bin/perl

# The code itself is Copyright (C) 2017, by David J. Cutler.
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

if ( @ARGV != 2 ) {
  print "\n Usage: ${0} directory genome.sdx \n";
  exit(1);
}

chdir("$ARGV[0]") || die "\nCan't Change Directory to $ARGV[0] \n";

opendir( DIR, "." ) || die "Cannot open directory";
@files = grep { !/mfile$/ } grep { /fastq/ } readdir(DIR);
close(DIR);
my $tot_files = 0;
foreach my $f (@files) {
  @fields = split( '\.', $f );
  $count{ $fields[0] }++;
  $tail = "";
  for ( $i = 1; $i < @fields; $i++ ) {
    $tail .= ".$fields[$i]";
  }
  $tot_files++;
}
@filess = sort ( keys %count );
open( FILE1, ">$ARGV[0].file1.txt" )
  || die "\n Can't open $ARGV[0].file1.txt for writing \n";
open( FILE2, ">$ARGV[0].file2.txt" )
  || die "\n Can't open $ARGV[0].file2.txt for writing \n";
for ( $i = 0; $i < @filess; $i++ ) {
  $_ = $filess[$i];
  s/_1_/_2_/;
  s/_R1_/_R2_/;
  my $s1 = $_;
  for ( $j = $i + 1; $j < @filess; $j++ ) {
    $_ = $filess[$j];
    s/_1_/_2_/g;
    s/_R1_/_R2_/g;
    my $s2 = $_;
    #print "\n $filess[$i] $s1 $files[$j] $s2 \n";
    if ( $filess[$i] eq $s2 || $filess[$j] eq $s1 ) {
      print("\n Matching $filess[$i] with $filess[$j] \n");
      $matches{ $filess[$i] } = $filess[$j];
      $matches{ $filess[$j] } = $filess[$i];
      $tot_files--;
    }
  }
}
my %done;
my $last_i    = 0;
my $pair_type = "sa";

for ( $i = 0; $i < @filess; $i++ ) {
  if ( !$done{ $filess[$i] } ) {
    if ( !$matches{ $filess[$i] } ) {
      $done{ $filess[$i] } = 1;
      print FILE1 "$filess[$i]$tail\n";
      $last_i = $i;
    }
    else {
      $pair_type                       = "pa";
      $done{ $filess[$i] }             = 1;
      $done{ $matches{ $filess[$i] } } = 1;
      print FILE1 "$filess[$i]$tail\n";
      print FILE2 "$matches{$filess[$i]}$tail\n";
      $last_i = $i;
    }
  }
}
close(FILE1);
close(FILE2);
open( FILE, ">$ARGV[0].sh" );
print FILE
  "\n\n/home/dcutler/Software/bin/pemapper  $ARGV[0] $ARGV[1] $pair_type $ARGV[0].file1.txt ";
if ( $pair_type eq "pa" ) {
  print FILE "$ARGV[0].file2.txt 500 0 ";
}
else {
  system("rm $ARGV[0].file2.txt");
}
print FILE "N 0.85 24 200000000\n\n";
close(FILE);
system("qsub  -v USER -v PATH -cwd -q pe.q -o /dev/null -j y $ARGV[0].sh");

