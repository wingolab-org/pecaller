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
use vars
  qw(@fields %tran_map @filehandles %ts %tv %transition %transver @files $i $j $last @samples);

if ( @ARGV != 3 ) {
  die "\n Usage: ${0} SNPFILE threshold_to_call call_rate\n";
}

open( FILE, "$ARGV[0]" )
  || die "\n Could not open $ARGV[0] which should be the snp file \n";
$_ = <FILE>;
print;
chomp;
@fields = split('\t');
my $i, $j;
$i = 0;

for ( $j = 6; $j < @fields; $j += 2, $i++ ) {
  $samples[$i] = $fields[$j];
}
my $tot_samples = $i;
while (<FILE>) {
  my $line = $_;
  chomp;
  @fields = split('\t');
  my $this_calls = 0;
  my $diff_ref   = 0;
  for ( $j = 6; $j < @fields; $j += 2 ) {
    if ( ( $fields[ $j + 1 ] >= $ARGV[1] ) && ( $fields[$j] ne "N" ) ) {
      $this_calls++;
      if ( $fields[$j] ne $fields[2] ) {
        $diff_ref++;
      }
    }
  }
  if ( ( $diff_ref > 0 ) && ( $this_calls / $tot_samples >= $ARGV[2] ) ) {
    print $line;
  }
}
