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
  qw(@fields $site_code $snp_id %tran_map %snp_type @filehandles %rs %rs_type %ts %tv %transition %transver @files $i $j $last @samples);

if ( @ARGV != 3 ) {
  die "\n Usage: ${0} SNPFILE ANNOTATION_FILE TYPE_OF_VARIANT[SNP,LOW,MESS,DEL,INS]\n";
}

open( FILE, "$ARGV[1]" )
  || die "\n Could not open $ARGV[1] which should be the annotation file \n";
my $type = uc( $ARGV[2] );
#print "\nSearching for $type \n";
$tran_map{"AG"} = 1;
$tran_map{"GA"} = 1;
$tran_map{"CT"} = 1;
$tran_map{"TC"} = 1;

$_      = <FILE>;
@fields = split('\t');
for ( $i = 0; $i < @fields; $i++ ) {
  if ( $fields[$i] eq "annotation_type" ) {
    $site_code = $i;
  }
  elsif ( $fields[$i] eq "snp_id" ) {
    $snp_id = $i;
  }
}
while (<FILE>) {
  chomp;
  #my $line = $_;
  @fields = split('\t');
  my $temp = "$fields[0]\_$fields[1]";
  my @sfields = split( '\;', $fields[$site_code] );
  my %sc_hash;
  foreach my $i (@sfields) {
    $sc_hash{$i} = 1;
  }
  @sfields = ( sort ( keys %sc_hash ) );
  my $sc = $sfields[0];
  for ( $i = 1; $i < @sfields; $i++ ) {
    $sc .= ";$sfields[$i]";
  }
  $snp_type{$temp} = $sc;
  #$_ = $fields[$site_code];
  #if(/NA/)
  #{
  #	print "Everything messed up.  line = $line\n site_code is $site_code \n which is $fields[$site_code]\n";
  #	exit();
  #}
  $_ = $fields[$snp_id];
  if (/^rs/) {
    $rs_type{$temp} = 1;
  }
  else {
    $rs_type{$temp} = 0;
  }

}
close(FILE);
open( FILE, "$ARGV[0]" )
  || die "\n Could not open $ARGV[0] which should be the snp file \n";
$_ = <FILE>;
chomp;
@fields = split('\t');
my $i, $j;
$samples[0] = "ALL";
$i = 1;

for ( $j = 6; $j < @fields; $j += 2, $i++ ) {
  $samples[$i] = $fields[$j];
}

while (<FILE>) {
  chomp;
  @fields = split('\t');
  $_      = uc( $fields[5] );
  #print "\n This one is a $_ \n";
  if (/$type/) {
    #print "\n Doing stuff \n";
    my $spot = "$fields[0]\_$fields[1]";
    if ( exists( $snp_type{$spot} ) ) {
      #$_ = $snp_type{$spot};
      #if(/NA/)
      #{
      #	print "\n Something is fucked up with $spot which has type $snp_type{$spot} and rs $rs_type{$spot} \n";
      #}
      my $this_type = $snp_type{$spot};
      $_ = $fields[2];
      my $is_trans = 0;
      my $is_rs    = $rs_type{$spot};
      $rs{$this_type}[0] += $is_rs;
      $i = 0;
      if ( ( $fields[3] eq "A,G" )
        || ( $fields[3] eq "C,T" )
        || ( $fields[3] eq "G,A" )
        || ( $fields[3] eq "T,C" ) )
      {
        $is_trans = 1;
      }
      elsif ( exists( $tran_map{"$fields[2]$fields[3]"} ) ) {
        $is_trans = 1;
      }
      if ($is_trans) {
        $ts{$this_type}[0]++;
      }
      else {
        $tv{$this_type}[0]++;
      }
      $i = 1;
      for ( $j = 6; $j < @fields; $j += 2, $i++ ) {
        if ( $fields[$j] ne $fields[2] && $fields[$j] ne "N" ) {
          $rs{$this_type}[$i] += $is_rs;
          if ($is_trans) {
            $ts{$this_type}[$i]++;
          }
          else {
            $tv{$this_type}[$i]++;
          }
        }
      }
    }
  }
}
my @tot_types = ( sort ( keys %tv ) );
print "\nSample";
foreach my $i (@tot_types) {
  print
    "\t$i\_transistions\t$i\_transversions\t$i\_ratio\t$i\_in_dbsnp\t$i\_dbsnp_ratio";
}
for ( $i = 0; $i < @samples; $i++ ) {
  print "\n$samples[$i]";
  $j = 0;
  foreach my $ii (@tot_types) {
    my $tot      = $ts{$ii}[$i] + $tv{$ii}[$i];
    my $rs_ratio = 0;
    if ( $tot > 0 ) {
      $rs_ratio = $rs{$ii}[$i] / $tot;
    }
    my $j = 1.0;
    if ( $tv{$ii}[$i] > 0 ) {
      $j = $ts{$ii}[$i] / $tv{$ii}[$i];
    }
    print "\t$ts{$ii}[$i]\t$tv{$ii}[$i]\t$j\t$rs{$ii}[$i]\t$rs_ratio";
  }
}
print "\n";
