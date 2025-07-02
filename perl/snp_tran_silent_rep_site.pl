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
  qw(@fields %a_count %by_freq_ts %by_freq_tv %by_freq_rs $site_code $snp_id %tran_map %snp_type @filehandles %rs_ts %rs_tv %by_freq_rs_ts %by_freq_rs_tv %rs %rs_type %ts %tv %transition %transver @files $i $j $last @samples);

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
$a_count{"A"}   = 2;
$a_count{"C"}   = 2;
$a_count{"G"}   = 2;
$a_count{"T"}   = 2;
$a_count{"I"}   = 2;
$a_count{"D"}   = 2;

$a_count{"R"} = 1;
$a_count{"Y"} = 1;
$a_count{"S"} = 1;
$a_count{"W"} = 1;
$a_count{"K"} = 1;
$a_count{"M"} = 1;
$a_count{"H"} = 1;
$a_count{"E"} = 1;

$_      = <FILE>;
@fields = split('\t');
my $max_freq = 0;
for ( $i = 0; $i < @fields; $i++ ) {
  if ( $fields[$i] eq "refSeq.exonicAlleleFunction" ) {
    $site_code = $i;
  }
  elsif ( $fields[$i] eq "dbSNP146.name" ) {
    $snp_id = $i;
  }
  elsif ( $fields[$i] eq "dbSNP.name" ) {
    $snp_id = $i;
  }
  elsif ( $fields[$i] eq "dbSNP" ) {
    $snp_id = $i;
  }
}
if ( $site_code < 1 ) {
  die "\n Could not find refSeq.exonicAlleleFunction \n";
}
if ( $snp_id < 1 ) {
  die "\n Could not find dbSNP146.name \n";
}
while (<FILE>) {
  chomp;
  #my $line = $_;
  @fields = split('\t');
  my $temp    = "$fields[0]\_$fields[1]";
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
        if ($is_rs) {
          $rs_ts{$this_type}[0]++;
        }
      }
      else {
        $tv{$this_type}[0]++;
        if ($is_rs) {
          $rs_tv{$this_type}[0]++;
        }
      }
      $i = 1;
      my $this_count = 0;
      for ( $j = 6; $j < @fields; $j += 2, $i++ ) {
        if ( ( $fields[$j] ne $fields[2] )
          && ( $fields[$j] ne "N" )
          && ( $fields[ $j + 1 ] >= 0.95 ) )
        {
          $this_count += $a_count{ $fields[$j] };
          $rs{$this_type}[$i] += $is_rs;
          if ($is_trans) {
            $ts{$this_type}[$i]++;
            if ($is_rs) {
              $rs_ts{$this_type}[$i]++;
            }
          }
          else {
            $tv{$this_type}[$i]++;
            if ($is_rs) {
              $rs_tv{$this_type}[$i]++;
            }
          }
        }
      }
      if ( $this_count > $max_freq ) {
        $max_freq = $this_count;
      }
      $by_freq_rs{$this_type}[$this_count] += $is_rs;
      if ($is_trans) {
        $by_freq_ts{$this_type}[$this_count]++;
        if ($is_rs) {
          $by_freq_rs_ts{$this_type}[$this_count]++;
        }
      }
      else {
        $by_freq_tv{$this_type}[$this_count]++;
        if ($is_rs) {
          $by_freq_rs_tv{$this_type}[$this_count]++;
        }
      }
    }
  }
}
my @tot_types = ( sort ( keys %tv ) );
print "Sample";
foreach my $i (@tot_types) {
  print
    "\t$i\_transistions\t$i\_transversions\t$i\_ratio\t$i\_in_dbsnp\t$i\_dbsnp_ratio\t$i\_dbsnp_ts\t$i\_dbsnp_tv\t$i\_dsnp_ratio";
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
    my $jj = 0;
    if ( $rs_tv{$ii}[$i] > 0 ) {
      $jj = $rs_ts{$ii}[$i] / $rs_tv{$ii}[$i];
    }
    print
      "\t$ts{$ii}[$i]\t$tv{$ii}[$i]\t$j\t$rs{$ii}[$i]\t$rs_ratio\t$rs_ts{$ii}[$i]\t$rs_tv{$ii}[$i]\t$jj";
  }
}
for ( $i = 1; $i <= $max_freq; $i++ ) {
  print "\nNonRef_Allele_Count_$i";
  $j = 0;
  foreach my $ii (@tot_types) {
    my $tot      = $by_freq_ts{$ii}[$i] + $by_freq_tv{$ii}[$i];
    my $rs_ratio = 0;
    if ( $tot > 0 ) {
      $rs_ratio = $by_freq_rs{$ii}[$i] / $tot;
    }
    my $j = 1.0;
    if ( $by_freq_tv{$ii}[$i] > 0 ) {
      $j = $by_freq_ts{$ii}[$i] / $by_freq_tv{$ii}[$i];
    }
    my $jj = 1.0;
    if ( $by_freq_rs_tv{$ii}[$i] > 0 ) {
      $jj = $by_freq_rs_ts{$ii}[$i] / $by_freq_rs_tv{$ii}[$i];
    }
    print
      "\t$by_freq_ts{$ii}[$i]\t$by_freq_tv{$ii}[$i]\t$j\t$by_freq_rs{$ii}[$i]\t$rs_ratio\t$by_freq_rs_ts{$ii}[$i]\t$by_freq_rs_tv{$ii}[$i]\t$jj";
  }
}
print "\n";
