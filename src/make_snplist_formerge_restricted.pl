#!/usr/bin/perl

#  The code itself is Copyright (C) 2017, by David J. Cutler.
#  This library is free software; you can redistribute it and/or
#  modify it under the terms of the GNU Lesser General Public
#  License as published by the Free Software Foundation; either
#  version 2.1 of the License, or (at your option) any later version.
#  This library is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  Lesser General Public License for more details.
#  You should have received a copy of the GNU Lesser General Public
#  License along with this library; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

use IO::Zlib;
use strict vars;
use vars
  qw(@fields %covered_base %snp_count %good_count @good_list @bad_list %chr_num @sort_snp @ssort_snp $f $i $j @sample_names $dirname %chr %pos);

$dirname = ".";

if ( @ARGV != 2 ) {
  print "\n!!! Do not use. Use make_snplist_formerge.pl instead !!!\n";
  print "\n Usage: ${0} sdx_file outname \n\n ";
  exit(1);
}

open( FILE, "$ARGV[0]" )
  || die "\nCan't open file $ARGV[0] which should contain genome_sdx_file\n";
my $chr_count = <FILE> +0;
my $ii        = 0;
for ( my $i = 0; $i < $chr_count; $i++ ) {
  $_ = <FILE>;
  chomp;
  @fields = split('\t');
  $chr_num{"$fields[1]"} = $ii;
  # print "\n Chromosome $fields[1] has offset $cc \n";
  $ii++;
}
print "\n Finished reading genome \n";
close(FILE);
print "\n Reading the bed files \n";
opendir( DIR, $dirname ) || die "Cannot open directory $dirname";
my @files = sort grep { /\.bed$/ } readdir(DIR);
close(DIR);
my $bed_count = @files +0;
my $f_count   = 0;

foreach my $f (@files) {
  print "\n Working $f \n";
  open( FILE, "$f" ) || die "\n Can not open $f for reading \n";
  while (<FILE>) {
    chomp;
    @fields = split('\t');
    for ( my $i = $fields[1]; $i <= $fields[2]; $i++ ) {
      my $name = "$fields[0]\_$i";
      if ( $f_count == 0 ) {
        $covered_base{$name}++;
      }
      elsif ( exists( $covered_base{$name} ) ) {
        if ( $covered_base{$name} == $f_count ) {
          $covered_base{$name}++;
        }
      }
    }
  }
  $f_count++;
  close(FILE);
}
print "\n Finished reading the bedfiles \n";
foreach my $i ( keys %covered_base ) {
  if ( $covered_base{$i} < $bed_count ) {
    delete( $covered_base{$i} );
  }
}

opendir( DIR, $dirname ) || die "Cannot open directory $dirname";
@files = sort grep { /\.snp$/ } readdir(DIR);
close(DIR);
foreach my $f (@files) {
  open( FILE, "$f" ) || die "\n Can't open $f for reading \n";
  print "\n About to read $f \n";
  <FILE>;
  my $ii = 0;
  while (<FILE>) {
    @fields = split('\t');
    my $name = "$fields[0]\_$fields[1]";
    if ( exists( $covered_base{$name} ) ) {
      $chr{$name} = $fields[0];
      $pos{$name} = $fields[1];
      if ( ( $fields[5] ne "LOW" ) && ( $fields[5] ne "MESS" ) ) {
        $good_count{$name}++;
      }
      if ( exists( $snp_count{$name} ) ) {
        $snp_count{$name}++;
      }
      else {
        $snp_count{$name} = 1;
      }
    }
    $ii++;
    if ( $ii % 100000 == 0 ) {
      print "\n Read $ii lines \n";
    }
  }
  close(FILE);
}
print "\n About to sort \n";
foreach my $k ( keys %snp_count ) {
  if ( $good_count{$k} > 0 ) {
    push @good_list, $k;
  }
  else {
    push @bad_list, $k;
  }
}

@sort_snp = sort { compare_snp() } (@good_list);
print "\n About to print good list \n";
open( FILE, ">$ARGV[1].good.bed" )
  || die "\n Can not open $ARGV[1].good.bed for writing \n";
my $start = $sort_snp[0];
my $end   = $sort_snp[0];
shift(@sort_snp);
foreach my $i (@sort_snp) {
  if ( ( $chr{$i} eq $chr{$start} ) && ( $pos{$i} - $pos{$end} == 1 ) ) {
    $end = $i;
  }
  else {
    print FILE "$chr{$start}\t$pos{$start}\t$pos{$end}\n";
    $start = $i;
    $end   = $i;
  }
}
print FILE "$chr{$start}\t$pos{$start}\t$pos{$end}\n";
close(FILE);
@sort_snp = sort { compare_snp() } (@bad_list);
print "\n About to print bad list\n";
open( FILE, ">$ARGV[1].bad.bed" )
  || die "\n Can not open $ARGV[1].bad.bed for writing \n";
$start = $sort_snp[0];
$end   = $sort_snp[0];
shift(@sort_snp);

foreach my $i (@sort_snp) {
  if ( ( $chr{$i} eq $chr{$start} ) && ( $pos{$i} - $pos{$end} == 1 ) ) {
    $end = $i;
  }
  else {
    print FILE "$chr{$start}\t$pos{$start}\t$pos{$end}\n";
    $start = $i;
    $end   = $i;
  }
}
print FILE "$chr{$start}\t$pos{$start}\t$pos{$end}\n";
close(FILE);

sub compare_snp ($a $b) {
  if ( $chr_num{ $chr{$a} } < $chr_num{ $chr{$b} } ) {
    return -1;
  }
  if ( $chr_num{ $chr{$b} } < $chr_num{ $chr{$a} } ) {
    return 1;
  }
  if ( $pos{$a} < $pos{$b} ) {
    return -1;
  }
  if ( $pos{$b} < $pos{$a} ) {
    return 1;
  }
  return 0;
}

