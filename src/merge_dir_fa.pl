#!/usr/bin/env perl
# Name:           merge_dir_fa.pl
# Date Created:   Tue Sep  1 11:47:42 2015
# Date Modified:  Tue Sep  1 11:47:42 2015
# By:             TS Wingo
#
# Description:    takes a directory with *.fa.gz files; reads those files
#                 and writes a merged fa file in the natural order of the
#                 chromosomes with the unplaced chromosomes at the end;
#                 use `grep '>' out_ext.*.fa` to check the order of the
#                 chromosomes;

use 5.10.0;
use strict;
use warnings;
use Getopt::Long;
use IO::Uncompress::Gunzip qw($GunzipError);
use Path::Tiny;
use Time::localtime;
use DDP;

# objects
my $now_timestamp = sprintf( "%d-%02d-%02d",
  ( localtime->year() + 1900 ),
  ( localtime->mon() + 1 ),
  localtime->mday() );

# variables
my ( @chrs, %chrFasta, %chrFastaCount, %chrFastaPrinted );
my ( $dir_name, $out_ext, $chr_list );

# get options
die
"Usage: $0 
  -d <dir with fasta files> 
  -c <chr_list (e.g., 1-19,M,X,Y)> 
  -o <base name for output fa file>\n"
  unless GetOptions(
  'c|chr_list=s' => \$chr_list,
  'd|dir=s'      => \$dir_name,
  'o|out=s'      => \$out_ext,
  ) and $chr_list
    and $dir_name
    and $out_ext;

my $chrs_aref = ProcChromList($chr_list);
PrintVaildChrList( $chrs_aref );

# read directory with fa files
my $dir = Path::Tiny->new($dir_name);
for my $file ( $dir->children ) {
  next unless $file->basename =~ m/\.fa\.gz\z/;
  ( my $name = $file->basename ) =~ s/\.fa\.gz\z//;
  my $fhz = new IO::Uncompress::Gunzip $file->absolute->stringify or die "$GunzipError: $!\n";
  local $/;
  $chrFasta{$name}      = <$fhz>;
  $chrFastaCount{$name} = length $chrFasta{$name};
}
PrintFaFiles( \%chrFasta );

# open file for writing
my $fa_fh = IO::File->new( "$out_ext.$now_timestamp.fa", 'w' ) || die "$!\n";

# print named chromosomes
for my $chr (@$chrs_aref) {
  if ( exists $chrFasta{$chr} ) {
    print {$fa_fh} $chrFasta{$chr};
    $chrFastaPrinted{$chr}++;
  }
  else {
    die "ERROR: Did not find expected chr '$chr'\n";
  }
}

# print unplaced chromosomes
for my $chr ( sort keys %chrFasta ) {
  next if exists $chrFastaPrinted{$chr};
  print {$fa_fh} $chrFasta{$chr};
}

%chrFastaPrinted = ();

# print out what we wrote
say "Wrote the following to $out_ext.$now_timestamp.fa:";
for my $chr ( @$chrs_aref, sort keys %chrFastaCount ) {
  next if exists $chrFastaPrinted{$chr};
  say join "\t", $chr, $chrFastaCount{$chr};
  $chrFastaPrinted{$chr}++;
}

sub PrintFaFiles {
  my $href = shift;
  say "Found files";
  say "\t" . $_ for sort keys %$href;
}

sub PrintVaildChrList {
  my $aref = shift;
  say "Chromomsome List:";
  say "\t" . $_ for @$aref;
  say "Will add unplaced and other chromsomes to the end of this list.";
}

sub ProcChromList {
  my $list = shift;
  my @chrs;
  my @chr_entries = split /\,/, $list;
  for my $chr (@chr_entries) {
    if ( $chr =~ m/(\d+)\-(\d+)/ ) {
      if ($1 > $2) {
        say "ERROR: expected range to be ascending, but got '$1'-'$2' from '$list'";
        exit(1);
      }
      for ( my $i = $1 ; $i <= $2 ; $i++ ) {
        push @chrs, "chr$i";
      }
    }
    elsif ( $chr =~ m/\A[\d|M|X|Y]+\z/  ) {
      push @chrs, "chr$chr";
    }
    else {
      say "unrecognized chr: $chr in expression: $list";
    }
  }
  if (@chrs) {
    return \@chrs;
  }
  else {
    say "failed to process list: '$list'";
    exit(1);
  }
}
