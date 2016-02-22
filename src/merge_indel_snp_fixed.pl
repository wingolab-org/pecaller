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

use IO::Zlib;
use strict vars;
use vars qw(@fields %chr_num @sort_snp @ssort_snp $f $i $j %ins_a_count %chr_offset %ins_consensus %ins_needed @sample_names $dirname @lines @chr @pos @type);

$dirname = ".";

if(@ARGV != 4) 
{
	print "\n Usage: ${0} sdx_file snp_file directory_indel_files outname \n\n "; 
	exit(1);
}

open(FILE,"$ARGV[0]") || die "\nCan't open file $ARGV[0] which should contain genome_sdx_file\n";
my $chr_count = <FILE> + 0;
my $cc = 0;
my $ii = 0;
for(my $i = 0;$i<$chr_count;$i++)
{
	$_ = <FILE>;
    	chomp;
	@fields = split('\t');
	$chr_offset{"$fields[1]"} = $cc;
	$chr_num{"$fields[1]"} = $ii;
	# print "\n Chromosome $fields[1] has offset $cc \n";
	$ii++;
	# $cc+=15;
}
close(FILE);
print "\n Finished Reading Genome File \n";
open(FILE,"$ARGV[1]") || die "\nCan't open file $ARGV[1] which should contain sequencing data\n";
$_ = <FILE>;
chomp;
my $header = $_;
@fields = split("\t");
for(my $i = 6;$i<@fields;$i+=2)
{
	$sample_names[$i] = $fields[$i];
}
my $line_count = 0;
my $TYPE_SNP = 0;
my $TYPE_DEL = 1;
my $TYPE_INS = 2;
while(<FILE>)
{
	chomp;
	$lines[$line_count] = $_;
	@fields = split('\t');
	$chr[$line_count] = $fields[0];
	$pos[$line_count] = $fields[1];
	$type[$line_count] = $TYPE_SNP;
	my $name = "$fields[0]\_$fields[1]";
	if( ($fields[5] eq "INS") || ($fields[5] eq "DENOVO_INS") )
	{
		# print "\n Setting $name as needing indel \n";
		$type[$line_count] = $TYPE_INS;
		$ins_needed{$name} = 1;
	}
	elsif( ($fields[5] eq "DEL") || ($fields[5] eq "DENOVO_DEL") )
	{
		$type[$line_count] = $TYPE_DEL;
	}
	elsif( ($fields[5] eq "MULTIALLELIC") || ($fields[5] eq "DENOVO_MULTIALLELIC") )
	{
		my @sfields = split('\,',$fields[3]);
		foreach $i (@sfields)
		{
			if($i eq "I")
			{
				if($type[$line_count] != $TYPE_DEL)
				{
					$type[$line_count] = $TYPE_INS;
				}
				$ins_needed{$name} = 1;
			}
			elsif($i eq "D") 
			{
				$type[$line_count] = $TYPE_DEL;
			}
		}
	}
	$line_count++;
	if($line_count % 100000 == 0)
	{
		print "\n Read $line_count lines of the SNP file \n";
		#close(FILE);
	}
}
close(FILE);
for(my $i = 6;$i<@sample_names;$i+=2)
{
	my $all_lines;
	my $fh = new IO::Zlib;
    	$fh->open("$ARGV[2]/$sample_names[$i].indel.txt.gz", "rb")  || die "\n Can not open $sample_names[$i].indel.txt.gz \n";
	print "\n Working on file $sample_names[$i].indel.txt.gz \n";
	my @lines = $fh->getlines;
	chomp @lines;
	my $jj = @lines + 0;
	print "\n Found a total of $jj lines \n";
	# while(($_ = <gzFILE>) && $line_count < 20000)
	for(my $jj = 1;$jj<@lines;$jj++)
	{
		@fields = split('\t',$lines[$jj]);
		my $j = $fields[1] - $chr_offset{"$fields[0]"};
		my $name = "$fields[0]\_$j";
		#print "\n Testing $name to see if we need to save the insertion";
		if(exists($ins_needed{$name}))
		{
			#print "\n Saving data for $name \n";
			for($j = 7;$j<@fields;$j++)
			{
				my $jj = $fields[$j];
				$ins_a_count{$name}{"$jj"}++;
			}			
		}

	}
}
print "\n Making Consensus insertions \n";
foreach my $i (keys %ins_needed)
{
	if(!exists($ins_a_count{$i}))
	{
		print "\n This is impossible.  We fail to find insertion $i \n";
	}
	my $top = -50;
	foreach my $j (keys %{$ins_a_count{$i}})
	{
		if($ins_a_count{$i}{$j} > $top)
		{
			$ins_consensus{$i} = $j;
			$top = $ins_a_count{$i}{$j};
		}
	}
}
for(my $i = 0;$i<@pos;$i++)
{
	$sort_snp[$i] = $i;
}
print "\n About to Sort SNPs \n";

@ssort_snp = sort {compare_snp()} (@sort_snp);
open(FILE,">$ARGV[3]") || die "\nCan not open $ARGV[3] for writing \n"; 
print "\n Writing output \n";
print FILE "$header\n";
for($i = 0; $i < @ssort_snp;$i++)
{
	$j = @ssort_snp[$i];
	if($type[$j] == $TYPE_DEL)
	{
		@fields = split("\t",$lines[$j]);
		my $allele = 1;
		my $k = $i+1;
		my $name = "$fields[0]\_$fields[1]";
		while( ($pos[$ssort_snp[$k]] - $pos[$ssort_snp[$k-1]] == 1) && ($type[$ssort_snp[$k]] == $TYPE_DEL))
		{
			$allele++;
			$k++;
		}
		my $old_3 = $fields[3];
		$allele = "-$allele";
		# print "\n\tFor $name we have  deletion and I'm going to substitute in $allele into $fields[3]";
		$fields[3] =~ s/D/$allele/;
		if(exists($ins_consensus{$name}))
		{
			$allele = "+$ins_consensus{$name}";
			# print "\n\tFor $name we have insertion and I'm going to substitute in $allele into $fields[3]";
			$fields[3] =~ s/I/$allele/;
		}
		# print "\nFinally About to substitute $fields[3] for $old_3 \n";
		$lines[$j] =~ s/$old_3/$fields[3]/;
		$i = $k-1;
	}
	elsif($type[$j] == $TYPE_INS)
	{
		@fields = split("\t",$lines[$j]);
		my $old_3 = $fields[3];
		my $name = "$fields[0]\_$fields[1]";
		my $allele = "";
		if(exists($ins_consensus{$name}))
		{
			$allele = "+$ins_consensus{$name}";
			$fields[3] =~ s/I/$allele/;
			# print "\n\tFor $name we have insertion and I'm going to substitute in $allele into $fields[3]\n";
		}
		# print "\nFinally About to substitute $fields[3] for $old_3 \n";
		$lines[$j] =~ s/$old_3/$fields[3]/;
	}
	print FILE "$lines[$j]\n";
	
}

sub compare_snp ($a $b)
{
    if($chr_num{$chr[$a]} < $chr_num{$chr[$b]})
    {
	return -1;
    }
    if($chr_num{$chr[$b]} < $chr_num{$chr[$a]})
    {
	return 1;
    }
    if($pos[$a] < $pos[$b])
    {
	return -1;
    }
    if($pos[$b] < $pos[$a])
    {
	return 1;
    }
    return 0;
}
