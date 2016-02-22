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
use vars qw(@fields %tran_map @filehandles %ts %tv %transition %transver @files $i $j $last @samples);

if(@ARGV != 1)
{
	die "\n Usage: ${0} SNPFILE\n";
}

open(FILE,"$ARGV[0]") || die "\n Could not open $ARGV[0] which should be the snp file \n";
$_ = <FILE>;
chomp;
@fields = split('\t');
my $i, $j;
$i = 0;
for($j=6;$j<@fields;$j+=2,$i++)
{
	$samples[$i] = $fields[$j];
}
$tran_map{"AG"} = 1;
$tran_map{"GA"} = 1;
$tran_map{"CT"} = 1;
$tran_map{"TC"} = 1;

while(<FILE>)
{
	chomp;
	@fields = split('\t');
		my $is_trans = 0;
		if( ($fields[3] eq "A,G") || ($fields[3] eq "C,T") || ($fields[3] eq "G,A") || ($fields[3] eq "T,C"))
		{
			$is_trans = 1;
		}
		elsif(exists($tran_map{"$fields[2]$fields[3]"}))
		{
			$is_trans = 1;
		}
		if($is_trans)
		{
			$transition{$fields[5]}++;
		}
		else
		{
			$transver{$fields[5]}++;
		}
		#print "$fields[3] $is_trans $transition $transver\n";
		$i=0;
		for($j=6;$j<@fields;$j+=2,$i++)
		{
			if($fields[$j] ne $fields[2] && $fields[$j] ne "N")
			{
				if($is_trans)
				{
					$ts{$fields[5]}[$i]++;
				}	
				else
				{
					$tv{$fields[5]}[$i]++;
				}
			}
		}
}
print "Category";
my @types = sort (keys %transver);
foreach my $i (@types)
{
	print "\t$i\_Transitions\t$i\_Transversion\t$i\_ratio";
}
print "\nALL";
foreach my $i (@types)
{
	$j = 1;
	$transition{$i} += 0;
	$transver{$i} += 0;
	if($transver{$i} > 0)
	{
		$j = $transition{$i} / $transver{$i};
	}
	print "\t$transition{$i}\t$transver{$i}\t$j";
}
for($i=0;$i<@samples;$i++)
{
	print "\n$samples[$i]";
	foreach my $k (@types)
	{
		$j = 1;
		$tv{$k}[$i] += 0;
		$ts{$k}[$i] += 0;
		if($tv{$k}[$i] > 0)
		{
			$j = $ts{$k}[$i] / $tv{$k}[$i];
		}
		print "\t$ts{$k}[$i]\t$tv{$k}[$i]\t$j";
	}
}
print "\n";
