#!/usr/bin/perl
use strict;
use warnings;

# Read in sam files
open(my $fi,"<","data.Sc_map1_unq.sam") or die "Could not open input";

# Specify output files
open(my $fo,">","output.fq") or die "Could not open output";

# Loop through mapped reads
while(my $line=<$fi>)
{	chomp($line);
	my @a=split(/\t/,$line);
	$a[5]=~s/M//;
	my $read_length=$a[5];
	my $num_mismatch=chop $a[$#a];
	my $read_name=$a[0];
	my $strand=$a[1];
	my $read=$a[9];
	my $score=$a[10];
	my $mismatch=$a[12];

	if($num_mismatch>0)
	{	my @mm_split=split(/[A-Z]/,$mismatch);
		if($strand==0 && substr($mismatch,5,1)==0)	# Check if the first character on the positive strand is a mismatch 
		{	$read_length--;
			$read=substr($read,1,$read_length);		# Delete the first nucleotide of the read
			$score=substr($score,1,$read_length);
		}
		elsif($strand==16 && $mm_split[$#mm_split]==0)	# Check if the last character on the negative strand read is a mismatch 
		{	$read_length--;
			$read=substr($read,0,$read_length);
			$score=substr($score,0,$read_length);
		}
	}		

	if($strand==16) # Check if the read maps to the negative strand
	{	$read=reverse($read);	# Reverse the negative strand
		$read=~tr/ATGC/TACG/;	# Complement of the negative strand
		$score=reverse($score);
	} 
	
	print $fo "@","$read_name\n","$read\n+\n$score\n";
}
close($fi);
close($fo);
