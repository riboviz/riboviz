#!/usr/bin/perl
use strict;
use warnings;

# Read in sam files
open(my $fi,"<",$ARGV[0]) or die "Could not open input";

# Specify output files
open(my $fo,">",$ARGV[1]) or die "Could not open output";

my $trim_pos=0;
my $trim_neg=0;
# Loop through mapped reads
while(my $line=<$fi>)
{	chomp($line);
	my @a=split(/\t/,$line);
	if(substr($line,0,1) eq "@")
	{	print $fo "$line\n";
		next;
	}elsif($a[1]==4)
	{	next;
	}

	$a[5]=~s/M//;
	my $read_length=$a[5];
	my $num_mismatch=chop $a[$#a];
	my $read_name=$a[0];
	my $strand=$a[1];
	my $read=$a[9];
	my $score=$a[10];
	my $mismatch=$a[12];
	my $xa_tag=$a[11];

	if($num_mismatch>0)
	{	my @mm_split=split(/[A-Z]/,$mismatch);
		if($strand==0 && substr($mismatch,5,1)==0)	# Check if the first character on the positive strand is a mismatch 
		{	$read_length--;
			$read=substr($read,1,$read_length);		# Delete the first nucleotide of the read
			$score=substr($score,1,$read_length);
			$mismatch=substr($mismatch,0,5).substr($mismatch,7,(length($mismatch)-7));
			$num_mismatch--;
			$xa_tag=substr($xa_tag,0,(length($xa_tag)-1)).$num_mismatch;
			$trim_pos++;
		}
		elsif($strand==16 && $mm_split[$#mm_split]==0)	# Check if the last character on the negative strand read is a mismatch 
		{	$read_length--;
			$read=substr($read,0,$read_length);
			$score=substr($score,0,$read_length);
			$mismatch=substr($mismatch,0,(length($mismatch)-2));
			$num_mismatch--;
			$xa_tag=substr($xa_tag,0,(length($xa_tag)-1)).$num_mismatch;
			$trim_neg++;
		}
	}		

	print $fo "$read_name\t$strand\t$a[2]\t$a[3]\t$a[4]\t",$read_length."M","\t$a[6]\t$a[7]\t$a[8]\t$read\t$score\t$xa_tag\t$mismatch\t",$a[13].$num_mismatch,"\n";
}
print "+ strand trimmed = $trim_pos\n";
print "- strand trimmed = $trim_neg\n";
close($fi);
close($fo);
