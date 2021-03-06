#!/usr/bin/perl

#################### Usage ###########################
###### perl PTM_summarizer.pl ref.fasta PSM.txt ######
######################################################

use strict;
use warnings;

my $ref_fasta=$ARGV[0];
my $infile=$ARGV[1];

my %hash=();
my %hashCount=();
my %hashLine=();
my %hash4upset=();
my %hashUniqPro=();
my @allM="";
#my @allM_Pro="";

my $ref_id="";
my %ref_Seq=();
open(REF,"$ref_fasta") or die "Could not open the file:$!\n";
while(<REF>)
{
	chomp;
	$_=~s/\r//g;
	if(/^>/)
	{
		$ref_id = (split /\|/)[1]; ## Header Parse for Master protein accession
		$ref_id =~s/\>//;
	}
	else
	{
		$ref_Seq{$ref_id}.="$_";
	}
}
close REF;

open(IN,"$infile") or die "Could not open the file:$!\n";
while(<IN>)
{
	chomp;
	$_=~s/\r//g;
	$_=~s/\"//g;
	unless(/^Annotated/)
	{
		my $line=$_;
		my ($Seq,$mas,$status,$fs,$ptmRS)=(split /\t/)[3,7,15,26,-1];#Sequence, Master Protein, mz, First scan, ptmRS
		#if(10==10)#if($status eq "Selected")
		{
			$mas=~s/\;\s/\;/g;
			my $annSeq=(split /\./,$Seq)[1];
			$ptmRS=~s/\s//g;
			my @allMod=(split /\;/,$ptmRS);
			foreach my $i (@allMod)
			{
				if($i=~/([A-Z])([0-9]*?)\((.*?)\)\:([0-9]*)/)
				{
					my ($aa,$pos,$ty,$rs)=($1,$2,$3,$4);
					if($rs >= 75)
					{
						my $k=$annSeq."_".$mas."_".$aa.$pos;
						unless(($ty eq "Oxidation") || (($ty eq "Deamidated") && (($aa eq "N")||($aa eq "Q"))))
						{
							foreach my $m (split /\;/,$mas)
							{
								#print "$m\n";	
								my $Proseq = $ref_Seq{$m};
								my $ann = $annSeq; $ann=~tr/[a-z]/[A-Z]/;
								my $index = (index ($Proseq, $ann) + $pos);
								#print OUT "$annSeq\t$mas\t$m\t$ty\t$aa\t$index\n";
								my $uniqPro=$m."_".$ty."_".$aa."_".$index;
								my $ty_aa = $ty."_".$aa;
								$hashUniqPro{$uniqPro}="$annSeq\t$m\t$ty\t$aa\t$ty_aa\t$index";
								#push(@allM_Pro,$ty);
							}
							unless(exists $hash{$k}{$ty})
							{						
								$hash{$k}{$ty}="";
								$hashCount{$ty}{$aa}++;
								$hashLine{$ty}{$aa}.="$_\n";
								my $mod=$ty."_".$aa;
								my $ty2=$ty;$ty2=~s/\+HeavyK//g;
								$hash4upset{$annSeq}{$ty2}.="";
								push(@allM,$ty2);
							}							
						}							
					}	
				}
			}
		}	
	}
}	
close IN;

my $outun = $infile; $outun =~s/\.txt$/\_UniqueProteinSite.txt/;
open(OUTUNIQPRO,">$outun") or die "Could not create the file:$!\n";
foreach my $k1 (keys %hashUniqPro)
{
	print OUTUNIQPRO "$hashUniqPro{$k1}\n";
}
close OUTUNIQPRO;

my $out1 = $infile; $out1 =~s/\.txt$/\_summary.txt/;
open(OUT1,">$out1") or die "Could not create the file:$!\n";
foreach my $k1 (keys %hashCount)
{
	print OUT1 "$k1\n";
	foreach my $k2 (keys %{$hashCount{$k1}})
	{
		print OUT1 "$k1\t$k2\t$hashCount{$k1}{$k2}\n";
		my $file=$k1."_".$k2.".txt";
		open(OUT,">$file") or die "Could not create the file:$!\n";
		print OUT "$hashLine{$k1}{$k2}\n";
		close OUT;
	}
	#print OUT1 "\n";
}
close OUT1;

my @unique = do { my %seen; grep { !$seen{$_}++ } @allM };
open(OUT2,">ForUpset.txt") or die "Could not create the file:$!\n";
print OUT2 "Peptide";
foreach (@unique)
{
	print OUT2 "\t$_";
}
print OUT2 "\n";
foreach my $k1 (keys %hash4upset)
{
	print OUT2 "$k1";
	foreach my $k2 (@unique)
	{
		if((exists $hash4upset{$k1}) && (exists $hash4upset{$k1}{$k2}))
		{
			print OUT2 "\t1";
		}
		else
		{
			print OUT2 "\t0";
		}
	}
	print OUT2 "\n";
}
close OUT2;
exit;
