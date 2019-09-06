#!/bin/perl
use strict;
use warnings;
use POSIX;
#use Statistics::Basic qw(:all);
#This is for calculating Dstat in inversion windows
#This takes a DA file (Derived-Ancestral), the names of samples 1 2 and 3, and then the genomic location of the inversions

my $sample1 = $ARGV[0];
my $sample2 = $ARGV[1];
my $sample3 = $ARGV[2];
my $species = $ARGV[3]; #annuus, argophyllus, petiolaris
my $chr_chosen = $ARGV[4]; #chr
my $mds_chosen = $ARGV[5]; #mds

my $region_file = "/scratch/gowens/wild_gwas/wild_gwas_2018/MDS_outliers/Ha412HO/all/Ha412HO_inv.v3.inversions.regions.v1.txt";
my %inv_regions;
open REGION, $region_file;
while(<REGION>){
	chomp;
	if ($. == 1){next;}
	my @a = split(/\t/,$_);
	if (($a[3] eq $chr_chosen) and ($species eq $a[4]) and ($a[7] eq $mds_chosen)){
		$inv_regions{$a[3]}{$a[0]}{'start'} = $a[1];
		$inv_regions{$a[3]}{$a[0]}{'end'} = $a[2];
	}
}
close REGION;
my $block_size = 5000000; #Block jackknife size
my @used_cols;
my %samples;
my %hash;
my %inv_hash;
while(<STDIN>){
  chomp;
  my @a = split(/\t/,$_);
  if ($. == 1){
    foreach my $i (2..$#a){
      if ($a[$i] eq $sample1){
	push(@used_cols, $i);
	$samples{1} = $i;
      }elsif ($a[$i] eq $sample2){
        push(@used_cols, $i);
	$samples{2} = $i;
      }elsif ($a[$i] eq $sample3){
        push(@used_cols, $i);
        $samples{3} = $i;
      }
    }
  }else{
    if ($a[$samples{3}] ne "B"){next;}
    my $chr = $a[0];
    my $pos = $a[1];
    my $block = floor($pos/$block_size);
    if ($inv_regions{$chr}){
      foreach my $i (1..20){
	if ($inv_regions{$chr}{$i}{'start'}){
	  if (($pos >= $inv_regions{$chr}{$i}{'start'}) & ($pos <=$inv_regions{$chr}{$i}{'end'})){
	     $inv_hash{"$a[$samples{1}]$a[$samples{2}]$a[$samples{3}]"}++;
	     next;
	  }
	}
      }
    }
    $hash{$chr.$block}{"$a[$samples{1}]$a[$samples{2}]$a[$samples{3}]"}++;
  }
}
#Run Jackknife
my @jackknife_values;
foreach my $skipped_block (sort keys %hash){
  my $abba = 0;
  my $baba = 0;
  foreach my $block (sort keys %hash){
    if ($block eq $skipped_block){next;}
    if ($hash{$block}{'ABB'}){
      $abba+=$hash{$block}{'ABB'};
    }
    if ($hash{$block}{'BAB'}){
      $baba+=$hash{$block}{'BAB'};
    }
  }
  my $dstat = ($abba - $baba) / ($abba + $baba);
  push(@jackknife_values, $dstat);
}
my $mean_jackknife = Mean(\@jackknife_values);
my $var_jackknife = variance(\@jackknife_values);

unless($inv_hash{'ABB'}){$inv_hash{'ABB'} = 0};
unless($inv_hash{'BAB'}){$inv_hash{'BAB'} = 0};
my $inv_Dstat = "NA";
if ( ($inv_hash{'ABB'} + $inv_hash{'BAB'}) > 0){
  $inv_Dstat = ($inv_hash{'ABB'}  - $inv_hash{'BAB'}) / ($inv_hash{'ABB'} + $inv_hash{'BAB'});
}
print "$chr_chosen\t$mds_chosen\t$species\t$sample1\t$sample2\t$sample3\t$mean_jackknife\t$var_jackknife\t$inv_Dstat\t$inv_hash{'ABB'}\t$inv_hash{'BAB'}\n";


#Variance and mean subroutines

sub sum
{
        my ($arrayref) = @_;
        my $result;
        foreach(@$arrayref) { $result+= $_; }
        return $result;
}

sub Mean {
    my ($arrayref) = @_;
    my $result;
    foreach (@$arrayref) { $result += $_ }
    return $result / @$arrayref;
}


sub variance
{
	my $mean = Mean(@_);
	return (sum [ map { ($_ - $mean)**2 } @{$_[0]}  ] ) / $#{$_[0]};
}



