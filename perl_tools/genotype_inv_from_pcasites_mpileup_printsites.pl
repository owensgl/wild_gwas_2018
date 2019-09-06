#!/bin/perl
use strict;
use warnings;

#This script takes an fst file, it picks out the sites that are above an Fst threshold and records what alleles are for each group. It then uses a vcf file to record how many of each type of allele each sample has. 
my $mds_fst_file = $ARGV[0]; #List of fst outliers
my $mpileup = $ARGV[1]; #samtools mpileup file
my $sample = $ARGV[2]; #Sample name
#my $min_freq = 0.95;
my $Ha412_remap = "TRUE"; #Set to TRUE to output the remapped positions

open (MDS, "gunzip -c $mds_fst_file |");

my %allele_1;
my %allele_2;
my %ha412_chr;
my %ha412_pos;
while(<MDS>){
  chomp;
  if ($. == 1){next;}
  my @a = split(/\t/,$_);
  my $chr = $a[2];
  my $pos = $a[3];
  my $ha412_chr = $a[0];
  my $ha412_pos = $a[1];
  $ha412_chr{$chr.$pos}= $ha412_chr;
  $ha412_pos{$chr.$pos}= $ha412_pos;
  my $inv_0_allele = $a[7];
  my $inv_2_allele = $a[8];
  $allele_1{$chr}{$pos} = $inv_0_allele;
  $allele_2{$chr}{$pos} = $inv_2_allele;
}
close MDS;
open(IN, "cat $mpileup |");

my %sample;
my %counts;
my $site_count = 0;
print "XRQchr\tXRQpos\tchr\tpos\tsample\tgenotype\tdepth";
while(<IN>){
  chomp;
  my @a = split(/\t/,$_);
  my $chr = $a[0];
  my $pos = $a[1];
  my $bases = uc $a[4];
  $site_count++;
  if ($site_count % 100000 == 0){print STDERR "$chr $pos processed...\n";}
  unless ($allele_1{$chr}{$pos}){next;}
  my $allele_1_count = 0;
  my $allele_2_count = 0;
  my @bases = split(//,$bases);
  foreach my $i (0..$#bases){
	  if ($bases[$i] eq $allele_1{$chr}{$pos}){
		  $allele_1_count++;
	}elsif ($bases[$i] eq $allele_2{$chr}{$pos}){
		$allele_2_count++;
	}
  }
  my $genotype;
  if (($allele_1_count >= 2) and ($allele_2_count >= 2)){
	  $genotype = "01";
  }elsif ($allele_1_count >= 2){
	  $genotype = "00";
  }elsif ($allele_2_count >= 2){
	  $genotype = "11";
  }else{
	  next;
  }
  my $dp = $allele_2_count + $allele_1_count;
  print "\n$chr\t$pos\t$ha412_chr{$chr.$pos}\t$ha412_pos{$chr.$pos}\t$sample\t$genotype\t$dp";
}
  
