#!/bin/perl
use strict;
use warnings;
use List::Util 'shuffle';
#This script takes the population info file and inversion genotypes.
#It pulls out populations that have two samples each homozygote, and selects a random set of them.

my $info_file = "/scratch/gowens/wild_gwas/wild_gwas_2018/sample_info_apr_2018.tsv";
my $inv_file = $ARGV[0];
my $samples_in_vcf = $ARGV[1];
my $chosen_species = $ARGV[2];
my $chosen_type = "wild";
open INFO, $info_file;

my %pop;
my %pop_list;
while(<INFO>){
  chomp;
  if ($. == 1){next;}
  my @a = split(/\t/,$_);
  my $sample = $a[2];
  my $species = $a[3];
  my $pop = $a[4];
  my $type = $a[5];
  unless ($species =~ m/$chosen_species/){next;}
  if ($type ne $chosen_type){next;}
  if ($pop eq "NA"){next;}
  $pop{$sample} = $pop;
  $pop_list{$pop}++;
}
close INFO;

my %genotype;
open INV, $inv_file;
while(<INV>){
  chomp;
  if ($. == 1){next;}
  my @a = split(/\t/,$_);
  my $sample = $a[0];
  my $perc_1 = $a[3];
  my $genotype = $a[5];
  if (($perc_1 < 0.9) and ($perc_1 > 0.1)){
    next;
  }
  $genotype{$sample} = $genotype;
}

my %sampled;
open SAMPLED, $samples_in_vcf;
while(<SAMPLED>){
  chomp;
  $sampled{$_}++;
}
close SAMPLED;	

my $zero_printed = 0;
my $two_printed = 0;
foreach my $sample ( shuffle keys %genotype ) {
  if ($sampled{$sample}){
    if (($zero_printed < 1) and ($genotype{$sample} eq "0") and ($pop{$sample})){
      print "$sample\n";
      $zero_printed++;
    }elsif (($two_printed < 1) and ($genotype{$sample} eq "2") and ($pop{$sample})){
      print "$sample\n";
      $two_printed++;
    }
  }
}

