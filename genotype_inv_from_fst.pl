#!/bin/perl
use strict;
use warnings;

#This script takes a piped in fst file, it picks out the sites that are above an Fst threshold and records what alleles are for each group. It then uses a vcf file to record how many of each type of allele each sample has. 
#Assumes that the vcf is biallelic. 
my $mds_fst_file = $ARGV[0];
my $gzvcf = $ARGV[1];
my $min_fst = 0.95;
open MDS, $mds_fst_file;

my %allele_1;
my %allele_2;
while(<MDS>){
  chomp;
  if ($. == 1){next;}
  my @a = split(/\t/,$_);
  my $chr = $a[0];
  my $pos = $a[1];
  my $fst = $a[6];
  my $mds = $a[9];
  my $inv_0_allele = $a[7];
  my $inv_2_allele = $a[8];
  if ($fst < $min_fst){next;}
  $allele_1{$chr}{$pos}{$mds} = $inv_0_allele;
  $allele_2{$chr}{$pos}{$mds} = $inv_2_allele;
}
close MDS;
open(IN, "gunzip -c $gzvcf |");

my %sample;
my %counts;
my $site_count = 0;
while(<IN>){
  chomp;
  if ($_ =~ m/^##/){
    next;
  }
  my @a = split(/\t/,$_);
  if ($_ =~ m/^#/){
    print STDERR "Loading vcf file now...\n";
    foreach my $i (9..$#a){
      $sample{$i} = $a[$i];
    }
    next;
  }
  my $chr = $a[0];
  my $pos = $a[1];
  $site_count++;
  if ($site_count % 100000 == 0){print STDERR "$chr $pos processed...\n";}
  unless ($allele_1{$chr}{$pos}){next;}
  my $ref = $a[3];
  my $alt = $a[4];
  my @alts = split(/,/,$alt);
#  print "Found outlier at $chr $pos\n";
  foreach my $mds (sort keys $allele_1{$chr}{$pos}){
    #Check to make sure that the alleles match;
    #Find out which allele in this site matches with the inversion alleles
    my $allele_1;
    my $allele_2;
    if ($ref eq $allele_1{$chr}{$pos}{$mds}){
      $allele_1 = 0;
    }elsif ($ref eq $allele_2{$chr}{$pos}{$mds}){
      $allele_2 = 0;
    }
    foreach my $x (0..$#alts){
       if ($alts[$x] eq $allele_1{$chr}{$pos}{$mds}){
         $allele_1 = ($x+1);
       }elsif ($alts[$x] eq $allele_2{$chr}{$pos}{$mds}){
         $allele_2 = ($x+1);
       }
     }
    foreach my $i (9..$#a){
      my @fields = split(/:/,$a[$i]);
      my $genotype = $fields[0];
      if (($genotype eq '.') or ($genotype eq './.')){next;}
      if (($genotype eq "$allele_1\/$allele_2") or  ($genotype eq "$allele_2\/$allele_1")){
	$counts{$sample{$i}}{$mds}{1}++;
      }elsif ($genotype eq "$allele_1\/$allele_1"){
	$counts{$sample{$i}}{$mds}{0}++;        
      }elsif ($genotype eq "$allele_2\/$allele_2"){
	$counts{$sample{$i}}{$mds}{2}++;
      }
    }
  }  
}
print "name\tmds_coord\t00\t01\t11";
foreach my $sample (sort keys %counts){
  foreach my $mds (sort keys %{$counts{$sample}}){
    foreach my $i (0..2){
      unless($counts{$sample}{$mds}{$i}){$counts{$sample}{$mds}{$i} = 0;}
    }
    print "\n$sample\t$mds\t$counts{$sample}{$mds}{0}\t$counts{$sample}{$mds}{1}\t$counts{$sample}{$mds}{2}";
  }
}

