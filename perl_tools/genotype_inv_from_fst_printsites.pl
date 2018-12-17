#!/bin/perl
use strict;
use warnings;

#This script takes an fst file, it picks out the sites that are above an Fst threshold and records what alleles are for each group. It then uses a vcf file to record how many of each type of allele each sample has. 
my $mds_fst_file = $ARGV[0]; #List of fst outliers
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
  my $mds = $a[10];
  my $inv_0_allele = $a[8];
  my $inv_2_allele = $a[9];
  if ($fst < $min_fst){next;}
  $allele_1{$chr}{$pos}{$mds} = $inv_0_allele;
  $allele_2{$chr}{$pos}{$mds} = $inv_2_allele;
}
close MDS;
open(IN, "gunzip -c $gzvcf |");

my %sample;
my %counts;
my $site_count = 0;
print "chr\tpos\tmds_coord\tsample\tgenotype\tdepth";
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
      my $dp = $fields[2];
      if (($genotype eq '.') or ($genotype eq './.')){next;}
      if (($genotype eq "$allele_1\/$allele_2") or  ($genotype eq "$allele_2\/$allele_1")){
	$counts{$sample{$i}}{$mds}{1}++;
        print "\n$chr\t$pos\t$mds\t$sample{$i}\t01\t$dp";
      }elsif ($genotype eq "$allele_1\/$allele_1"){
	print "\n$chr\t$pos\t$mds\t$sample{$i}\t00\t$dp";      
      }elsif ($genotype eq "$allele_2\/$allele_2"){
	print "\n$chr\t$pos\t$mds\t$sample{$i}\t11\t$dp";
      }
    }
  }
}  
