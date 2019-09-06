#!/bin/perl
use strict;
use warnings;

#This script takes a VCF, a list of outgroup samples, and outputs only sites where the outgroups are fixed for a single allele. It then polarizes all SNPs as A (ancestral) or B (derived). Requires biallelic sites only.
#In this version the outgroup is in the vcf itself
my $min_out_sampled = 0.45;

my $popfile = $ARGV[0]; #List of outgroup samples
my %outgroups;
open POP, $popfile;
while(<POP>){
  chomp;
  $outgroups{$_}++;
}
close POP;

my %sample;
my @outgroup_columns;
while(<STDIN>){
  chomp;
  my @a = split(/\t/,$_);
  if ($_ =~ m/^#CHR/){
    print "chr\tpos";
    foreach my $i (9..$#a){
      $sample{$i} = $a[$i];
      if ($outgroups{$a[$i]}){
        push(@outgroup_columns,$i);
      }else{
        print "\t$a[$i]";
      }
    }
    next;
  }
  if ($_ =~ m/^#/){next;}
  my $alts = $a[4];
  my @alts = split(/,/,$alts);
  if ($alts[1]){next;} #Skips sites with more than two alleles
  my %outgroup_alleles;
  my $outgroup_present= 0;
  my $outgroup_absent = 0;
  foreach my $i (@outgroup_columns){
    if (($a[$i] eq '.') or ($a[$i] eq '.:0,0')){
      $outgroup_absent++;
    }else{
      my @infos = split(/:/,$a[$i]);
      if ($infos[0] eq './.'){
        $outgroup_absent++;
        next;
      }
      my @genos = split(/\//,$infos[0]);
      foreach my $a (0..1){
        $outgroup_alleles{$genos[$a]}++;
      }
      $outgroup_present++;
    }
  }
  my $outgroup_called_percent = $outgroup_present / ($outgroup_absent + $outgroup_present);
  if ($outgroup_called_percent < $min_out_sampled){next;} #Skip rows with too many missing outgroups
  my $n_outgroup_alleles = scalar keys %outgroup_alleles;
  if ($n_outgroup_alleles != 1){next;} #Skip rows with polymorphic outgroups
  my @tmp = sort keys %outgroup_alleles;
  my $outgroup_allele = $tmp[0];

  print "\n$a[0]\t$a[1]";
  #Now to figure out if the samples have derived or ancestral
  foreach my $i (9..$#a){
    if ($outgroups{$sample{$i}}){next;} #Don't analyze outgroup samples.
    if (($a[$i] eq '.') or ($a[$i] eq '.:0,0')){
      print "\tN";
    }else{
      my @infos = split(/:/,$a[$i]);
      if ($infos[0] eq './.'){
        print "\tN";
        next;
      }
      my @genos = split(/\//,$infos[0]);
      if (($genos[0] == $outgroup_allele ) and ($genos[1] == $outgroup_allele )){
        print "\tA";
      }elsif (($genos[0] != $outgroup_allele ) and ($genos[1] != $outgroup_allele)){
        print "\tB";
      }else{
        print "\tH";
      }
    }
  }
}
