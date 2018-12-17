#!/bin/perl
use warnings;
use strict;

#This takes the fst.txt.gz (including extra two columns for XRQ, and outputs the sites above a threshold for Fst and 0/1 heterozygosity)

my $min_het = 0.5;
my $min_fst = 0.8;

while(<STDIN>){
  chomp;
  if ($. == 1){
    print "$_";
  }else{
    my @a = split(/\t/,$_);
    my $fst = $a[8];
    my $het = $a[12];
    if ($fst >= $min_fst){
      if ($min_het >= $min_het){
        print "\n$_";
      }
    }
  }
}
