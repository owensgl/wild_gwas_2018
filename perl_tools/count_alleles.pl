#!/bin/perl
use strict;
use warnings;
my $indel_count = 0;
my $non_indel_count = 0;
my %allele_count;
while(<STDIN>){
  chomp;
  if ($_ =~ m/^#/){
    next;
  }
  my @a = split(/\t/,$_);
  my $ref = $a[3];
  my @alt = split(/,/,$a[4]);
  #check if indel
  my $indel;
  if (length($ref) > 1){
    $indel++;
  }
  foreach my $i (0..$#alt){
    if (length($alt[$i]) > 1){
      $indel++;
    }
    if ($alt[$i] eq '*'){
      $indel++;
    }
  }
  if ($indel){
    $indel_count++;
  }else{
    $non_indel_count++;
  }
  $allele_count{($#alt+1)}++;
  
}
print "type\tcount";
print "\nindels\t$indel_count";
print "\nSNPs\t$non_indel_count";
foreach my $i (1..1000){
  if ($allele_count{$i}){
    print "\nalelles.$i\t$allele_count{$i}";
  }
}
