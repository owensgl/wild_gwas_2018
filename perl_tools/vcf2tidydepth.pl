#!/bin/perl
use warnings;
use strict;
#VCF to depth in tidy format

my %sample;
print "chr\tpos\tsample\tdepth";
while(<STDIN>){
  chomp;
  if ($_ =~ m/##/){next;}
  if ($_ =~ m/#/){
    my @a = split(/\t/,$_);
    foreach my $i (9..$#a){
      $sample{$i} =$a[$i];
    }
  }else{
    my @a = split(/\t/,$_);
    my $chr = $a[0];
    my $pos = $a[1];
    foreach my $i (9..$#a){
      my @fields = split(/:/,$a[$i]);
      my $depth = $fields[2]; 
      print "\n$chr\t$pos\t$sample{$i}\t$depth";
    }
  }
}
