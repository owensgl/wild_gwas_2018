#!/bin/perl
use warnings;
use strict;
#This takes the invOut and converts it to a format for R
my $chr = $ARGV[0];
my $inv_info;
print "inv_chr\tref1_start\tref1_end\tref2_start\tref2_end\tregion1_start\tregion1_end\tregion2_start\tregion2_end";
while(<STDIN>){
  chomp;
  if ($_ =~ m/^#/){
    my @a = split(' ',$_);
    if ($#a == 7){
      $inv_info = "$chr\t$a[2]\t$a[3]\t$a[6]\t$a[7]";
    }elsif ($#a == 5){
      $inv_info = "$chr\t$a[1]\t$a[2]\t$a[4]\t$a[5]";
    }
    next;
  }else{
    my @a = split(' ',$_);
    print "\n$inv_info\t$a[0]\t$a[1]\t$a[2]\t$a[3]";
  }
}
