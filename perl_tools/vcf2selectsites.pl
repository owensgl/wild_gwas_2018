#!/bin/perl
use strict;
use warnings;
my $site_file = $ARGV[0];

my %chosen_sites;
open FILE, $site_file;
while(<FILE>){
  chomp;
  my @a = split(/\t/,$_);
  $chosen_sites{"$a[0].$a[1]"}++;
}

while(<STDIN>){
  chomp;
  if ($. == 1){
    print "$_";
    next;
  }
  if ($_ =~ /^#/){
    print "\n$_";
  }else{
    my @a = split(/\t/,$_);
    if ($chosen_sites{"$a[0].$a[1]"}){
      print "\n$_";
    }
  }
}
