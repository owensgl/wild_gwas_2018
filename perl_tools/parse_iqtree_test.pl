#!/usr/bin/perl
#Takes the iqtree output and parses out tree support values

use strict;
use warnings;

my @fields;
my $counter;
my %results;
while(<STDIN>){
  chomp;
  if ($_ =~ /^deltaL/){goto PRINT;}
  if (($_ =~ m/Tree/) and ($_ =~ m/logL/)){
    @fields = split(' ',$_);
    $counter++;
    next;
  }
  if ($counter){
    if ($_ =~ /---/){next;}
    my @a = split(' ',$_);
    unless($a[1]){next;}
    foreach my $i (1..$#a){
      $results{$a[0]}{$i} = $a[$i];
    }
  }
}
PRINT:
print "$fields[0]\t$fields[1]\t$fields[2]\t$fields[3]\t$fields[3].conf\t";
print "$fields[4]\t$fields[4].conf\t$fields[5]\t$fields[5].conf\t";
print "$fields[6]\t$fields[6].conf\t$fields[7]\t$fields[7].conf";

foreach my $tree (sort keys %results){
  print "\n$tree";
  foreach my $i (1..12){
    print "\t$results{$tree}{$i}";
  }
}
