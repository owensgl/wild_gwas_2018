#!/bin/perl 
use warnings;
use strict;
use POSIX;

#This script takes a fst.txt script, and calculates windows  with high percentages of highly differentiated snps, and then extends around them to create large regions.
#It uses the frequency difference instead of fst.

my $mds_coord = $ARGV[0];
my $min_freq_dif = 0.75;
my $min_outlier_perc = 0.1;
my $spanning_range = 2000000;
my $window_size = 10000;

my %outliers;
my %total_sites;
my %outlier_chrs; #Count outliers per chromosome to pick main target.
while(<STDIN>){
  chomp;
  if ($. == 1){next;}
  my @a = split(/\t/,$_);
  my $chr = $a[0];
  my $pos = $a[1];
  my $window = floor($pos/$window_size)*$window_size;
  my $freq_dif = $a[7];
  my $status = 0;
  if ($freq_dif >= $min_freq_dif){
    $status = 1;
    $outlier_chrs{$chr}++;
  }
  $outliers{$chr}{$window}+= $status;
  $total_sites{$chr}{$window}++;
}
print STDERR "Loaded fst file...\n";
#Main_chr is the one with the most freq_dif outliers;
my @outlier_chrs = sort { $outlier_chrs{$b} <=> $outlier_chrs{$a} } keys %outlier_chrs;
my $main_chr = $outlier_chrs[0];



#Find all windows with more high freq_dif sites than threshold
my %outlier_windows;
foreach my $chr (sort keys %outliers){
  foreach my $window (sort {$a <=> $b} keys %{$outliers{$chr}}){
    my $perc_outlier = $outliers{$chr}{$window} / $total_sites{$chr}{$window};
    if ($perc_outlier >= $min_outlier_perc){
      $outlier_windows{$chr}{$window}++;
    }
  }
}

#Expand those regions by merging near windows into larger regions
my %outlier_regions; #Key = start, value = end;

foreach my $chr (sort keys %outlier_windows){
  foreach my $window (sort {$a <=> $b} keys %{$outlier_windows{$chr}}){
    my $start = $window;
    my $end = $window + $window_size;
    my $found_neighbour;
    foreach my $other_window (sort {$a<=>$b} keys %{$outlier_regions{$chr}}){
      if (($outlier_regions{$chr}{$other_window} + $spanning_range) >= $start){
        $outlier_regions{$chr}{$other_window} = $end;
        $found_neighbour++;
      }
    }
    unless($found_neighbour){
      $outlier_regions{$chr}{$start} = $end;
    }
  }
}

print "chr\tstart\tend\tclass\tmds-coord";
foreach my $chr (sort keys %outlier_regions){
  foreach my $window (sort {$a <=> $b} keys %{$outlier_regions{$chr}}){
    my $class = "secondary";
    if ($chr eq $main_chr){
      $class = "primary";
    }
    print "\n$chr\t$window\t$outlier_regions{$chr}{$window}\t$class\t$mds_coord";
  }
}





