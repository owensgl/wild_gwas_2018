#!/bin/perl
use warnings;
use strict;

#This script takes the list of mds outlier windows, and picks windows in there by merging neighbouring windows.

my $mds_file = $ARGV[0]; #e.g. Anomalus.tranche90.snp.gwas.90.bi.500.mds_cluster_windows.txt

open MDS, $mds_file;
my $current_mds = "NA";
my $previous_n;
my $previous_end;
while(<MDS>){
  chomp;
  if ($. == 1){next;}
  my @a = split(/\t/,$_);
  my $chr = $a[0];
  my $start = $a[1];
  my $end = $a[2];
  my $n = $a[6];
  my $mds_coord = $a[7];
  if ($current_mds eq "NA"){
    print "$chr:$start";
    $current_mds = $mds_coord;
    $previous_n = $n;
    $previous_end = $end;
  }elsif ($current_mds eq $mds_coord){
    if (($n - 1) == $previous_n){
      $previous_n = $n;
      $previous_end = $end;
    }else{
      print "-$previous_end\t$current_mds\n";
      print "$chr:$start";
      $previous_n = $n;
      $previous_end = $end;
    }
  }else{
    print "-$previous_end\t$current_mds\n";
    print "$chr:$start";
    $current_mds = $mds_coord;
    $previous_n = $n;
    $previous_end = $end;
  }
}
print "-$previous_end\t$current_mds\n";
