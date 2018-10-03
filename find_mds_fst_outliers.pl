#!/bin/perl 
use warnings;
use strict;

#This script takes mds outlier file, and an Fst file, and pulls out high fst sites within those windows.

my $mds_file = $ARGV[0]; #MDS outlier file mds_cluster_windows.txt
my $fst_file_prefix = $ARGV[1]; #FST outlier prefix. e.g. Annuus.tranche90.snp.env.90.bi
my $min_fst = 0.90;

open MDS, $mds_file;

my %outlier_start;
my %outlier_end;
while(<MDS>){
  chomp;
  if ($. == 1 ){next;}
  my @a = split(/\t/,$_);
  my $chr = $a[0];
  my $start = $a[1];
  my $end = $a[2];
  my $n = $a[6];
  my $mds = $a[7];
  $outlier_start{$mds}{$chr}{$n} = $start;
  $outlier_end{$mds}{$chr}{$n} = $end;
}
close MDS;
print "chr\tpos\tN1\tN2\tFstNum\tFstDenom\tFst\tINV_0_allele\tINV_2_allele\tmds_coord";
foreach my $mds (sort keys %outlier_start){
  my $filename = "$fst_file_prefix.$mds.fst.txt.gz";
  open(IN, "zcat -c $filename |");
  while(<IN>){
    chomp;
    if ($_ =~ /^chr/){next;}
    my @a = split(/\t/,$_);
    my $chr = $a[0];
    my $pos = $a[1];
    my $fst = $a[6];
    if ($fst < $min_fst){next;}
    if ($outlier_start{$mds}{$chr}){
      foreach my $n (sort {$a <=> $b} keys %{$outlier_start{$mds}{$chr}}){
        if (($outlier_start{$mds}{$chr}{$n} <= $pos) and ($outlier_end{$mds}{$chr}{$n} >= $pos)){
          print "\n$_\t$mds";
        }
        
      }
    }
  }
  close IN;
}
