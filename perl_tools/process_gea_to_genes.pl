#!/bin/perl
#This takes the gff from Ha412HO, pulls out regions (plus surround)
#Then it takes the GEA file and pulls out the stats for each GEA variable
use strict;
use warnings;
use POSIX;
my $gff_file = $ARGV[0];
my $gea_file = $ARGV[1];
my $column_file = $ARGV[2]; #e.g. column.names.annuus.snps.varout.txt
my $surrounding_bases = 2000;
my $info_columns = 3; #2 for soil, 3 for climate
open(GFF, "gunzip -c $gff_file |");

my %gene_start;
my %gene_chr;
my %gene_end;
my %gene_window;
my $window_size = 100000; #For speed up purposes;
while(<GFF>){
  chomp;
  if ($_ =~ m/^#/){next;}
  my @a = split(/\t/,$_);
  if ($a[2] ne "gene"){next;}
  my $chr = $a[0];
  my $start = $a[3] - $surrounding_bases;
  my $end = $a[4] + $surrounding_bases;
  my @infos = split(/;/,$a[8]);
  my $name = $infos[0];
  $name =~ s/ID=gene://g;
  my $window_1 = floor($start/$window_size);
  my $window_2 = floor($end/$window_size);
  $gene_start{$name} = $start;
  $gene_end{$name} = $end;
  $gene_chr{$name} = $chr;
  push( @{ $gene_window {"$chr.$window_1" } }, "$name");
  if ($window_1 != $window_2) {
    push( @{ $gene_window {"$chr.$window_2" } }, "$name");
  }
#  print STDERR "Processing $name in $chr.$window_1\n";
}
close GFF;
my @columns;
open (COL, $column_file);
while(<COL>){
  chomp;
  @columns = split(/\t/,$_);
}
close COL;
open(GEA, "gunzip -c $gea_file |");
my %gea_total_sites;
my %gea_total_hits;
my %gea_top_hit;
my %gea_second_hit;
my $old_chr = 0;
while(<GEA>){
  chomp;
  my @a = split(' ',$_);
  my @infos = split(/__/,$a[0]);
  my $chr = $infos[0];
  if ($chr ne $old_chr){
    print STDERR "Processing $chr...\n";
    $old_chr = $chr;
  }
  my $pos = $infos[1];
  my $window = floor($pos/$window_size);
  #Find all possible genes in that window.
  unless($gene_window{"$chr.$window"}){next;} #Skip if no genes around there.
  my @possible_genes = @{ $gene_window{"$chr.$window"} };
  foreach my $gene (@possible_genes){
    if (($gene_start{$gene} <= $pos) and ($gene_end{$gene} >= $pos)){
      foreach my $i ($info_columns..$#a){
        $gea_total_sites{$i}{$gene}++;
        if ($a[$i] >= 10){
          $gea_total_hits{$i}{$gene}++;
        }
        unless($gea_top_hit{$i}{$gene}){
          $gea_top_hit{$i}{$gene} = $a[$i];
          $gea_second_hit{$i}{$gene} = 0;
        }elsif ($a[$i] > $gea_top_hit{$i}{$gene}){
          $gea_second_hit{$i}{$gene} = $gea_top_hit{$i}{$gene};
          $gea_top_hit{$i}{$gene} = $a[$i];
        }elsif ($a[$i] > $gea_second_hit{$i}{$gene}){
          $gea_second_hit{$i}{$gene} = $a[$i];
	}
      }
      next;
    }
  }
}
print "variable\tgene\tchr\tstart\tend";
print "\ttotal_sites";
print "\ttotal_hits";
print "\ttop_hit";
print "\tsecond_hit";
#Now print out all the results;
foreach my $i ($info_columns..$#columns){
  foreach my $gene (sort keys %{$gea_total_sites{3}}){
    print "\n$columns[$i]\t$gene\t$gene_chr{$gene}\t$gene_start{$gene}\t$gene_end{$gene}";
    print "\t$gea_total_sites{$i}{$gene}";
    unless($gea_total_hits{$i}{$gene}){
      $gea_total_hits{$i}{$gene} = 0;
    }
    print "\t$gea_total_hits{$i}{$gene}";
    print "\t$gea_top_hit{$i}{$gene}";
    print "\t$gea_second_hit{$i}{$gene}";
  }
}


