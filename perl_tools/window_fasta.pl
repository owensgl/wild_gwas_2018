#!/bin/perl
use warnings;
use strict;
use POSIX;
use Parallel::ForkManager;
my $pm = new Parallel::ForkManager(15);
#This script breaks the vcf into windows and converts each window into a fasta
my $vcfgz =  $ARGV[0];
my $outprefix = $vcfgz;
my $vcfconverter = "/scratch/gowens/wild_gwas/wild_gwas_2018/perl_tools/vcf2fasta_basic.pl";
$outprefix =~ s/.vcf.gz//g;
my $window_size = 1000000;
my %chr;
open(IN, "gunzip -c $vcfgz |");
while(<IN>){
  chomp;
  if ($_ =~ m/##contig/){
    my $line = $_;
    $line =~ s/##contig=<ID=//g;
    my @a = split(/,/,$line);
    my $len = $a[1];
    $len =~ s/len=//g;
    $len =~ s/>//g;
    if ($a[0] =~ m/Chr00c/){next;}
    $chr{$a[0]} = $len;
  }elsif ($_ =~ m/#CHR/){
    goto PRINTOUT;
  }
}
PRINTOUT:
foreach my $chr (sort keys %chr){
  my $max_windows = floor($chr{$chr}/$window_size);
  foreach my $i (0..$max_windows){
    $pm->start and next;
    my $start = $i * $window_size;
    my $end = ($i+1) * $window_size;
    if ($end > $chr{$chr}){
      $end = $chr{$chr};
    }
   my $command = "bcftools view -r $chr:$start-$end $vcfgz -O v | perl $vcfconverter > $outprefix.$chr.$start.fasta";
   system($command);
   my $raxml_command = "/home/gowens/bin/standard-RAxML/raxmlHPC-PTHREADS-AVX -m ASC_GTRGAMMA --asc-corr lewis -s $outprefix.$chr.$start.fasta -T 2 -p 10920 -x 10920 -# 100 -n $outprefix.$chr.$start.1.tree";
   system($command);
   $pm->finish;
  }
}

$pm->wait_all_children;
