#!/bin/perl
use warnings;
use strict;
use Parallel::ForkManager;

#This script takes a list of sites, uses samtools to pull out the XRQ region, then blasts the HA412 genome and pulls out hits that are high enough.


my $ha412_ref = "/home/owens/ref/ha412_sunflower_26Apr2018_CHus7_top17scaf.fa";
my $xrq_ref = "/home/owens/ref/HanXRQr1.0-20151230.fa";

my $bases_surrounding=100; #Number of bases before and after the target site for blasting.
my $min_bit = 250;
my $counter = 0;
my $ncores = 10;
my $pm = new Parallel::ForkManager($ncores);


while(<STDIN>){
  chomp;
  my @a = split(/\t/,$_);
  $counter++;
  $pm->start and next;
  my $chr = $a[0];
  if ($chr eq "chr"){next;}
  my $pos = $a[1];
  my $start = $pos - $bases_surrounding;
  my $end = $pos + $bases_surrounding;

  system("samtools faidx $xrq_ref $chr:$start-$end > tmp.$counter.fa");
  my $blast_counter = 1;
  open CMD,'-|',"blastn -db $ha412_ref -query tmp.$counter.fa -outfmt 6" or die $@;
  my $line;
  while (defined($line=<CMD>)) {
      chomp($line);
      my @b = split(/\t/,$line);
      my $blast_chr = $b[1];
      my $blast_start_xrq = $b[6];
      my $blast_end_xrq = $b[7];
      my $blast_start_ha412 = $b[8];
      my $blast_end_ha412 = $b[9];
      my $bit_score = $b[11];
      if ($bit_score < $min_bit){next;}
      my $snp_position;
      if ($blast_end_ha412 > $blast_start_ha412){
        #it's going forward
        $snp_position = $blast_start_ha412 + (101 - $blast_start_xrq);
      }else{
        $snp_position = $blast_start_ha412 - (101 - $blast_start_xrq);
      }
      print "$chr\t$pos\t$blast_chr\t$snp_position\t$bit_score\n";
      goto NEXTLINE;
  }
  NEXTLINE:
  close CMD;
  system("rm tmp.$counter.fa");
  $pm->finish
}
$pm->wait_all_children;
