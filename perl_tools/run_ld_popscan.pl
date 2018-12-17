#!/bin/perl
use warnings;
use strict;
use POSIX;
use Parallel::ForkManager;

#This script runs vcftools for the whole dataset and then divides by population

my $popfile = $ARGV[1];
my $gzvcf = $ARGV[0];
my $prefix = $gzvcf;
$prefix =~ s/.vcf.gz//g;
my $processing_script = "/home/owens/bin/reformat/emerald2windowldcounts.pl";

my $ncores = 8;
my $pm = new Parallel::ForkManager($ncores);

my $chr_prefix = "Ha412TMPChr";
my $chr_n = 17;
my $mac = 4;
my %pops;
my %fh;
open POP, $popfile;
my $file = "$popfile.all.tmp";
open $fh{"all"},'>', $file or die "Can't open the output file: $!";

while(<POP>){
  chomp;
  my @a = split(/\t/,$_);
  my $sample = $a[0];
  my $pop = $a[1];
  $pops{$pop}++;
  print { $fh{"all"} } "$sample\n";
  if ($pop =~ m/_N/){next;} #Skip sketchy populations;
  unless ($fh{$pop}){
   $file = "$popfile.$pop.tmp";
   open $fh{$pop},'>', $file or die "Can't open the output file: $!";
  }
  print { $fh{$pop} } "$sample\n";
}
foreach my $pop (sort keys %pops){
  foreach my $chr (1..$chr_n){
    my $padded_chr =  sprintf("%02d", $chr);
    my $full_chr = $chr_prefix . $padded_chr;
    $pm->start and next;
    my $system = "/home/owens/bin/emeraLD/bin/emeraLD -i $gzvcf  --stdout --window 1000000000 --region $full_chr:1-1000000000 --include $popfile.$pop.tmp  --mac $mac";
    $system .= " | perl $processing_script > $prefix.$pop.$full_chr.window.ld";
    system($system);
    $pm->finish;
  }
}


$pm->wait_all_children;
system ("rm $popfile*tmp");
