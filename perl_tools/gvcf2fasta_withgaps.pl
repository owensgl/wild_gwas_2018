#!/bin/perl
use warnings;
use strict;
#GVCF into fasta
my $min_depth = 2;
my $sample;
my $gaps = 0;
my $header;
my $seq;
my $position = 1;
while(<STDIN>){
  chomp;
  if ($_ =~ m/^##/){next;}
  if ($_ =~ m/^#/){
    my @a = split(/\t/,$_);
    $sample = $a[9];
    $header .= ">$sample";
    next;
  }
  my @a = split(/\t/,$_);
  my $ref = $a[3];
  my $alt = $a[4];
  my @alts = split(/,/,$alt);
  my $genotype;
 #If we previously selected an alt allele shorter than the ref allele, then we have to introduce gaps.
  if ($gaps > 0){
    $genotype = "-";
    $gaps--;
    goto PRINTOUT;
  }
  my @field_names = split(/:/,$a[8]);
  #Skip weird rows
  unless (exists($field_names[2])){
    $genotype = "N";
    goto PRINTOUT;
  }
  if ($field_names[2] ne "DP"){
    $genotype = "N";
    goto PRINTOUT;
  }
  my $info = $a[9];
  my @infos = split(/:/,$info);
  my $depth = $infos[3];
  my $call = $infos[0];
  my @calls;
  if ($call =~ /\//){
    @calls = split(/\//,$call);
  }elsif ($call =~ /\|/){
    @calls = split(/\|/,$call);
  }
  my $rand_call = int(rand(2));
  my $called_value = $calls[$rand_call];
  if ($depth < $min_depth){
    $genotype = "N";
    goto PRINTOUT;
  }
  if (length($ref) == 1){
    if ($called_value == 0){
      $genotype = $ref;
    }else{
      $genotype = $alts[($called_value - 1)];
      if (length($genotype) > length($ref)){
        my $extra_gap = length($genotype) - length($ref);
        $header .= ":$position,$extra_gap";
      }
    }
  }else{
    if ($called_value == 0){
      my @bases = split(//,$ref);
      $genotype = $bases[0]; #Only print the first base of the ref, because the ref is represented by further lines
    }else{
      $genotype = $alts[($called_value - 1)];
      if (length($genotype) < length($ref)){
	$gaps = (length($ref) - length($genotype));
      }
      if (length($genotype) > length($ref)){
	my $extra_gap = length($genotype) - length($ref);
        $header .= ":$position,$extra_gap";
      }
    }
  }
  PRINTOUT:
  $seq .= "$genotype";
  $position++;
}
print "$header\n$seq";
