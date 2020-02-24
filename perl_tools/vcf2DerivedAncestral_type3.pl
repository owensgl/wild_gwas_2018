#!/bin/perl
use strict;
use warnings;

#This script takes a VCF, a list of outgroup samples, and outputs only sites where the outgroups are fixed for a single allele. It then polarizes all SNPs as A (ancestral) or # for derived. There can be multiple different derived alleles..
#In this version the outgroup is defined in a separate file

my $outfile = $ARGV[0]; #e.g. perennial_alleles.20190619.txt.gz
my $vcffile = $ARGV[1]; #The VCF to be processed
#First run through the VCF to load all sites that need processing
my $debug = 0; #Put to zero if you don't want to print the error messages
my $min_ancestral_sequenced = 1;
my %used_sites;
open(VCF, "gunzip -c $vcffile |");
while(<VCF>){
  chomp;
  my @a = split(/\t/,$_);
  if ($_ =~ m/^#/){
    next;
  }
  my $info = $a[7];
  my @infos = split(/\;/,$info);
  my $XRQ_info = $infos[$#infos];
  $XRQ_info =~ s/XRQ=//g;
  my $location = $XRQ_info;
  $used_sites{"$location"}++;
}
close VCF;

my %ancestral_state;
open(OUTGROUP, "gunzip -c $outfile |");
my $linecounter = 1;
while(<OUTGROUP>){
  chomp;
  if ($. == 1){
    next;
  }
  my @a = split(/\t/,$_);
  my $chr = $a[0];
  my $pos = $a[1];
  my $ref = $a[2];
  my $count = $a[3];
  if ($count < $min_ancestral_sequenced){next;}
  $linecounter++;
  if ($linecounter % 1000000 == 0){print STDERR "Processing outgroup file $chr $pos...\n";}
  unless ($used_sites{"$chr.$pos"}){next;}
  $ancestral_state{"$chr.$pos"} = $ref;
}
close OUTGROUP;

open(VCF2, "gunzip -c $vcffile |");
my %sample;
while(<VCF2>){
  chomp;
  my @a = split(/\t/,$_);
  if ($_ =~ m/^#CHR/){
    print "chr\tpos";
    foreach my $i (9..$#a){
      $sample{$i} = $a[$i];
      print "\t$a[$i]";
    }
    next;
  }
  if ($_ =~ m/^#/){next;}
  my $info = $a[7];
  my @infos = split(/\;/,$info);
  my $XRQ_info = $infos[$#infos];
  $XRQ_info =~ s/XRQ=//g;
  my $location = $XRQ_info;
  my $ref = $a[3];
  my $alts = $a[4];
  my @alts = split(/,/,$alts);
  my $found_ancestral;
  if ($ancestral_state{"$location"} eq $ref){
     $found_ancestral++;
  }
  foreach my $alt (@alts){
     if ($ancestral_state{"$location"} eq $alt){
       $found_ancestral++;
     }
  }
  unless ($found_ancestral){
    if ($debug){
      print STDERR "Ancestral allele not found $location\n";
    }
    next;
  }
  unless ($ancestral_state{"$location"}){
    if ($debug){
      print STDERR "Ancestral state not found\t$location\n";
    }
    next;
  }
  my $outgroup_allele;
  if ($ref eq $ancestral_state{"$location"}){
    $outgroup_allele = 0;
  }
  foreach my $i (0..$#alts){
    if ($alts[$i] eq $ancestral_state{"$location"}){
      $outgroup_allele = 1+$i;
    }
  }

  
  print "\n$a[0]\t$a[1]";
  #Now to figure out if the samples have derived or ancestral
  foreach my $i (9..$#a){
    if (($a[$i] eq '.') or ($a[$i] eq '.:0,0')){
      print "\tNN";
    }else{
      my @infos = split(/:/,$a[$i]);
      if ($infos[0] eq './.'){
        print "\tNN";
        next;
      }
      my @genos = split(/\//,$infos[0]);
      if ($genos[0] eq '.'){
	print "\tNN";
	next;
      }
      print "\t";
      foreach my $j (0..1){
	if ($genos[$j] eq $outgroup_allele){
	  print "A";
	}else{
	  print "$genos[$j]";
	}
      }
    }
  }
}
close VCF2;
