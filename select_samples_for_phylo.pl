#!/bin/perl
use warnings;
use strict;

#This script takes a list of samples and populations, as well as the list of genotype calls, and picks samples for phylogenetics of inversions.
#It picks one sample per genotype per population.

my $info_file = $ARGV[0];
my $genotype_file = $ARGV[1];
my $target_species = $ARGV[2];
my $samples_per_species = 5;

open INFO, $info_file; 
my %species;
my %pop;
my %species_hash;
my %pop_hash;
while(<INFO>){
  chomp;
  if ($. == 1){next;}
  my @a = split(/\t/,$_);
  my $name = $a[2];
  my $pop = $a[4];
  my $species = $a[3];
  my $type = $a[5];
  my $project = $a[6];
  #Skipping bad samples;
  if ($project eq "pangenome"){next;}
  if ($type ne "wild"){next;}
  if ($pop =~ m/_N/){next;}
  if ($pop =~ m/MK/){next;}
  if ($pop =~ m/IA/){next;}
  if ($pop =~ m/KS/){next;}
  if ($pop =~ m/Man/){next;}
  if ($pop =~ m/Pioneer/){next;}
  if ($pop =~ m/ND/){next;}
  if ($pop =~ m/SD/){next;}
  if ($pop =~ m/SK/){next;}
  $pop{$name} = $pop;
  $species{$name} = $species;
  $species_hash{$species}++;
  $pop_hash{$pop}++;
}

close INFO;

open GENO, $genotype_file;
my %mds_coords;
my %data;
while(<GENO>){
  chomp;
  if ($. == 1){next;}
  my @a = split(/\t/,$_);
  my $mds = $a[0];
  my $name = $a[1];
  my $genotype = $a[3];
  $mds_coords{$mds}++;
  $data{$mds}{$name} = $genotype;
}

#Loop through the mds_coordinates and output the samples I want to use.

foreach my $mds (sort keys %mds_coords){
  my $file_1 = "$target_species.$mds.phylosamples.names.txt";
  my $file_2 = "$target_species.$mds.phylosamples.info.txt";
  open my $fh_1, '>', $file_1;
  open my $fh_2, '>', $file_2;
  #Print target species samples;
  foreach my $pop (sort keys %pop_hash){
    foreach my $name ( keys %pop){
      if (($pop{$name} eq $pop) and (exists($data{$mds}{$name}))){
        if ($data{$mds}{$name} eq 0){
          print $fh_1 "$name\n";
          print $fh_2 "$name\t$species{$name}.0\t$pop\n";
          goto TRY2;
        }
      }
    }
    TRY2:
    foreach my $name ( keys %pop){
      if (($pop{$name} eq $pop) and ($data{$mds}{$name})){
        if ($data{$mds}{$name} eq 2){
          print $fh_1 "$name\n";
          print $fh_2 "$name\t$species{$name}.2\t$pop\n";
          goto DONE2;
        }
      }
    }
    DONE2:
    
  }
  #Print out other species;
  foreach my $species ( keys %species_hash){
    if ($species eq $target_species){next;}
    my $species_count;
    foreach my $name  (keys %species){
      if ($species{$name} eq $species){
        print $fh_1 "$name\n";
        print $fh_2 "$name\t$species\tNA\n";
        $species_count++;
        if ($species_count == $samples_per_species){
          goto NEXTSPECIES;
        }
      }
    }
    NEXTSPECIES:
  }
}


