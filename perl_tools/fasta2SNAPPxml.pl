#!/bin/perl
#This script takes a Fasta file and converts it to a BEASTv2 xml

use warnings;
use strict;

my $fasta = $ARGV[0];
my $pop_info = $ARGV[1];
my $chainlength="25000";
my $storeEvery_chain="100";
my $storeEvery_state="1000";


my $dataset=$fasta;
$dataset=~ s/.fasta//g;
$dataset=~ s/.fa//g;

my $header = "SNAPP_header.txt";
my $map = "SNAPP_map.txt";
my $tail = "SNAPP_tail.txt";
open POP, $pop_info;
my %pop;
my %pop_list;

while(<POP>){
  chomp;
  my @a = split(/\t/,$_);
  $pop{$a[0]} = $a[1];
  $pop_list{$a[1]}++;
}
close POP;

open FASTA, $fasta;
my %seq;
my $current_sample;
while(<FASTA>){
  chomp;
  if ($_ =~ m/^>/){
    my $sample = $_;
    $sample =~ s/\>//g;
    $current_sample = $sample;
  }else{
    $seq{$current_sample} = $_;
  }
}
#Print out the pre-determined header
my $header_text = `cat $header`;
print "$header_text";

#Print out the sequences for each population
print "<data id=\"$dataset\" name=\"rawdata\">\n";
foreach my $sample (sort keys %pop){
  print "<sequence id=\"seq_$sample\" taxon=\"$sample\" totalcount=\"4\" value=\"$seq{$sample}\"/>\n"
}
print "</data>\n";
#Print out math distributions
my $map_text = `cat $map`;
print "\n$map_text\n";

#Print out taxon information
print "<run id=\"mcmc\" spec=\"MCMC\" chainLength=\"$chainlength\" storeEvery=\"$storeEvery_chain\">\n";
print "<state id=\"state\" storeEvery=\"$storeEvery_state\">\n";
print "<stateNode id=\"Tree.$dataset\" spec=\"beast.util.ClusterTree\" clusterType=\"upgma\" nodetype=\"snap.NodeData\">\n";
print "<taxa id=\"snap.$dataset\" spec=\"snap.Data\" dataType=\"integerdata\">\n";
print "<rawdata idref=\"$dataset\"/>\n";
foreach my $pop (sort keys %pop_list){
  print "<taxonset id=\"$pop\" spec=\"TaxonSet\">\n";
  foreach my $sample (sort keys %pop){
    if ($pop eq $pop{$sample}){
      print "<taxon id=\"$sample\" spec=\"Taxon\"/>\n";
    }
  }
  print "</taxonset>\n";
}
print "</taxa>\n";

my $tail_text = `cat $tail | sed s/FILENAME/$dataset/g`;

print "$tail_text";
