#!/bin/perl
use strict;
use warnings;

my $prev_chr;
my $chunk_start;
my $prev_position;
my $max_distance = 500000;
my $first_line;
while(<STDIN>){
  chomp;
  my @a = split(/\t/,$_);
  unless ($a[0] =~ m/HanXRQChr/){next;}
  unless($prev_chr){
    $prev_chr = $a[0];
    $chunk_start = $a[1];
  }else{
    if ($prev_chr ne $a[0]){
      if ($first_line){
        print "\n$prev_chr\t$chunk_start\t$prev_position";
      }else{
        print "$prev_chr\t$chunk_start\t$prev_position";
	$first_line++;
      }
      $prev_chr = $a[0];
      $chunk_start = $a[1];
    }elsif (($a[1] - $prev_position > $max_distance)){
      if ($first_line){
        print "\n$prev_chr\t$chunk_start\t$prev_position";
      }else{
        print "$prev_chr\t$chunk_start\t$prev_position";
        $first_line++;
      }
      $prev_chr = $a[0];
      $chunk_start = $a[1];
    }
  }
  $prev_position=$a[1];
}
if ($first_line){
        print "\n$prev_chr\t$chunk_start\t$prev_position";
}else{
        print "$prev_chr\t$chunk_start\t$prev_position";
        $first_line++;
}

