#!/usr/bin/perl

use strict;
use warnings;

open(FILE, $ARGV[0]);

my @FileContents = <FILE>;

close(FILE);

#print sprintf("%.2f", 13.0115) . "\n";

foreach my $Line (@FileContents) {

  my @split = split /\s+/, $Line;

  my @save;

  my $count = 0;

  foreach my $Element (@split) {

    unless($Element =~ /:/) { push @save, $Element; next }     
 
    my @spl = split /:/, $Element;

    $count++;
  
    if($count > 2) {  
  
      push @save, sprintf("%.3f", $spl[1]); 

    }

    else {

     push @save, $spl[1];
  
    }

  } 

  print join(" ", @save) . "\n";
   

}
