#!/usr/bin/perl

use strict;
use warnings;

my @FF = qw(top_all22_prot top_all27_na top_all35_carb top_all35_ethers top_all36_cgenff top_all36_lipid);

system("rm molecule.txt");

foreach my $FF1 (@FF) {

  foreach my $FF2 (@FF) {

    system("./TestIncrements.t $FF1 $FF2");

  }

}
