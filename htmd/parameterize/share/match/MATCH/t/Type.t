#!/usr/bin/env perl

use strict;
use warnings;

use Test::More tests => 3;
use lib '../lib'; 

BEGIN { use_ok('Type') };


my $Type = Type->New("CT1", 4, "C.4(H)");

$Type->Initiate;

my $Atoms = $Type->getAtoms;

foreach my $Atom (@$Atoms) {
	
	#print $Atom->getState . "\n";
	
}

my $Type2 = Type->New("CT1", 4, "!C.4(C.4%6,6(C.4%4(N.3%4,4,4)))");

$Type2->Initiate;

$Atoms = $Type2->getAtoms;

foreach my $Atom (@$Atoms) {
	
	#print $Atom->getState . "\n";
	#print join(" ", @{$Atom->getRings}) . "\n";
	
}


my $TestType = Type->New("CT2", 4, "C.4(H)(H)");

$TestType->Initiate;

is($Type->DoesLookUpTableMatch($TestType->getLookUpTable), 1, 'Types Match');


$TestType = Type->New("CG3C52", 4, "C.4%5([^C]%5)(H.1)");

$TestType->Initiate;

my $Type3 = Type->New("CT1", 4, "C.4%5(H)(H)");

$Type3->Initiate;

is($TestType->DoesLookUpTableMatch($Type3->getLookUpTable), 0, 'Types Do No Match');

$TestType = Type->New("HGR62", 1, "H.1(C.3%6(=C.3%6,*(F|CL|BR|I.1)))");

$TestType->Initiate;

foreach my $Atom (@{$TestType->getLookUpTable->getNodes}) {
	
  print $Atom->Stringify . "\n";

}

my $TestType2 = Type->New("Test", 1, "H.1(C.3%6(C.3%6(C.3%6)(F.1)))");

$TestType2->Initiate;

#is($TestType2->DoesLookUpTableMatch($TestType->getLookUpTable), 1, 'Types Match');

is($TestType->getLookUpTable->AreAllNodesSharedBetween([$TestType2->getLookUpTable]), 1, 'Types Match');


