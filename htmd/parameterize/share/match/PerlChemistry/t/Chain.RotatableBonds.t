#!/usr/bin/perl -w

use strict;
use warnings;

use lib $ENV{'PerlChemistry'} . "/lib";

use Test::More tests => 1;
use BaseObject ':vars';
use Parameters;
use MoleculeFileHandler;

BEGIN { use_ok('Chain') };

my $DefaultParameters = Parameters->New;

$DefaultParameters->Initiate;

$Parameters = $DefaultParameters;

my $MoleculeFileHandler = MoleculeFileHandler->New($ARGV[0]);


my $Structures = $MoleculeFileHandler->BuildObjectsFromFile;

$Structures->[0]->Initiate;

my @RotBonds = $Structures->[0]->CalculateRotatableBonds;

print scalar(@RotBonds) . "\n";

foreach my $Bond (@RotBonds) {
	
	print $Bond->getPrimaryName . " " . $Bond->getSecondaryName . "\n";
	
}


