#!/usr/bin/env perl

use strict;
use warnings;

use Test::More tests => 1;

use lib '../../lib'; 

BEGIN { use_ok('Type') };

use Type ':test';
use Atom;

subtest "Type->GetAtomInformationFromString" => sub {
	
	plan tests => 5;
	
	my $Atom = Atom->New;
	
	GetAtomInformationFromString($Atom, "C.3");
	
	is($Atom->getState, "C.3", 'State recored correctly');
	
	is($Atom->getNumOfBonds, 3, 'Num Of Bonds recored correctly');
	
	is_deeply($Atom->getRings, [], 'Rings correct');
	
	my $Atom2 = Atom->New;
	
	GetAtomInformationFromString($Atom2, "N.3%5,6,7");
	
	is($Atom2->getState, "N.3", 'State recored correctly');
	
	is_deeply($Atom2->getRings, [5,6,7], 'Rings correct');
	
	
};