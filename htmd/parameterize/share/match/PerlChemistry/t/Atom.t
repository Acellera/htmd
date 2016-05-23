#!/usr/bin/perl -w

use strict;
use warnings;

use lib $ENV{'PerlChemistry'} . "/lib";
use Test::More tests => 2;

subtest "Atom Package" => sub {

BEGIN { use_ok('Atom') };

use Atom;
use Parameters ':all';
use BaseObject ':vars';

#Deal With Initial Parameter Setup, This will be handled in main .pl script
my $DefaultParameters = Parameters->New;

$DefaultParameters->Initiate;

#Set Global Parameters
$Parameters = $DefaultParameters;

my $TestAtom = Atom->New;

$TestAtom->setName("N1");

is($TestAtom->getName, "N1", 'AutoLoad works');

my $TestAtom2 = Atom->New({ State => "C.3"});

is($TestAtom2->getState, "C.3", 'Fields correctly adds to Internal variables of Atom Object');

subtest "AtomInstance->Initate" => sub {
  plan tests => 3;

  my $TestAtom = Atom->New({ Type => "CG2O1" });

  $TestAtom->Initiate;

	is($TestAtom->getElement, "C", 'Element works');
	
	is($TestAtom->getElementWeight, "12.01100", 'ElementWeight works');
	
	is($TestAtom->getState, "C.0", 'State works');
		

  #is($TestAtom->getIsNormalBondingNum->{0}, -4, 'IsNormalBondingNum works');
  	
};

subtest "AtomInstance->Stringify" => sub {
    plan tests => 4;

    my $TestAtom2 = Atom->New({ State => "C.3"});

  	is($TestAtom2->Stringify, "C.3", 'Stringify Works');

		$TestAtom2->setNoRing(1);

		is($TestAtom2->Stringify, "!C.3", 'Stringify Works No Rings');

		my $TestRingAtom = Atom->New({ State => "C.3" });

    $TestRingAtom->setRingAtom(1);
		$TestRingAtom->setRings([6]);

		is($TestRingAtom->Stringify, "C.3%6", 'Stringify Works with Rings');

		$TestRingAtom->setRings([6,5]);

		is($TestRingAtom->Stringify, "C.3%6,5", 'Stringify Works with Rings');
};

subtest "AtomInstance->AddBond" => sub {
	
	plan tests => 13;
	
	use Bond;
	
	my $Atom1 = Atom->New({Name => "C1"});
	my $Atom2 = Atom->New({Name => "C2"});

	$Atom1->Initiate;
	$Atom2->Initiate;
	
	my $Bond = Bond->New($Atom1, $Atom2, 1);
	
	is($Atom1->getNumOfNotSolvedBonds, 1, 'getNumofNotSolvedBonds Works');
	
	is($Atom1->getSumofBondTypes, 1, 'getSumofBondTypes Works');
	
	is($Atom1->getState, "C.1", 'State works');
	
	is($Atom1->getProtonationState, -3, 'ProtonationState works');
		
	isnt($Atom1->IsBondedTo($Atom2), 0, '$Atom1 is Bonded to $Atom2');
	
	is($Atom1->IsBondedTo(Atom->New), 0, '$Atom1 not Bonded to random Atom');
		
	$Bond->setIncrement([.1, -.1]);
	
	$Bond->ApplyIncrement(-.1, .1);
	
	is($Atom1->getNumOfNotSolvedBonds, 0, 'getNumofNotSolvedBonds Works');
	
	$Bond->UnBondAtoms;
		
	is($Atom1->getNumOfNotSolvedBonds, 0, 'getNumofNotSolvedBonds Works');
	
	is($Atom1->getSumofBondTypes, 0, 'getSumofBondTypes Works');
	
	is($Atom1->getState, "C.0", 'State works');
	
	is($Atom1->getProtonationState, -4, 'ProtonationState works');
	
	is($Atom1->IsBondedTo($Atom2), 0, '$Atom1 is not Bonded to $Atom2');
	
	is($Atom1->getCharge, 0, 'Remove Bond cleanup worked');
		
};

subtest "AtomInstance->AmountofChargeToAllowedProtonationState" => sub {
	
	plan tests => 1;
	
	my $Atom1 = Atom->New({Name => "H1"});
	my $Atom2 = Atom->New({Name => "C2"});

	$Atom1->Initiate;
	$Atom2->Initiate;
		
	my $Bond = Bond->New($Atom1, $Atom2, 1);
		
	is($Atom1->AmountofChargeToAllowedProtonationState, 0, "Should be neutral");
	
	
};

#Test Copy, EXPAND! when Bond objects are working

my $Test = Atom->New({ State => "C.3", Rings => [5,6] });

my $CopiedTest = $Test->Copy;

is($Test->Stringify, $CopiedTest->Stringify, 'Test Stringify is the same as Copied');

};