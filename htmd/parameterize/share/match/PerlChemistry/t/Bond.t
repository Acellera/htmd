#!/usr/bin/perl -w

use strict;
use warnings;

use lib $ENV{'PerlChemistry'} . "/lib";
use Test::More tests => 2;

use BaseObject ':vars';

use Parameters;

my $DefaultParameters = Parameters->New;

$DefaultParameters->Initiate;

$Parameters = $DefaultParameters;

subtest "Bond Package" => sub {

BEGIN { use_ok('Bond') };

use Bond ':all';

my $Atom1 = Atom->New({Name => "C1"});
my $Atom2 = Atom->New({Name => "C2"});

$Atom1->Initiate;
$Atom2->Initiate;

my $Bond = Bond->New($Atom1, $Atom2, 1);

is($Bond->getPrimaryName, "C1", 'Test Bond Autoloading function');

is($Bond->IsIncrementSolved, 0, 'Increment Should Not be solved');

subtest "BondInstance->ApplyIncrement" => sub {
 
  plan tests => 4;
 
  $Bond->setIncrement([.1, -.1]);

  is($Bond->IsIncrementSolved, 1, 'Increment is Solved');

  $Bond->ApplyIncrement(-.1, .1);

  is($Atom1->getCharge, -.1, 'Application of Charge worked');

  $Bond->UnSetIncrement;

  is($Atom1->getCharge, 0, 'Removal of Charge worked');

  is($Bond->IsIncrementSolved, 0, 'Increment is Not Solved');
  
};

subtest "BondInstance->SwapPrimaryAtom" => sub {

  plan tests => 2;
  
  $Bond->setIncrement([.1, -.1]);
  
  $Bond->SwapPrimaryAtom;

  is($Bond->getPrimary, $Atom2, 'Atom swap worked');

  is_deeply($Bond->getIncrement, [-0.1, 0.1], 'Increment swap worked');

  $Bond->UnSetIncrement;

};

is_deeply($Bond->getBondingPartner($Atom1), $Atom2, 'getBondingPartner works');

$Bond->MakeThisAtomPrimary($Atom1);

is_deeply($Bond->getPrimary, $Atom1, 'Testing MakeThisAtomPrimary 1');

$Bond->MakeThisAtomPrimary($Atom2);

is_deeply($Bond->getPrimary, $Atom2, 'Testing MakeThisAtomPrimary 2');

};

