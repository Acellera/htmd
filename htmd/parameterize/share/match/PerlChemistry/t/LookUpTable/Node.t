#!/usr/bin/perl -w

use strict;
use warnings;

use lib '../../lib';
use Test::More tests => 6;

BEGIN { use_ok('LookUpTable::Node') };

use LookUpTable::Node ':all';
use Atom;

my $Atom = Atom->New( { State => "C.3", Level => 0 });

my $Node = LookUpTable::Node->New($Atom);

my $Atom2 = Atom->New( { State => "C.4", Level => 0 });

my $Node2 = LookUpTable::Node->New($Atom2);

$Node->AddParent($Node2);

is($Node->getAtom, $Atom, 'Autoload test');

is($Node->getAtomState, 'C.3', 'Autoload test 2');

subtest "Test StringCompare" => sub {
	
	plan tests => 6;
		
	is($Node->StringCompare($Node2), 1, 'Compared Correctly');
	
	$Node2->setString("[^N].3");
	
	is($Node->StringCompare($Node2), 2, 'Not Nitrogen Works'); #may be interesting to implement this later
	
	is($Node2->StringCompare($Node), 2, 'Not Nitrogen Works'); #may be interesting to implement this later
	
	$Node2->setString("N.3");
	
	is($Node->StringCompare($Node2), 0, 'Incorrect Element');
	
	$Node2->setString("*");
	
	is($Node->StringCompare($Node2), 1, 'Wild Card works');
	
	$Node2->setString("C.4|3");
	
	is($Node->StringCompare($Node2), 2, 'Compared Correctly');	
	
};

subtest "Test FullStateCompare" => sub {
	
	plan tests => 18;
	
	$Node->setString("C.3"); $Node2->setString("C.3");
	
	$Node->setRings(["6N"]); $Node2->setRings(["6A"]);
	
	is($Node->FullStateCompare($Node2), 0, 'Aromatic Test');
	
	
	$Node->setRings(["6N"]); $Node2->setRings(["6N"]);
	
	is($Node->FullStateCompare($Node2), 3, 'Aromatic Test 2');
	
	
	$Node->setRings(["6N"]); $Node2->setRings(["6"]);
	
	is($Node->FullStateCompare($Node2), 3, 'Aromatic Test 3');
	
	
	$Node->setRings(["6"]); $Node2->setRings(["6A"]);
	
	is($Node->FullStateCompare($Node2), 3, 'Aromatic Test 4');
	
		
	$Node->setRings(["6"]); $Node2->setRings(["6"]);
	
	is($Node->FullStateCompare($Node2), 3, 'Compared Correctly');
	
	$Node2->setRings([5]);
	
	is($Node->FullStateCompare($Node2), 0, 'Difference Rings');
	
	$Node->setRings([6,6]); $Node2->setRings([6,5]);
	
	is($Node->FullStateCompare($Node2), 0, 'Difference Rings 2');
	
	$Node->setRings([6]);
	
	is($Node->FullStateCompare($Node2), 3, 'Same Ring');
	
	$Node->setRings([6,5]); $Node2->setRings([6]);
	
	is($Node->FullStateCompare($Node2), 0, 'Not Same Ring');
	
	$Node->setRings([5]); $Node2->setRings([5,5]);
	
	is($Node->FullStateCompare($Node2), 3, 'Not Same Ring');
	 
	$Node->setRings([]); $Node2->setRingAtom(0);
	$Node->setNoRing(1); $Node2->setRings([]);
	
	is($Node->FullStateCompare($Node2), 2, 'NoRing Set');
	
	$Node->setRings([]); $Node2->setRingAtom(1);
	$Node->setNoRing(1); $Node2->setRings([5]);
		
	is($Node->FullStateCompare($Node2), 0, 'NoRing Set 2');
	
	$Node->setRingAtom(1);
	
	is($Node->FullStateCompare($Node2), 0, 'NoRing Set 3');
	
	$Node->setRingAtom(1); 	$Node2->setNoRing(0);
	
	is($Node->FullStateCompare($Node2), 0, 'NoRing Set 4');
	
	$Node->setString("N.2"); $Node2->setString("N.3");
	$Node->setRingAtom(1); 	$Node2->setNoRing(1);
	$Node->setRings([5]); 	$Node2->setRings([5]);
	
	is($Node->FullStateCompare($Node2), 0, 'Test!');
	
	$Node->setString("N.3"); $Node2->setString("N.3");
	$Node->setRingAtom(1); 	 $Node2->setRingAtom(1); 	 
	$Node->setRings([5,"*"]);  $Node2->setRings([5,6]);
	$Node->setNoRing(0);     $Node2->setNoRing(0);
	
	is($Node->FullStateCompare($Node2), 4, 'Test!');
	
	$Node->setString("N.3"); $Node2->setString("N.3");
	$Node->setRingAtom(0); 	 $Node2->setRingAtom(0); 	 
	$Node->setRings([]);  $Node2->setRings([]);
	$Node->setNoRing(0);     $Node2->setNoRing(0);
	$Node->setMissingCharge(1); $Node2->setMissingCharge(0);	
	
	is($Node->FullStateCompare($Node2), 0, 'Charge Test');
	
	$Node->setMissingCharge(1); $Node2->setMissingCharge(.33);	
	
	is($Node->FullStateCompare($Node2), 2, 'Charge Test');
	
};

#subtest "Test StringMatchesTo" => sub {
	
	#plan tests => 1;
	
	#$Node2->setString("C.3");
	
	#is($Node->StringMatchesTo($Node2), 1, 'Matches');
	
	
	
#};

$Node->setRingAtom(1);
$Node2->setRingAtom(1);
$Node2->setNoRing(0);


$Node->setRings([5,6]); $Node2->setRings([6]);

print $Node->FullMatchesTo($Node2) . "\n";


