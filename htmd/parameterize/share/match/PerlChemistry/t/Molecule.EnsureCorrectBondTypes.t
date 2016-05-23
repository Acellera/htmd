#!/usr/bin/perl -w

use strict;
use warnings;
use Carp;

use lib '../lib';
use Test::More tests => 1;

BEGIN { use_ok('Molecule') };

use Molecule ':all';
use MoleculeFileHandler ':all';
use Parameters ':all';
use BaseObject ':vars';

#Deal With Initial Parameter Setup, This will be handled in main .pl script
my $DefaultParameters = Parameters->New;

$DefaultParameters->Initiate;

#Set Global Parameters
$Parameters = $DefaultParameters;

#Build all molecules from CGENFF topology file, look at all molecules that are residues and find all Carbons that are in rings. This is easy due to 
#the naming convention, all Ring atoms have R in their Type name. This will print out the name of all atoms that have R in their name but do not belong
#to any rings

my $start = time();

my $MoleculeFileHandler = MoleculeFileHandler->New($ENV{'MATCH'} . "/resources/top_all36_cgenff_new/top_all36_cgenff.rtf");
#my $MoleculeFileHandler = MoleculeFileHandler->New("../../MATCH/resources/top_all36_lipid/top_all36_lipid.rtf");


my $Structures = $MoleculeFileHandler->BuildObjectsFromFile;

my @Chains = grep { ! $_->getMolecule(0)->getIsPatch} @$Structures;

foreach my $Chain (@Chains) {
	
	#next if $Chain->getMolecule(0)->getName ne "NADP";
			
	print $Chain->getMolecule(0)->getName .  "\n";
	
	$Chain->getMolecule(0)->CleanEnsureCorrectBondTypes;
	#$Chain->getMolecule(0)->NewEnsureCorrectBondTypes;
	
	foreach my $Atom (@{$Chain->getMolecule(0)->getAtoms}) {
		
		#print $Atom->getName . " " . $Atom->getProtonationState . " " . $Atom->getSumofBondTypes . "\n";
		#$TotalCharge += $Atom->getProtonationState;	
		
	}
	
	
	#$Chain->getMolecule(0)->Initiate;
	
	#print "Reduce Charge!\n";
	
	$Chain->getMolecule(0)->NewReduceAtomicChargesToLowestPossible;		
		
	my $TotalCharge  = 0;
	
	foreach my $Atom (@{$Chain->getMolecule(0)->getAtoms}) {
		
		#print $Atom->getName . " " . $Atom->getProtonationState . " " . $Atom->getSumofBondTypes . "\n";
		$TotalCharge += $Atom->getProtonationState;	
		
	}
		
	if($TotalCharge != $Chain->getMolecule(0)->getCharge) {print $Chain->getMolecule(0)->getName . " " .  $TotalCharge . " " . $Chain->getMolecule(0)->getCharge . "\n" }
	
	#last;
		
}

my $end = time();

print "time : " . ($end - $start) . "\n";