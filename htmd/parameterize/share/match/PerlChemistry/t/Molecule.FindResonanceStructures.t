#!/usr/bin/perl -w

use strict;
use warnings;
use Carp;

use lib $ENV{'PerlChemistry'} . "/lib";

use Molecule ':all';
use MoleculeFileHandler ':all';
use Parameters ':all';
use BaseObject ':vars';

#Deal With Initial Parameter Setup, This will be handled in main .pl script
my $DefaultParameters = Parameters->New;

$DefaultParameters->Initiate;

#Set Global Parameters
$Parameters = $DefaultParameters;

#Build all molecules from CGENFF topology file

my $MoleculeFileHandler = MoleculeFileHandler->New("../resources/top_all36_cgenff/top_all36_cgenff.rtf");

my $Structures = $MoleculeFileHandler->BuildObjectsFromFile;

my @Chains = grep { ! $_->getMolecule(0)->getIsPatch} @$Structures;

foreach my $Chain (@Chains) {
	
	#SetUp Atom LookUpTables
	#next if $Chain->getMolecule(0)->getName ne "DMSO";
	
	my $Atoms = $Chain->getMolecule(0)->getAtoms;
	
	foreach my $Atom (@$Atoms) {
				
		my $AtomLookUpTable = LookUpTable->New;
		
		$AtomLookUpTable->Initiate($Atom);
	  
	  $Atom->setLookUpTable($AtomLookUpTable);	  
		
	}
	
	$Chain->getMolecule(0)->NewEnsureCorrectBondTypes;		
	
	$Chain->getMolecule(0)->ReduceAtomicChargesToLowestPossible;	
	
	$Chain->getMolecule(0)->FindResonanceStructures;
	
	foreach my $Atom (@$Atoms) {
		
	#	print $Atom->getName . " " .  $Atom->getSumofBondTypes . " " . $Atom->getMissingCharge . " " . $Atom->getProtonationState . "\n";
		
		
	}
	
	
		
}