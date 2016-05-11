#!/usr/bin/env perl

use strict;
use warnings;

use Test::More tests => 1;
use lib $ENV{'MATCH'} . "/lib"; 
use Carp;

use MATCHBaseObject;
use BaseObject ':vars';

use MATCHParameters;
use MoleculeFileHandler;

BEGIN { use_ok('AtomTyper') };

#Deal With Initial Parameter Setup, This will be handled in main .pl script
my $DefaultParameters = MATCHParameters->New;

$DefaultParameters->Initiate;

$Parameters = $DefaultParameters;

#Grab Test Molecule

my $MoleculeFileHandler = MoleculeFileHandler->New("MoleculeParameterizerResources/1bjy_ligand_correct.pdb");

my $Structures = $MoleculeFileHandler->BuildObjectsFromFile;

my $TestMolecule = $Structures->[0]->getMolecule(0);

$TestMolecule->Initiate;

my $AtomTyper = AtomTyper->New;

$AtomTyper->Initiate;

#$AtomTyper->TypeAtomsInMolecule($TestMolecule);

chdir $ENV{'PerlChemistry'} . "/t/MoleculeFileHandler/TestMol2Files";

my @Mol2Files = <*>;

foreach my $Mol2File (@Mol2Files) {
	
	 $MoleculeFileHandler = MoleculeFileHandler->New($Mol2File);
	
	 $Structures = $MoleculeFileHandler->BuildObjectsFromFile;
	
	 print $Structures->[0]->getFileName . "\n";
	
	 my $Molecule = $Structures->[0]->getMolecule(0);
		
	 $Molecule->Initiate;
	
	 $AtomTyper->TypeAtomsInMolecule($Molecule);
	
}