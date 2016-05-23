#!/usr/bin/perl -w

use strict;
use warnings;
use Carp;

use lib '../lib';
use Test::More tests => 4;

BEGIN { use_ok('Molecule') };

use Molecule ':all';
use MoleculeFileHandler ':all';
use Parameters ':all';
use BaseObject ':vars';

use TestingFunctions ':func';

#Deal With Initial Parameter Setup, This will be handled in main .pl script
my $DefaultParameters = Parameters->New;

$DefaultParameters->Initiate;

#Set Global Parameters
$Parameters = $DefaultParameters;

my $MoleculeFileHandler = MoleculeFileHandler->New("MoleculeFileHandler/TestPDBFiles/1byj_ligand_correct.pdb");

my $Structures = $MoleculeFileHandler->BuildObjectsFromFile;

my ($Molecule) = @{$Structures->[0]->getMolecules};

subtest "Copy" => sub {
	
	plan tests => 3;
	
	my $CopiedMolecule = $Molecule->Copy;
	
	is($Molecule->getName, $CopiedMolecule->getName, 'Named Match');
	
	my $Atoms = $Molecule->getAtoms;
	my $CopiedAtoms = $CopiedMolecule->getAtoms;
	
  my %AtomNames = map { $_->getName => $_ } @$Atoms;
  my %CopiedAtomNames = map { $_->getName => $_ } @$CopiedAtoms;

  my $AtomsMatch = 1;
  my $BondsMatch = 1;

  while ( my ($AtomName, $Atom) = each %AtomNames ) {
		
		my $CopiedAtom =  $CopiedAtomNames{$AtomName};
		
	  unless(is_atom_copied($Atom,$CopiedAtom)) {
		
		  $AtomsMatch = 0; last;
		
	  }
	
	  my $Bonds = $Atom->getBonds;
	  my $CopiedBonds = $CopiedAtom->getBonds;
	
	  foreach my $Bond (@$Bonds) {
		
		  my $Found = 0;
		
		  foreach my $CopiedBond (@$CopiedBonds) {
			
			  next if $Bond->getBondingPartner($Atom)->getName ne $CopiedBond->getBondingPartner($CopiedAtom)->getName;
			
			  $Found = 1; 
			
		    if($Bond->getType != $CopiedBond->getType) {
			
			    $BondsMatch = 0; last;
			
		    }
		
		    if(scalar(@{$Bond->getIncrement}) > 0 && scalar(@{$CopiedBond->getIncrement}) ) {
						
			     $BondsMatch = 0 if $Bond->getIncrement->[0] != $CopiedBond->getIncrement->[0]; last;
			
		    }
			
		  }
		
		  $BondsMatch = 0 if !$Found;
		
	  }
	
  }

  is($AtomsMatch, 1, 'All Atom Properties Match Up');
  is($BondsMatch, 1, 'All Bond Properties Match Up');
  

	
	
};






