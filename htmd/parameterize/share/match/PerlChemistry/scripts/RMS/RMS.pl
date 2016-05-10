#!/usr/bin/perl -w

use strict;
use warnings;
use Carp;
use Cwd qw(abs_path);

use lib '../../lib';
use MoleculeFileHandler::PdbFileHandler ':all';
use Parameters ':all';
use BaseObject ':vars';

#Deal With Initial Parameter Setup, This will be handled in main .pl script
my $DefaultParameters = Parameters->New;

my ($Path) = (abs_path($0) =~ /(\S+)\/\S+$/);

my ($DefaultParametersPath) = ($Path =~ /(\S+)\/\S+$/);

$DefaultParametersPath = "/Users/skullnite/Dropbox/Public/Work/MATCH_RELEASE/PerlChemistry/resources/DefaultParameters.par";

$DefaultParameters->Initiate($Path, $DefaultParametersPath);

#Set Global Parameters
$Parameters = $DefaultParameters;

my $MoleculeFileHandler = MoleculeFileHandler->New("1eht_ligand_correct.pdb");

my $Structures = $MoleculeFileHandler->BuildObjectsFromFile;

my $CorrectNamed = $Structures->[0];

RemoveHydrogens($CorrectNamed->getMolecule(0));

$MoleculeFileHandler = MoleculeFileHandler->New("1eht_ligand_icm.pdb");

$Structures = $MoleculeFileHandler->BuildObjectsFromFile;

my $IncorrectNamed = $Structures->[0];

RemoveHydrogens($IncorrectNamed->getMolecule(0));

$CorrectNamed->Initiate; $IncorrectNamed->Initiate;

my $CorrectNamedAtoms = $CorrectNamed->getMolecule(0)->getAtoms;

my $IncorrectNamedAtoms = $IncorrectNamed->getMolecule(0)->getAtoms;

my %Used;

my $RMSD = 0;

foreach my $IncorrectAtom (@$IncorrectNamedAtoms) {
			
	foreach my $CorrectAtom (@$CorrectNamedAtoms) {
		
		if(! exists $Used{$CorrectAtom->getName} && $CorrectAtom->DoesLookUpTableMatch($IncorrectAtom->getLookUpTable)) {
		
			$Used{$CorrectAtom->getName}++;
			
      my $DiffInDistance = 0;

      my $IncorrectAtomCoor = $IncorrectAtom->getCartesianCoordinates;
      my $CorrectAtomCoor = $CorrectAtom->getCartesianCoordinates;

      $DiffInDistance += ($IncorrectAtomCoor->[$_] - $CorrectAtomCoor->[$_]) ** 2 foreach (0,1,2);

      $RMSD += $DiffInDistance;
			
		}
		
	}	

}

$RMSD /= @$IncorrectNamedAtoms;

$RMSD = sqrt($RMSD);

print $RMSD . "\n";

#$IncorrectNamed->WriteToPdbFile("1eht_ligand_icm_renamed");


sub RemoveHydrogens { 
	
	my $Molecule = shift;
	
	my $Atoms = $Molecule->getAtoms;
	
	my $AtomCounter = 0;
	
	while ($AtomCounter < @$Atoms) {

		if($Atoms->[$AtomCounter]->getElement eq 'H') {

			foreach my $Bond (@{$Atoms->[$AtomCounter]->getBonds}) {

		    $Bond->UnBondAtoms;			
			}

			splice(@$Atoms, $AtomCounter, 1);

		}

		else { $AtomCounter++ }

	}
	
	
}