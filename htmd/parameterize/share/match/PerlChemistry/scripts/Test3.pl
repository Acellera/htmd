#!/usr/bin/perl

use strict;
use warnings;

use lib $ENV{'PerlChemistry'} . "/lib";

use Parameters ':all';
use BaseObject ':vars';
use MoleculeFileHandler ':all';

use ChemicalSpaceMatch;

#Deal With Initial Parameter Setup, This will be handled in main .pl script
my $DefaultParameters = Parameters->New;

$DefaultParameters->Initiate;

#Set Global Parameters
$Parameters = $DefaultParameters;

$Parameters->setCheckElementBondingNumber(0);

my $PDBFileHandler = MoleculeFileHandler->New($ARGV[0]);

my $PDBStructures = $PDBFileHandler->BuildObjectsFromFile;

my $LigandChain;

foreach my $Structure(@$PDBStructures) {
	
	if($Structure->getSegmentId eq "LIG0") {
		
		$LigandChain = $Structure; last;
		
	}
	
}

my $Molecule = $PDBStructures->[0]->getMolecule(0);

$Molecule->Initiate();

my $Atoms = $Molecule->getAtoms;


my @Associations;
my %AssociationsDependencies;
my %HashedAssociations;

my %AcceptedAssociationsSwaps;
my %ContainedInADependecy;

my %AssociationLevel;

#Setup ChemicalSpaceMatches

my @HeavyAtoms =  @$Atoms;

my @ChemicalSpaceMatches;

foreach my $i (0 .. @HeavyAtoms - 1) {
	
	$HeavyAtoms[$i]->getLookUpTable->SetAtomLevelsUsingNodeLevels;
		
	foreach my $j ($i+1 .. @HeavyAtoms - 1) {
		
		if($HeavyAtoms[$i]->DoesLookUpTableMatch($HeavyAtoms[$j]->getLookUpTable)) {
		
		  my $Level = $HeavyAtoms[$j]->getLevel;
		
		  my $CurrentChemicalSpaceMatch = ChemicalSpaceMatch->New($HeavyAtoms[$i], $HeavyAtoms[$j], $Level);
		
		  foreach my $ChemicalSpaceMatch (@ChemicalSpaceMatches) {
			
			  $CurrentChemicalSpaceMatch->CheckForPrerequisite($ChemicalSpaceMatch);
			
			  $ChemicalSpaceMatch->CheckForPrerequisite($CurrentChemicalSpaceMatch);
			  
		  }
		
		  push @ChemicalSpaceMatches, $CurrentChemicalSpaceMatch;
		
		  $HashedAssociations{$CurrentChemicalSpaceMatch->Stringify} = $CurrentChemicalSpaceMatch;
		
		}
		
	}
	
}


my @Files = @ARGV;

shift @Files;

my %BestAssociation;
my $CoordinateHash;


foreach my $File (@Files) {

  my $PDBFileHandler2 = MoleculeFileHandler->New($File);

	$CoordinateHash = $PDBFileHandler2->ExtractCoordinatesHash;
	
  %BestAssociation = ();	

	foreach my $i (0 .. @HeavyAtoms - 1) {
	
	  $BestAssociation{$HeavyAtoms[$i]->getName} = [ [CalculateSquaredDistance($HeavyAtoms[$i]->getCartesianCoordinates, getNewCoords($HeavyAtoms[$i])), $HeavyAtoms[$i]] ];

  }	


my %UsedChemicalSpaceMatches;

foreach my $ChemicalSpaceMatch (@ChemicalSpaceMatches) {
	
	next if $ChemicalSpaceMatch->getIsAPrerequisite == 1;
	
  my $Prerequisites = $ChemicalSpaceMatch->getPrerequisites;
  
  my ($CurrentSquaredDistanceTotal, $NewSquaredDistanceTotal) = (0, 0);
  
  foreach my $PrerequisiteChemicalSpaceMatch ($ChemicalSpaceMatch, @$Prerequisites) {
		
	  $PrerequisiteChemicalSpaceMatch->AddAPrerequisiteOf($ChemicalSpaceMatch);
		
	  my $FirstDistance  = getSqrdDist($PrerequisiteChemicalSpaceMatch->getPAtom, getNewCoords($PrerequisiteChemicalSpaceMatch->getSAtom)); 
	  my $SecondDistance = getSqrdDist($PrerequisiteChemicalSpaceMatch->getSAtom, getNewCoords($PrerequisiteChemicalSpaceMatch->getPAtom));
	  
	  $NewSquaredDistanceTotal += $FirstDistance + $SecondDistance;

	  $CurrentSquaredDistanceTotal += getBestDistance($PrerequisiteChemicalSpaceMatch->getPAtom) + getBestDistance($PrerequisiteChemicalSpaceMatch->getSAtom);
		
  }

  next if $NewSquaredDistanceTotal > $CurrentSquaredDistanceTotal;
	
	my @AtomsToCheck;
	
	foreach my $PrerequisiteChemicalSpaceMatch ($ChemicalSpaceMatch, @$Prerequisites) {
	
    if(getBestAssociation($PrerequisiteChemicalSpaceMatch->getPAtom) ne $PrerequisiteChemicalSpaceMatch->getPAtom && getBestAssociation($PrerequisiteChemicalSpaceMatch->getPAtom) ne $PrerequisiteChemicalSpaceMatch->getSAtom) {
			
			if(exists $HashedAssociations{$PrerequisiteChemicalSpaceMatch->getPAtom->getName . " " . getBestAssociation($PrerequisiteChemicalSpaceMatch->getPAtom)->getName}) {
				
			  push @AtomsToCheck, $PrerequisiteChemicalSpaceMatch->getPAtom;
			
			}
			
			else {
				
			  push @AtomsToCheck, getBestAssociation($PrerequisiteChemicalSpaceMatch->getPAtom);
			  
			}
			
	  }
	
	  if(getBestAssociation($PrerequisiteChemicalSpaceMatch->getSAtom) ne $PrerequisiteChemicalSpaceMatch->getSAtom && getBestAssociation($PrerequisiteChemicalSpaceMatch->getSAtom) ne $PrerequisiteChemicalSpaceMatch->getPAtom) {

				if(exists $HashedAssociations{$PrerequisiteChemicalSpaceMatch->getSAtom->getName . " " . getBestAssociation($PrerequisiteChemicalSpaceMatch->getSAtom)->getName}) {

				  push @AtomsToCheck, $PrerequisiteChemicalSpaceMatch->getSAtom;

				}

				else {

				  push @AtomsToCheck, getBestAssociation($PrerequisiteChemicalSpaceMatch->getSAtom);

				}
		}
		
	}
	
	my %Seen; 
			
  my @ChemicalSpaceMatchesToCheck = grep { ! $Seen{$_} ++ && $UsedChemicalSpaceMatches{$_}  }
	 																  map  { @{$HashedAssociations{$_->getName . " " . getBestAssociation($_)->getName }->getAPrerequisitesOf} } @AtomsToCheck;
	
	#print "REMOVED: ";
		
	foreach my $ChemicalSpaceMatchToCheck (@ChemicalSpaceMatchesToCheck) {

    #print $ChemicalSpaceMatchToCheck->getPAtom->getName . " => " . $ChemicalSpaceMatchToCheck->getSAtom->getName . " ";

		shiftBestDistance($ChemicalSpaceMatchToCheck->getPAtom); shiftBestDistance($ChemicalSpaceMatchToCheck->getSAtom);

  }

  #print "\nACCEPTED: ";

	foreach my $PrerequisiteChemicalSpaceMatch ($ChemicalSpaceMatch, @$Prerequisites) {

		#print $PrerequisiteChemicalSpaceMatch->getPAtom->getName . " => " . $PrerequisiteChemicalSpaceMatch->getSAtom->getName . " ";

	  my $FirstDistance  = getSqrdDist($PrerequisiteChemicalSpaceMatch->getPAtom, getNewCoords($PrerequisiteChemicalSpaceMatch->getSAtom)); 
		my $SecondDistance = getSqrdDist($PrerequisiteChemicalSpaceMatch->getSAtom, getNewCoords($PrerequisiteChemicalSpaceMatch->getPAtom));
		  
	  setBestDistance($PrerequisiteChemicalSpaceMatch->getPAtom, $FirstDistance, $PrerequisiteChemicalSpaceMatch->getSAtom); 
		setBestDistance($PrerequisiteChemicalSpaceMatch->getSAtom, $SecondDistance, $PrerequisiteChemicalSpaceMatch->getPAtom);

	}

	#print "\n";
	  
	
	$UsedChemicalSpaceMatches{$ChemicalSpaceMatch} = 1; 
		
}

my $BestRMSD = 0;
my $RMSD = 0;

foreach my $Atom (@HeavyAtoms) {
	
	#print $Atom->getName . " => " . getBestAssociation($Atom)->getName . "\n";
	
	$BestRMSD += getBestDistance($Atom);
	
	$RMSD += CalculateSquaredDistance($Atom->getCartesianCoordinates, getNewCoords($Atom));
	
}

$BestRMSD /= @HeavyAtoms;
$RMSD /= @HeavyAtoms;

print sqrt($RMSD) . " " . sqrt($BestRMSD) . "\n";

}

sub getBestAssociation {
	
	my $Atom = shift;
	
	return $BestAssociation{$Atom->getName}->[0]->[1];
	
}

sub shiftBestDistance {
	
	my $Atom = shift;
	
	shift @{$BestAssociation{$Atom->getName}};
	
}

sub setBestDistance {
	
	my ($RefAtom, $Distance, $Atom) = @_;
	
	unshift @{$BestAssociation{$RefAtom->getName}}, [$Distance, $Atom];
	
}

sub getBestDistance {
	
	my $Atom = shift;
	
	return $BestAssociation{$Atom->getName}->[0]->[0];
	
}

sub getSqrdDist {
	
	my ($Atom, $Vect) = @_;
	
	return CalculateSquaredDistance($Atom->getCartesianCoordinates, $Vect); 
	
}

sub CalculateSquaredDistance {
	
	my ($Vector1, $Vector2) = @_;
	
	return ($Vector1->[0] - $Vector2->[0]) ** 2 + ($Vector1->[1] - $Vector2->[1]) ** 2 + ($Vector1->[2] - $Vector2->[2]) ** 2; 
	
}

sub getNewCoords {
	
	my $Atom = shift;
	
	return $CoordinateHash->{$Atom->getChain->getName . "-" . $Atom->getMolecule->getName . "-" . $Atom->getChain->getSegmentId . "-" . $Atom->getName};
	
}	
