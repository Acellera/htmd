#!/usr/bin/perl

use strict;
use warnings;

use lib $ENV{'PerlChemistry'} . "/lib";

use Parameters ':all';
use BaseObject ':vars';
use MoleculeFileHandler ':all';

#Deal With Initial Parameter Setup, This will be handled in main .pl script
my $DefaultParameters = Parameters->New;

$DefaultParameters->Initiate;

#Set Global Parameters
$Parameters = $DefaultParameters;

my $PDBFileHandler = MoleculeFileHandler->New($ARGV[0]);

my $PDBStructures = $PDBFileHandler->BuildObjectsFromFile;

my $Molecule = $PDBStructures->[0]->getMolecule(0);

$Molecule->Initiate();

my $Atoms = $Molecule->getAtoms;

my $PDBFileHandler2 = MoleculeFileHandler->New($ARGV[1]);

my $CoordinateHash = $PDBFileHandler2->ExtractCoordinatesHash;

my @Associations;
my %AssociationsDependencies;
my %HashedAssociations;

my %AcceptedAssociationsSwaps;
my %ContainedInADependecy;

my %BestAssociation;
my %AssociationLevel;

my @HeavyAtoms =  @$Atoms;

foreach my $i (0 .. @HeavyAtoms - 1) {
		
	$HeavyAtoms[$i]->getLookUpTable->SetAtomLevelsUsingNodeLevels;
	
	$BestAssociation{$HeavyAtoms[$i]->getName} = [ [CalculateSquaredDistance($HeavyAtoms[$i]->getCartesianCoordinates, getNewCoords($HeavyAtoms[$i])), $HeavyAtoms[$i]] ];
		
	foreach my $j ($i + 1 .. @HeavyAtoms - 1) {
		
		my $Level = $HeavyAtoms[$j]->getLevel;
	  	
		if($HeavyAtoms[$i]->DoesLookUpTableMatch($HeavyAtoms[$j]->getLookUpTable)) {
			
			$AssociationsDependencies{$HeavyAtoms[$i]->getName . "-" . $HeavyAtoms[$j]->getName} = [];
			
			$AssociationLevel{$HeavyAtoms[$i]->getName . "-" . $HeavyAtoms[$j]->getName} = $Level;
			
			foreach my $Association (@Associations) {
			
			  if((($HeavyAtoms[$i]->IsBondedTo($Association->[0]) && $HeavyAtoms[$j]->IsBondedTo($Association->[1])) || ($HeavyAtoms[$i]->IsBondedTo($Association->[1]) && $HeavyAtoms[$j]->IsBondedTo($Association->[0]))) && 	$AssociationLevel{$Association->[0]->getName . "-" . $Association->[1]->getName} < $Level) {
				
				  $AssociationsDependencies{$HeavyAtoms[$i]->getName . "-" . $HeavyAtoms[$j]->getName} = [$Association, @{$AssociationsDependencies{$Association->[0]->getName . "-" . $Association->[1]->getName}}] ;
				
				  $ContainedInADependecy{$Association->[0]->getName . "-" . $Association->[1]->getName} = 1;
				
			  }	
				
			}
			
			push @Associations, [$HeavyAtoms[$i], $HeavyAtoms[$j]];		
						
		}
					
	}
	
}

my @HeadAssociations = sort { @{$AssociationsDependencies{$a->[0]->getName . "-" . $a->[1]->getName}} <=> @{$AssociationsDependencies{$b->[0]->getName . "-" . $b->[1]->getName}} } grep { ! exists $ContainedInADependecy{$_->[0]->getName . "-" . $_->[1]->getName} } @Associations; 

foreach my $Association (@HeadAssociations) {
	
	$HashedAssociations{$Association->[0]->getName . "-" . $Association->[1]->getName } = $Association;
	
	my @CurrentAssociations = ($Association, @{$AssociationsDependencies{$Association->[0]->getName . "-" . $Association->[1]->getName}});
	
	my $OldTotal = 0;
	my $NewTotal = 0;
	
  foreach my $CurrentAssociation (@CurrentAssociations) {
		
		$ContainedInADependecy{$CurrentAssociation->[0]->getName . "-" . $CurrentAssociation->[1]->getName} = $Association->[0]->getName . "-" . $Association->[1]->getName;
		
	  my ($FirstDistance, $SecondDistance) = (CalculateSquaredDistance($CurrentAssociation->[0]->getCartesianCoordinates, getNewCoords($CurrentAssociation->[1])), CalculateSquaredDistance($CurrentAssociation->[1]->getCartesianCoordinates, getNewCoords($CurrentAssociation->[0])));
	  
	  $NewTotal += $FirstDistance + $SecondDistance;
	
	  $OldTotal += getBestDistance($CurrentAssociation->[0]) + getBestDistance($CurrentAssociation->[1]);
	
  }

  if($NewTotal < $OldTotal) {
	
	  my @AtomsToCheck;
	
	  foreach my $CurrentAssociation (@CurrentAssociations) {
		
		  if(getBestAssociation($CurrentAssociation->[0]) ne $CurrentAssociation->[0] && getBestAssociation($CurrentAssociation->[0]) ne $CurrentAssociation->[1] ) {
			
			  push @AtomsToCheck, $CurrentAssociation->[0];
			
		  }
		
		  if(getBestAssociation($CurrentAssociation->[1]) ne $CurrentAssociation->[1] && getBestAssociation($CurrentAssociation->[1]) ne $CurrentAssociation->[0] ) {
			
			  push @AtomsToCheck, getBestAssociation($CurrentAssociation->[1]);
			
		  }
		
	  }
	
	  my %Seen; 
	
	  my @AssociationsToCheck = map  { $HashedAssociations{ $ContainedInADependecy{ $_->getName . "-" . getBestAssociation($_)->getName } } } 
															grep { ! $Seen{ $ContainedInADependecy{ $_->getName . "-" . getBestAssociation($_)->getName } }++  } @AtomsToCheck;  
			
	  print "REMOVED: ";
		
		foreach my $AssociationToCheck (@AssociationsToCheck) {
			
			print $AssociationToCheck->[0]->getName . " => " . $AssociationToCheck->[1]->getName . " ";
			
			shiftBestDistance($AssociationToCheck->[0]); shiftBestDistance($AssociationToCheck->[1]);
			
		}
		
		print "\nACCEPTED: ";
		
		foreach my $CurrentAssociation (@CurrentAssociations) {
			
			print $CurrentAssociation->[0]->getName . " => " . $CurrentAssociation->[1]->getName . " ";
			
		  my ($FirstDistance, $SecondDistance) = (CalculateSquaredDistance($CurrentAssociation->[0]->getCartesianCoordinates, getNewCoords($CurrentAssociation->[1])), CalculateSquaredDistance($CurrentAssociation->[1]->getCartesianCoordinates, getNewCoords($CurrentAssociation->[0])));
			
			setBestDistance($CurrentAssociation->[0], $FirstDistance, $CurrentAssociation->[1]); 
			setBestDistance($CurrentAssociation->[1], $SecondDistance, $CurrentAssociation->[0]);
			
		}
		
		print "\n";
			
  }
	
	
}


foreach my $Atom (@HeavyAtoms) {
	
	print $Atom->getName . " => " . getBestAssociation($Atom)->getName . "\n";
	
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

sub CalculateSquaredDistance {
	
	my ($Vector1, $Vector2) = @_;
	
	return ($Vector1->[0] - $Vector2->[0]) ** 2 + ($Vector1->[1] - $Vector2->[1]) ** 2 + ($Vector1->[2] - $Vector2->[2]) ** 2; 
	
}

sub getNewCoords {
	
	my $Atom = shift;
	
	return $CoordinateHash->{$Atom->getChain->getName . "-" . $Atom->getMolecule->getName . "-" . $Atom->getChain->getSegmentId . "-" . $Atom->getName};
	
}


