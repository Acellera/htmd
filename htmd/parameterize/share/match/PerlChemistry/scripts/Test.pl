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

my %AssociationLevels;
my %HasMultipleAssociations;

my @HeavyAtoms =  @$Atoms;

foreach my $Atom (@HeavyAtoms) {
	
	$HasMultipleAssociations{$Atom} = 0;
	
}

foreach my $i (0 .. @HeavyAtoms - 1) {
		
	$HeavyAtoms[$i]->getLookUpTable->SetAtomLevelsUsingNodeLevels;
		
	foreach my $j ($i + 1 .. @HeavyAtoms - 1) {
		
		my $Level = $HeavyAtoms[$j]->getLevel;
	  	
		if ($HeavyAtoms[$i]->DoesLookUpTableMatch($HeavyAtoms[$j]->getLookUpTable)) {
			
			push @{$AssociationLevels{$Level}}, [$HeavyAtoms[$i], $HeavyAtoms[$j]];
			
			$HasMultipleAssociations{$HeavyAtoms[$i]} = $Level;
			$HasMultipleAssociations{$HeavyAtoms[$j]} = $Level;
			
			
		}
			
	}
	
}

my $Count = 0;
my $Done = 0;

#"$ChainName-$ResidueName-$SegmentId-" . $_->{'Name'

my %SavedDistances;

my $RMSD = 0;

my %BestAssociation;

foreach my $Atom (@HeavyAtoms) {
	
	my $Coords1 = $Atom->getCartesianCoordinates;
	
	my $Coords2 = $CoordinateHash->{$Atom->getChain->getName . "-" . $Atom->getMolecule->getName . "-" . $Atom->getChain->getSegmentId . "-" . $Atom->getName};

  my $SquaredDistance =  CalculateSquaredDistance($Coords1, $Coords2);

  $BestAssociation{$Atom->getName} = [[$SquaredDistance, $Atom]];
	
	$SavedDistances{$Atom->getName}{$Atom->getName} = $SquaredDistance;
	
	$RMSD += $SquaredDistance;
	
}

$RMSD *= 1/@HeavyAtoms;

print sqrt($RMSD) . "\n";

my $BestRMSD = 0;
my $DoesNotExists = 0;

while (!$Done) {

  $Count++;

  $Done = 1 if $DoesNotExists > 5;

  unless(exists $AssociationLevels{$Count}) {
	
	  $DoesNotExists++; next
	
  }

  $DoesNotExists = 0;

  my $CurrentLevelAssociations = $AssociationLevels{$Count};
  
  foreach my $Association (@$CurrentLevelAssociations) {
		
	  my ($FirstDistance, $SecondDistance) = (CalculateSquaredDistance($Association->[0]->getCartesianCoordinates, getNewCoords($Association->[1])), CalculateSquaredDistance($Association->[1]->getCartesianCoordinates, getNewCoords($Association->[0])));
	
	  $SavedDistances{$Association->[0]->getName}{$Association->[1]->getName} = $FirstDistance;
	  $SavedDistances{$Association->[1]->getName}{$Association->[0]->getName} = $SecondDistance;
	  
	  if(getBestDistance($Association->[0]) + getBestDistance($Association->[1]) >  $FirstDistance + $SecondDistance) {
		
		  setBestDistance($Association->[0], $FirstDistance, $Association->[1]);
		  setBestDistance($Association->[1], $SecondDistance, $Association->[0]);
		
	  }
	
	  else {
		
      #if(getBestAssociation($Association->[0]) eq $Association->[0]) {
	
	      print $Association->[0]->getName . " " . $Association->[1]->getName . "\n";

        my @PossibleAssociationReversalAtoms = grep { $HasMultipleAssociations{$_} < $  } @{$Association->[0]->getBondedAtoms};


				#print getBestDistance($Association->[0]) . " " .  getBestDistance($Association->[1]) . " " . $FirstDistance . " " . $SecondDistance . "\n";
	
	      #print "made it!\n";
	 
	      #my @CurrentAtoms = 
	
	      #while () {
		
	      #}
	
	      #exit 1 if $Count > 2;
	
     # }
		
	  }
	
	
  }

}

foreach my $Atom (@HeavyAtoms) {
	
	$BestRMSD += getBestDistance($Atom);
	
}

$BestRMSD *= 1/@HeavyAtoms;

print sqrt($BestRMSD) . "\n";


sub getBestAssociation {
	
	my $Atom = shift;
	
	return $BestAssociation{$Atom->getName}->[0]->[1];
	
}

sub setBestDistance {
	
	my ($RefAtom, $Distance, $Atom) = @_;
	
	unshift @{$BestAssociation{$Atom->getName}}, [$Distance, $Atom];
	
}

sub getBestDistance {
	
	my $Atom = shift;
	
	return $BestAssociation{$Atom->getName}->[0]->[0];
	
}

sub getNewCoords {
	
	my $Atom = shift;
	
	return $CoordinateHash->{$Atom->getChain->getName . "-" . $Atom->getMolecule->getName . "-" . $Atom->getChain->getSegmentId . "-" . $Atom->getName};
	
}

sub CalculateSquaredDistance {
	
	my ($Vector1, $Vector2) = @_;
	
	return ($Vector1->[0] - $Vector2->[0]) ** 2 + ($Vector1->[1] - $Vector2->[1]) ** 2 + ($Vector1->[2] - $Vector2->[2]) ** 2; 
	
}
