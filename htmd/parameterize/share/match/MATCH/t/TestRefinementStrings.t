#!/usr/bin/env perl

use 5.010000;
use strict;
use warnings;
use Carp;
use Cwd qw(abs_path);

#Load MATCH Libraries 
use lib $ENV{'MATCH'} . "/lib"; 

use MATCHBaseObject;
use BaseObject ':vars';

use MATCHParameters;
use MATCHFunctions ':func';
use AtomCharger;
use AtomTyper;
use AtomTypeSubstituter;
use Storable;

our $VERSION = '0.01';

#Exported Variables


#Non Exported Variables


# Preloaded methods go here.

srand(time ^ ($$ + ($$ << 15)));

#Deal With Initial Parameter Setup, This will be handled in main .pl script
my $DefaultParameters = MATCHParameters->New;

$DefaultParameters->Initiate("../resources/" . $ARGV[0] . ".par");

$Parameters = $DefaultParameters;

my $AtomTypeSubstituter = AtomTypeSubstituter->New;

$AtomTypeSubstituter->Initiate;

my $AtomTyper = AtomTyper->New;

$AtomTyper->Initiate($AtomTypeSubstituter);

my $AtomCharger = AtomCharger->New;

$AtomCharger->Initiate($AtomTypeSubstituter);

my $Chains =  retrieve("StoredTopologys/" . $ARGV[0] . ".dat");

#my $Chains = retrieve("../scripts/Chains.dat");

my @Atoms = map { @{$_->getAtoms} } map { @{$_->getMolecules} } @$Chains;

my $PrimaryType = $ARGV[1];
my $SecondaryType = $ARGV[2];

my @TypeAtoms = grep { $_->getType eq $PrimaryType } @Atoms;

my @IncrementMolecules;

foreach my $Atom (@TypeAtoms) {
	
	my $Bonded = $Atom->getBondedAtoms;
	
	foreach my $BondedAtom (@$Bonded) {
		
		push @IncrementMolecules, $Atom->getMolecule if $BondedAtom->getType eq $SecondaryType;
		
	}
	
}

my %SeenMolecule;

my @UniqueIncrementMolecules = grep { ! $SeenMolecule{$_->getName } ++  } @IncrementMolecules;

my $Sum = 0;



foreach my $Molecule (@UniqueIncrementMolecules) {
	
	#next if $AtomTyper->TypeAtomsInMolecule($Molecule) == 0;
	
  my $CopiedMolecule = $Molecule->Copy;

  foreach my $Atom (@{$CopiedMolecule->getAtoms}) {
	
	  $Atom->setCharge(0);
	
  } 


  next if $AtomCharger->ChargeAtomsInMolecule($CopiedMolecule) != 1;

  my $TotalForMolecule = 0;

  print $Molecule->getName ."\n";
 
  my %AtomChargeHash = map { $_->getName => $_->getCharge } @{$Molecule->getAtoms};
  my %CopiedAtomChargeHash = map { $_->getName => $_->getCharge } @{$CopiedMolecule->getAtoms};
  
  my @SavedPrintOuts;

  foreach my $Atom (@{$CopiedMolecule->getAtoms}) {
	
	  if($Atom->getType eq $PrimaryType || $Atom->getType eq $SecondaryType) {
		
		  my $Value = $AtomChargeHash{$Atom->getName} - $CopiedAtomChargeHash{$Atom->getName} ;
		
		  if ($Value > 0.0000001 || $Value < -0.00000001  ) {
			
			  push(@SavedPrintOuts, [$Atom, $Atom->getName . " " . $Atom->getType . " " . $Value]);
			
		  }
		
	  }
	
  }

  my %RealPrintOuts;

  foreach my $i (0 .. @SavedPrintOuts-1) {
	
	  foreach my $j ($i+1 .. @SavedPrintOuts-1) {
		
		  if($SavedPrintOuts[$i]->[0]->getType ne $SavedPrintOuts[$j]->[0]->getType && $SavedPrintOuts[$i]->[0]->IsBondedTo($SavedPrintOuts[$j]->[0])) {
			
			  $RealPrintOuts{$SavedPrintOuts[$i]->[0]} = $SavedPrintOuts[$i];
			  $RealPrintOuts{$SavedPrintOuts[$j]->[0]} = $SavedPrintOuts[$j];
			
			}
			
		}
			
  }
  
#	next if keys %RealPrintOuts < 2;
	
	print $Molecule->getName . "\n";
	
	foreach my $Atom (@{$CopiedMolecule->getAtoms}) {

	#  if($Atom->getType eq $PrimaryType || $Atom->getType eq $SecondaryType) {

		  my $Value = $AtomChargeHash{$Atom->getName} - $CopiedAtomChargeHash{$Atom->getName} ;

		  if ($Value > 0.0000001 || $Value < -0.00000001  ) {

			  print $Atom->getName . " " . $Atom->getType . " " . $Value . "\n";
		  
		 }

	 # }

  }
 
	
	foreach my $Value (values %RealPrintOuts) {
		
	#	print $Value->[1] . "\n";
		
		my @split = split(/\s+/, $Value->[1]);
				
		$Sum += abs($split[-1]);
		
	}	
		
}

print "Total Charge Missing: " . $Sum . "\n";


