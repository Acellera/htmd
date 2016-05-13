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

my $AtomCharger = AtomCharger->New;

$AtomCharger->Initiate($AtomTypeSubstituter);

#my $Chains = SetupMoleculesFromTopology("../resources/top_all22_prot/top_all22_prot.inp",
# 																		    "../resources/top_all22_prot/top_all22_prot.patches");

my $Chains = retrieve("StoredTopologys/" . $ARGV[0] . ".dat");

my @Atoms = map { @{$_->getAtoms} } map { @{$_->getMolecules} } @$Chains;

my @Molecules = map { $_->getMolecule(0) } @$Chains;

my $Count = 1;

foreach my $Molecule (@Molecules) {
	
	my $CopiedMolecule = $Molecule->Copy;
	
	foreach my $Atom (@{$CopiedMolecule->getAtoms}) {
	
	  $Atom->setCharge(0);
	
  } 
  
  next if $CopiedMolecule->getName ne $ARGV[1];

  print $CopiedMolecule->getName . "\n";
  
  next if $AtomCharger->ChargeAtomsInMolecule($CopiedMolecule) != 1;

  my %AtomChargeHash = map { $_->getName => $_->getCharge } @{$Molecule->getAtoms};
  my %CopiedAtomChargeHash = map { $_->getName => $_->getCharge } @{$CopiedMolecule->getAtoms};

  my @SavedPrintOuts;

  my $Total = 0;

  #print $CopiedMolecule->getName . "\n";

  my $Next = 0;


  foreach my $Atom (@{$CopiedMolecule->getAtoms}) {
		
	  #$Next = 1 if @{$Atom->getRings} > 1;
			
	  my $Value = $AtomChargeHash{$Atom->getName} - $CopiedAtomChargeHash{$Atom->getName} ;
	
		
  	print $Atom->getName . " " . $Atom->getType . " " . $Value .  "\n" if (abs($Value) > 0.01);
		
		$Total += abs($Value);
		
  }
#  next if $Next;

  next if $Total < 0.10;

  print $CopiedMolecule->getName . "  " . $Total . " " . "\n";
  

  $Count++;
		
}


