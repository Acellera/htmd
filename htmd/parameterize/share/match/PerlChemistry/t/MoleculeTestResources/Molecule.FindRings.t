#!/usr/bin/perl -w

use strict;
use warnings;
use Carp;

use lib $ENV{'PerlChemistry'} . '/lib';
use Test::More tests => 4;

BEGIN { use_ok('Molecule') };

use Molecule ':all';
use MoleculeFileHandler ':all';
use Parameters ':all';
use BaseObject ':vars';

#use TestingFunctions ':func';

#Deal With Initial Parameter Setup, This will be handled in main .pl script
my $DefaultParameters = Parameters->New;

$DefaultParameters->Initiate;

#Set Global Parameters
$Parameters = $DefaultParameters;

my $FileName = "RingTest1.pdb";

my $Name = substr($FileName, 0, -4);

my $MoleculeFileHandler = MoleculeFileHandler->New($ENV{'PerlChemistry'} . "t/MoleculeTestResources/FindRingExamples/$FileName");

my $Structures = $MoleculeFileHandler->BuildObjectsFromFile;

my ($Molecule) = @{$Structures->[0]->getMolecules};

$Molecule->Initiate;

my $MoleculeRings = $Molecule->getRings;

open(FILE, $ENV{'PerlChemistry'} . "t/MoleculeTestResources/FindRingExamples/$Name.dat");

my @FileContents = <FILE>;

close(FILE);

foreach my $Ring (@$MoleculeRings) {
	
	#print join(" ", map { $_->getName} @$Ring) . "\n";
	
}

foreach my $Line (@FileContents) {
	
	my @AtomNames = split /\s+/, $Line;
	
	my $Found = 0;
	
	foreach my $Ring (@$MoleculeRings) {
	
	  next if @$Ring != @AtomNames;
	  
    my $Matched = 0; 
	
	  foreach my $RingAtom (@$Ring) {
		
		  foreach my $AtomName (@AtomNames) {
			
			  $Matched++ if $RingAtom->getName eq $AtomName;
			
			}
		
	  }
	
	  if($Matched eq scalar(@AtomNames)) {
		
		  $Found = 1; last;
		
	  }	
		
	}
	
	croak "Cannot find Ring: " . join(" ", @AtomNames) . " in $FileName\n" if $Found == 0;
	
}