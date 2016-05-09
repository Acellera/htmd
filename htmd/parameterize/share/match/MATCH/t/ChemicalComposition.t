#!/usr/bin/env perl

use strict;
use warnings;

use lib $ENV{'MATCH'} . "/lib";  
use Carp;
use Storable;

use MATCHBaseObject;
use BaseObject ':vars';

use MATCHParameters;
use MATCHFunctions ':all';
use MoleculeFileHandler;

use AtomTyper;
use AtomTypeSubstituter;
use MoleculeParameterizer;

my $DefaultParameters = MATCHParameters->New;

$DefaultParameters->Initiate("../resources/" . $ARGV[0] . ".par");
#$DefaultParameters->Initiate;

$Parameters = $DefaultParameters;

my $Chains =  retrieve("StoredTopologys/" . $ARGV[0] . ".dat");

my $Count = 0;
my $Total = 0;

foreach my $Chain (@$Chains) {
	
	my $Atoms = $Chain->getMolecule(0)->getAtoms;
	
	foreach my $Atom (@$Atoms) {
		
    if($Atom->getElement eq "P" || $Atom->getElement eq "F" || $Atom->getElement eq "CL" || $Atom->getElement eq "BR" || $Atom->getElement eq "I") {
	
	    $Count++; last;
	
    }	
		
	}
	
	$Total++;
	
	
}

print $Count . " "  . $Total . "\n";