#!/usr/bin/perl -w

use strict;
use warnings;
use Cwd qw(abs_path);

use lib '../lib';
use Test::More tests => 2;

BEGIN { use_ok('LookUpTable') };

use LookUpTable;
use MoleculeFileHandler ':all';
use Parameters ':all';
use BaseObject ':vars';

#Deal With Initial Parameter Setup, This will be handled in main .pl script
my $DefaultParameters = Parameters->New;

my ($Path) = (abs_path($0) =~ /(\S+)\/\S+$/);

my ($DefaultParametersPath) = ($Path =~ /(\S+)\/\S+$/);

$DefaultParametersPath .= "/resources/DefaultParameters.par";

$DefaultParameters->Initiate($Path, $DefaultParametersPath);

#Set Global Parameters
$Parameters = $DefaultParameters;

#Run this test if something changes in LookUpTable::Stringify
subtest "Test Stringify" => sub { 
	
	chdir "MoleculeFileHandler/TestMol2Files";

	my @Mol2Molecules = <*>;

	chdir "../..";
	
	plan tests => scalar(@Mol2Molecules);

	foreach my $Mol2File (@Mol2Molecules) {

	  my $MoleculeFileHandler = MoleculeFileHandler->New("MoleculeFileHandler/TestMol2Files/" . $Mol2File);

	  my $Structures = $MoleculeFileHandler->BuildObjectsFromFile;
	
	  my ($Molecule) = @{$Structures->[0]->getMolecules};
	
	  my $Atoms = $Molecule->getAtoms;
	
	  my $LookUpTable = LookUpTable->New;

	  $LookUpTable->Initiate($Atoms->[0]);
	
	  my $ExpectedStringifyLength = -2; 

	  $ExpectedStringifyLength += length($_->Stringify) + 2 foreach (@$Atoms);

	  is($ExpectedStringifyLength, length($LookUpTable->Stringify), 'Stringify is expected length for ' . $Structures->[0]->getFileName); 

	}
	

};

