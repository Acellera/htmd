#!/usr/bin/perl -w

use strict;
use warnings;
use Cwd qw(abs_path);
use Benchmark;

use lib '../lib';
use Test::More tests => 4;

BEGIN { use_ok('LookUpTable') };

use LookUpTable;
use MoleculeFileHandler ':all';
use Parameters ':all';
use BaseObject ':vars';

srand(time ^ ($$ + ($$ << 15)));

#Deal With Initial Parameter Setup, This will be handled in main .pl script
my $DefaultParameters = Parameters->New;

my ($Path) = (abs_path($0) =~ /(\S+)\/\S+$/);

my ($DefaultParametersPath) = ($Path =~ /(\S+)\/\S+$/);

$DefaultParametersPath .= "/resources/DefaultParameters.par";

$DefaultParameters->Initiate($Path, $DefaultParametersPath);

#Set Global Parameters
$Parameters = $DefaultParameters;

my $MoleculeFileHandler = MoleculeFileHandler->New("MoleculeFileHandler/TestPDBFiles/1byj_ligand_correct.pdb");

my $Structures = $MoleculeFileHandler->BuildObjectsFromFile;

my ($Molecule) = @{$Structures->[0]->getMolecules};

my $Atoms = $Molecule->getAtoms;

foreach my $Atom (@$Atoms) {
	
	my $AtomLookUpTable = LookUpTable->New;
		
	$AtomLookUpTable->Initiate($Atom);
  
  $Atom->setLookUpTable($AtomLookUpTable);	 
	
}

my $CopiedMolecule = $Molecule->Copy;

my $CopiedAtoms = $CopiedMolecule->getAtoms;

foreach my $CopiedAtom (@$CopiedAtoms) {
	
	my $AtomLookUpTable = LookUpTable->New;
		
	$AtomLookUpTable->Initiate($CopiedAtom);
  
  $CopiedAtom->setLookUpTable($AtomLookUpTable);	 
	
}

subtest "Self Comparison of Copied LookUpTables" => sub {
		
	plan tests => scalar(@$Atoms);
	
	foreach my $Atom (@$Atoms) {
		
		foreach my $CopiedAtom (@$CopiedAtoms) {
			
			next if $Atom->getName ne $CopiedAtom->getName;
			
			my $ExpectedStringifyLength = -2; 

			$ExpectedStringifyLength += length($_->getAtom->Stringify) + 2 foreach (@{$Atom->getLookUpTable->getNodes});
			
			my @CompareTables = ( $CopiedAtom->getLookUpTable );
			
			my $Consensus = $Atom->getLookUpTable->DetermineConsensusLookUpTable(\@CompareTables);
			
			print $Atom->getLookUpTable->AreAllNodesSharedBetween(\@CompareTables);
												
			#is(length($Consensus->Stringify), $ExpectedStringifyLength, 'Stringify Test Pans out with ' . $Atom->getName . ' as the Head');
					
		}
		
	}
	
};

exit 1;

my $CopiedLookUpTable = $Atoms->[0]->getLookUpTable->Copy;

my $ExpectedStringifyLength = -2; 

$ExpectedStringifyLength += length($_->getAtom->Stringify) + 2 foreach (@{$Atoms->[0]->getLookUpTable->getNodes});
	
my $Consensus = $Atoms->[0]->getLookUpTable->DetermineConsensusLookUpTable($CopiedLookUpTable);
				
is(length($Consensus->Stringify), $ExpectedStringifyLength, 'Stringify Test Pans out with ' . $Atoms->[0]->getName . ' as the Head');





