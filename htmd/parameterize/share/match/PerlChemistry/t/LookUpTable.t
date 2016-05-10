#!/usr/bin/perl -w

use strict;
use warnings;
use Cwd qw(abs_path);

use lib '../lib';
use Test::More tests => 4;

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

my $MoleculeFileHandler = MoleculeFileHandler->New("MoleculeFileHandler/TestPDBFiles/1byj_ligand_correct.pdb");

my $Structures = $MoleculeFileHandler->BuildObjectsFromFile;

my ($Molecule) = @{$Structures->[0]->getMolecules};

my $Atoms = $Molecule->getAtoms;

subtest "Add Atom To LookUpTable" => sub {
	
	plan tests => 1;
	
	my $CopyMolecule = $Molecule->Copy;
	
	my $CopiedAtoms = $CopyMolecule->getAtoms;
	
	my $LookUpTable = LookUpTable->New;
	
	$LookUpTable->Initiate($CopiedAtoms->[0]);
	
	my $NewAtom = Atom->New({Name => "CL1", Element => "CL", Level => 1 });
	
	$NewAtom->Initiate;
	
	my $TestBond = Bond->New($CopiedAtoms->[0], $NewAtom, 1);
	
	$LookUpTable->AddAtom($NewAtom);
	
	my $NewNode;
	
	foreach (@{$LookUpTable->getNodes}) {
	
	  if($_->getAtomName eq $NewAtom->getName) { $NewNode = $_; last; }
		
	}
	
	$TestBond->UnBondAtoms;
	
	$LookUpTable->RemoveNode($NewNode);
	
	print $LookUpTable->Stringify . "\n";
	
	
};

exit 1;

subtest "BuildNodesFromAtoms" => sub {
	
	plan tests => 1;
	
	my $LookUpTable = LookUpTable->New;
	
	$LookUpTable::HARD_DEPTH_LIMIT = 100;

	$LookUpTable->BuildNodesFromAtoms($Atoms->[0]);
	
	my $Nodes = $LookUpTable->getNodes;
	
	is(@$Nodes, @$Atoms, 'Same number of Nodes and Atoms');
	
};

subtest "Test GridElements" => sub {
	
	plan tests => 1;
		
	my $LookUpTable = LookUpTable->New;
	
	$LookUpTable::HARD_DEPTH_LIMIT = 100;

	$LookUpTable->Initiate($Atoms->[0]);
		
	my $Nodes = $LookUpTable->getNodes;
	
	my @SortedNodes = sort { $b->getLevel <=> $a->getLevel } @$Nodes;
			
	
	is(@SortedNodes, @$Atoms, 'Count is correct');
	
	
	
};

subtest "ScoreChildrenUsingStateGrid" => sub {
	
	plan tests => 2;
			
	my %AtomNameHash = map { $_->getName => $_ } @$Atoms;
	
	my $O6LookUpTable = LookUpTable->New;

	$O6LookUpTable->Initiate($AtomNameHash{'O6'});
	
	my $O2LookUpTable = LookUpTable->New;

	$O2LookUpTable->Initiate($AtomNameHash{'O2'});
	
	my ($O6Node) = grep { $_->getLevel == 0 } @{$O6LookUpTable->getNodes};
	
	my ($O2Node) = grep { $_->getLevel == 0 } @{$O2LookUpTable->getNodes};
	
	is(ScoreChildrenUsingStateGrid($O6Node, $O2Node, 0), 1, 'Worked Correctly Level 0');
	
	is(ScoreChildrenUsingStateGrid($O6Node, $O2Node, 1), 3, 'Worked Correctly Level 1');
	
	
};


