#!/usr/bin/env perl

use strict;
use warnings;

use Test::More tests => 1;
use lib $ENV{'MATCH'} . "/lib"; 
use Carp;

use MATCHBaseObject;
use BaseObject ':vars';

use MATCHParameters;
use MoleculeFileHandler;

BEGIN { use_ok('MATCHer') };

#Deal With Initial Parameter Setup, This will be handled in main .pl script
my $DefaultParameters = MATCHParameters->New;

$DefaultParameters->Initiate;

$Parameters = $DefaultParameters;

my $MATCHer = MATCHer->New;

$MATCHer->Initiate;

#Grab Test Molecule

my $MoleculeFileHandler = MoleculeFileHandler->New("MoleculeParameterizerResources/1bjy_ligand_correct.pdb");

my $Structures = $MoleculeFileHandler->BuildObjectsFromFile;

my $TestMolecule = $Structures->[0]->getMolecule(0);

$TestMolecule->Initiate;
$TestMolecule->setName("UNK");

$MATCHer->BuildBothTopologyAndParameterFileForMolecule($TestMolecule);