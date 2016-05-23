#!/usr/bin/perl -w

use strict;
use warnings;

use lib '../lib';
use Test::More tests => 1;

BEGIN { use_ok('MoleculeFileHandler') };

use MoleculeFileHandler ':all';

my $MoleculeFileHandler = MoleculeFileHandler->New("TestMol2Files/multi_load.mol2");

