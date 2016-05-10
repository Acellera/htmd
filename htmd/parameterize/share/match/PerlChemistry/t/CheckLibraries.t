#!/usr/bin/perl -w

use strict;
use warnings;
use Carp;


my $PerlChemistryPath = $ENV{'PerlChemistry'};

croak("PerlChemistry path is not set:\nbash: export PerlChemistry=/Clone/PerlChemistry") unless defined $PerlChemistryPath;

use lib $ENV{'PerlChemistry'} . "/lib";

use BaseObject ':vars';

use Parameters;

use Test::More tests => 1;

my $DefaultParameters = Parameters->New;

$DefaultParameters->Initiate;

$Parameters = $DefaultParameters;

print $Parameters->getExitifNotInitiated . "\n";

is(system("$PerlChemistryPath/t/Atom.t"), 0, 'Atoms Package Loaded Correctly');

is(system("$PerlChemistryPath/t/Bond.t"), 0, 'Bond Package Loaded Correctly');
