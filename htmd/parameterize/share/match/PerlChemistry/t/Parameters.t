#!/usr/bin/perl -w

use strict;
use warnings;

use lib '../lib';
use Test::More tests => 2;

BEGIN { use_ok('Parameters') };

use Parameters ':all';

my $Parameters = Parameters->New;

ok(exists $ENV{'PerlChemistry'}, 'PerlChemistry is defined');

$Parameters->Initiate;

print $Parameters->getBondLengthFilePath . "\n";

