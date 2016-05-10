#!/usr/bin/env perl

use strict;
use warnings;

use Test::More tests => 1;
use lib $ENV{'MATCH'} . "/lib"; 

BEGIN { use_ok('MATCHParameters') };

my $Parameters = MATCHParameters->New;

$Parameters->Initiate;

print $Parameters->getBondLengthFilePath . "\n";
