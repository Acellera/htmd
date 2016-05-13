#!/usr/bin/env perl

use strict;
use warnings;

use lib '../lib';
use Test::More tests => 3;

BEGIN { use_ok('BaseObject') };

use BaseObject ':all';

ok(defined $Parameters, 'can load $Parameters');

#Make dummy internal variable to check to see if I can retreive it through AUTOLOAD
my $Self = { _Test => 1};
my $test = bless $Self, 'BaseObject';

ok($test->getTest == 1, 'AUTOLOAD function is operation');

#ok(PrintMessage("Love", 0), 'test PrintMessage');

#########################

# Insert your test code below, the Test::More module is use()ed here so read
# its man page ( perldoc Test::More ) for help writing this test script.

