#!/usr/bin/perl -w

use strict;
use warnings;

use lib '../../lib';
use Test::More tests => 9;

BEGIN { use_ok('LookUpTable::GridElement') };

use LookUpTable::GridElement;
use LookUpTable::Node;
use Atom;

my $Atom = Atom->New( { State => "C.3", Level => 0 });

my $Node = LookUpTable::Node->New($Atom);
 
$Node->setUsed(1);

my $GridElement = LookUpTable::GridElement->New;

$GridElement->AddNodeReference(\$Node);

is($GridElement->getCount, 1, 'Added correctly');

$GridElement->RemoveNodeReference(\$Node);

is($GridElement->getCount, 0, 'Deleted correctly');

is_deeply($GridElement->getNodeReferences, [ ], 'Everything cleaned out correctly');

$GridElement->AddNodeReference(\$Node);
$GridElement->AddNodeReference(\$Node);

my $CopiedGridElement = $GridElement->Copy;

is_deeply($GridElement->getNodeReferences, $CopiedGridElement->getNodeReferences, 'Copied Correctly');

is($GridElement->getCount, $CopiedGridElement->getCount, 'Count the Same');

$GridElement->Merge($CopiedGridElement);

is($GridElement->getCount, 1, 'Count Updated With Merge');

$GridElement->Update;

is($GridElement->getCount, 0, 'Count is Updated');

$Node->setUsed(undef);

$GridElement->Update;

is($GridElement->getCount, 1, 'Count is Updated');
