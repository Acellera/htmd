#!/usr/bin/perl -w

use strict;
use warnings;

use lib '../../lib';
use Test::More tests => 4;

BEGIN { use_ok('LookUpTable::Node') };

use Atom;


sub TESTgetConsensusAtomInformation {
	
	my $Nodes = shift;
	
	my @NodesWithRings = grep { $_->getAtomRings } @$Nodes;
		
	return { NoRing => @NodesWithRings == 0 ? $Nodes->[0]->getAtomElement !~ /CL|BR|H|F|I/ : 0, 
		       Rings  => [], 
		       String => $Nodes->[0]->getString } if @NodesWithRings != @$Nodes;
	
	my (%SeenRings, %SeenNodes, @ConsensusRings);
	
	foreach my $NodeWithRing (@NodesWithRings) {
	
	  foreach my $Ring (@{$NodeWithRing->getAtomRings}) {
		
		  $SeenRings{$Ring}++;
		  $SeenNodes{$Ring}{$NodeWithRing}++;
		
	  }
	
	}
	
	while ( my ($Ring, $Count) = each %SeenRings) {
				
		next if $Count < @NodesWithRings || (keys %{$SeenNodes{$Ring}}) != @$Nodes;
				
	  push @ConsensusRings, $Ring foreach ( 0 .. int ( $Count / @$Nodes ) - 1);
		
	}
	
	return { NoRing => 0, 
		       Rings  => \@ConsensusRings,
		       String => $Nodes->[0]->getString };	
	
}

my $Atom = Atom->New({ State => "C.3", Rings => [6], Element => 'C' });

my $Node1 = LookUpTable::Node->New($Atom);

my $Results = TESTgetConsensusAtomInformation([$Node1]);

is_deeply($Results->{Rings}, [6], 'Rings worked');

is($Results->{NoRing}, 0, 'NoRing Worked');

my $Atom2 = Atom->New({ State => "C.3", Rings => [5, 5], Element => 'C' });

my $Node2 = LookUpTable::Node->New($Atom2);

$Results = TESTgetConsensusAtomInformation([$Node1, $Node2]);

is_deeply($Results->{Rings}, [], 'Rings worked');

is($Results->{NoRing}, 0, 'NoRing Worked');







