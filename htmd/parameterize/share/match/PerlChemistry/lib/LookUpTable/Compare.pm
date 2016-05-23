package LookUpTable::Compare;

=head1 NAME


=head1 SYNOPSIS


=head1 DESCRIPTION


=head2 EXPORT


=head1 AUTHOR

Joseph Yesselman, E<lt>jyesselm@umich.eduE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2010 by Joseph Yesselman

This library is free software; you can redistribute it and/or modifycd
it under the same terms as Perl itself, either Perl version 5.10.0 or,
at your option, any later version of Perl 5 you may have available.

=head1 FUNCTIONS

=cut

#use 5.010000;
use strict;
use warnings;
use Carp;

require Exporter;

our @ISA = qw(Exporter);

our @EXPORT = qw(CompareLookUpTables AreAllNodesShared CompareTypes);

our $VERSION = '0.01';

#Exported Variables


#Non Exported Variables


# Preloaded methods go here.

=head2 CompareLookUpTables

Usage: $LookUpTableInstance->DetermineConsensusLookUpTable(@LookUpTableInstances);

Arguments:
  @LookUpTableInstances: All LookUpTables that need to be compared to $LookUpTableInstance

Synopsis:
  Create a new LookUpTable that is the consensus of all LookUpTables included in Arguements

=cut

sub CompareLookUpTables {
	
	my $CompareToLookUpTables = shift;
	
	my $ConsensusStateLookUpByAtomName;
	my $MatchCount = 0;
			
  #Reset all LookUpTables, this is important, pretty much the way Nodes are selected to be part of the consensus is by flagging the Used internal
  #  variable 
  $_->Reset() foreach (@$CompareToLookUpTables);

  my $ReferenceLookUpTable = shift @$CompareToLookUpTables;

  #Start with the HeadNode, which is always level 0, thus sorting by Level should put the HeadNode in the first element 
  my @ReferenceNodes = sort { $a->getLevel <=> $b->getLevel } @{$ReferenceLookUpTable->getNodes};

  $_->getAtom->setMatchNum(-1) foreach (@ReferenceNodes);

	my $NameToNodeHash = { map { $_->getAtomUniqueString => $_ }  @ReferenceNodes };

  my @NotProcessedYet = ([undef, $ReferenceNodes[0] ]);
 
  #This whole section is a Breath First search starting at the HeadNode of the ReferenceLookUpTable to find all Children that are shared across all the
  #  LookUpTables, Every time a Node is found to exist in all LookUpTables its Children are added to NotProcessedYet to also be considered.
  while(@NotProcessedYet != 0) {
	
	  my ($ParentNode, $CurrentNode) = @{shift @NotProcessedYet};
		
	  my $CurrentLevel = $CurrentNode->getLevel;
	  my $NotMatched = 0;
		
	  next if $CurrentNode->getUsed != 0; #If $CurrentNode has already been used, go to the next one
	  
	#  print $CurrentNode->getAtomName . " " . $CurrentNode->getString . "\n";
	  
	
	  my $CurrentSharedNodes = FindMatchingNodesWithStringMatch($ParentNode, $CurrentNode, $CompareToLookUpTables);
	
	  
	   foreach my $MatchingNodes (@$CurrentSharedNodes) {

		#	  print join(" ", map { $_->getAtomName} @$MatchingNodes) . "\n";

		  }
	
	  next unless @$CurrentSharedNodes;
	
	  #Find all Nodes that have the same path as $CurrentNode but also have the same Parent
	  my @CurrentNodeTwins = $ReferenceLookUpTable->FindSiblingNodesWithSameState($ParentNode, $CurrentNode);
	
	  my @ReferenceNodesWithSamePathandParent = ($CurrentNode, @CurrentNodeTwins);
	
	  my $MinAvailableSpots = @ReferenceNodesWithSamePathandParent;
	
	  foreach my $MatchingNodes (@$CurrentSharedNodes) { $MinAvailableSpots = @$MatchingNodes if $MinAvailableSpots > @$MatchingNodes }
		
		my $SortedKeysAndValues = GenerateBestNodeMatches(\@ReferenceNodesWithSamePathandParent, $CurrentSharedNodes, $NameToNodeHash, $MinAvailableSpots);

		for(my $i = 0; $i < $MinAvailableSpots; $i++) {
	 
	    my $ReferenceNodeName =shift @{$SortedKeysAndValues->[$i]};
			
			my @AssociatedNodes = @{$SortedKeysAndValues->[$i]};
			my $Node = $NameToNodeHash->{$ReferenceNodeName};
									
			$_->setUsed($Node) foreach ($Node, @AssociatedNodes);
			
			$Node->getAtom->setMatchNum($MatchCount);
			
			$MatchCount++;
										
			push(@NotProcessedYet, map { [$Node, $_] }@{$Node->getChildren});
		  		  
			$ConsensusStateLookUpByAtomName->{$Node->getAtom->getMatchNum} = getConsensusAtomInformation([$Node, @AssociatedNodes]);
    }
	
	}

	
  return $ReferenceLookUpTable->Copy unless keys %$ConsensusStateLookUpByAtomName;
			
	return $ReferenceLookUpTable->Copy($ConsensusStateLookUpByAtomName);
	
}

=head2 CompareLookUpTables

Usage: AreAllNodesShared($CompareToLookUpTables);

Arguments:
  $CompareToLookUpTables: All LookUpTables that need to be compared to $LookUpTableInstance

Synopsis:
  Create a new LookUpTable that is the consensus of all LookUpTables included in Arguements

=cut

sub AreAllNodesShared {
	
	my $CompareToLookUpTables = shift;
		
  #Reset all LookUpTables, this is important, pretty much the way Nodes are selected to be part of the consensus is by flagging the Used internal
  #  variable 
  $_->Reset() foreach (@$CompareToLookUpTables);

  my $ReferenceLookUpTable = shift @$CompareToLookUpTables;

  #Start with the HeadNode, which is always level 0, thus sorting by Level should put the HeadNode in the first element 
  my @ReferenceNodes = sort { $a->getLevel <=> $b->getLevel } @{$ReferenceLookUpTable->getNodes};

	my $NameToNodeHash = { map { $_->getAtomUniqueString => $_ }  @ReferenceNodes };

  my @NotProcessedYet = ([undef, $ReferenceNodes[0] ]);
 
  my $MatchedCount = 0;

  #This whole section is a Breath First search starting at the HeadNode of the ReferenceLookUpTable to find all Children that are shared across all the
  #  LookUpTables, Every time a Node is found to exist in all LookUpTables its Children are added to NotProcessedYet to also be considered.
  while(@NotProcessedYet != 0) {
	
	  my ($ParentNode, $CurrentNode) = @{shift @NotProcessedYet};
		
	  my $CurrentLevel = $CurrentNode->getLevel;
	  my $NotMatched = 0;
	
	  next if $CurrentNode->getUsed != 0; #If $CurrentNode has already been used, go to the next one
	  
		#print "Current: " . $CurrentNode->getAtomName . " " . $CurrentNode->Stringify .  "\n";
	 
	
	  #print "USED: " . join(" ", map { $_->getAtomName . " " . $_->Stringify } grep { $_->getUsed != 0} @ReferenceNodes) . "\n";
			
	  my $CurrentSharedNodes = FindMatchingNodesWithFullMatch($ParentNode, $CurrentNode, $CompareToLookUpTables);
	  
	  	
	  foreach my $MatchingNodes (@$CurrentSharedNodes) {
		
		 # print join(" ", map { $_->getAtomName .   " " . $_->Stringify . " " } @$MatchingNodes) . "\n";
		
	  }
	
	  #return $MatchedCount*2 / (@ReferenceNodes + @{$CompareToLookUpTables->[0]->getNodes}) unless @$CurrentSharedNodes;
	  return 0 unless @$CurrentSharedNodes;
	    
	
	  #Find all Nodes that have the same path as $CurrentNode but also have the same Parent
	  my @CurrentNodeTwins = $ReferenceLookUpTable->FindSiblingNodesWithSameFullState($ParentNode, $CurrentNode);
	
	  my @ReferenceNodesWithSamePathandParent = ($CurrentNode, @CurrentNodeTwins);
	
	  my $MinAvailableSpots = @ReferenceNodesWithSamePathandParent;
	
	  foreach my $MatchingNodes (@$CurrentSharedNodes) { $MinAvailableSpots = @$MatchingNodes if $MinAvailableSpots > @$MatchingNodes }
		
		my $SortedKeysAndValues = GenerateBestNodeMatches(\@ReferenceNodesWithSamePathandParent, $CurrentSharedNodes, $NameToNodeHash, $MinAvailableSpots);

    return 0 unless $MinAvailableSpots == @ReferenceNodesWithSamePathandParent;
 
		for(my $i = 0; $i < $MinAvailableSpots; $i++) {
	 
	    my $ReferenceNodeName =shift @{$SortedKeysAndValues->[$i]};
			
			my @AssociatedNodes = @{$SortedKeysAndValues->[$i]};
			my $Node = $NameToNodeHash->{$ReferenceNodeName};
														
			$_->setUsed($Node) foreach ($Node, @AssociatedNodes);
					
			push(@NotProcessedYet, map { [$Node, $_] }@{$Node->getChildren});
			
			$MatchedCount++;
    
    }

    #@NotProcessedYet = sort { $a->[1]->getLevel <=> $b->[1]->getLevel } @NotProcessedYet;	

  }
			
	return 1;
	
}


sub CompareTypes {
	
	my $CompareToLookUpTables = shift;
	
	my $ConsensusStateLookUpByAtomName;
	my $Score = 0;
		
  #Reset all LookUpTables, this is important, pretty much the way Nodes are selected to be part of the consensus is by flagging the Used internal
  #  variable 
  $_->Reset() foreach (@$CompareToLookUpTables);

  my $ReferenceLookUpTable = shift @$CompareToLookUpTables;

  #Start with the HeadNode, which is always level 0, thus sorting by Level should put the HeadNode in the first element 
  my @ReferenceNodes = sort { $a->getLevel <=> $b->getLevel } @{$ReferenceLookUpTable->getNodes};

	my $NameToNodeHash = { map { $_->getAtomUniqueString => $_ }  @ReferenceNodes };

  my @NotProcessedYet = ([undef, $ReferenceNodes[0] ]);
 
  my $MatchCount = 0;

  #This whole section is a Breath First search starting at the HeadNode of the ReferenceLookUpTable to find all Children that are shared across all the
  #  LookUpTables, Every time a Node is found to exist in all LookUpTables its Children are added to NotProcessedYet to also be considered.
  while(@NotProcessedYet != 0) {
	
	  my ($ParentNode, $CurrentNode) = @{shift @NotProcessedYet};
		
	  my $CurrentLevel = $CurrentNode->getLevel;
	  my $NotMatched = 0;
	
	  next if $CurrentNode->getUsed != 0; #If $CurrentNode has already been used, go to the next one
	  
	#	print "Current: " . $CurrentNode->getAtomName . " " . $CurrentNode->getAtom->Stringify .  "\n";
	 
	
	 # print "USED: " . join(" ", map { $_->getAtomName . " " . $_->Stringify } grep { $_->getUsed != 0} @ReferenceNodes) . "\n";
			
	  my $CurrentSharedNodes = FindMatchingNodesWithHybidMatch($ParentNode, $CurrentNode, $CompareToLookUpTables);
	  
	  	
	  foreach my $MatchingNodes (@$CurrentSharedNodes) {
		
		#  print join(" ", map { $_->getAtomName .   " " . $_->Stringify . " " } @$MatchingNodes) . "\n";
		
	  }
	
	  #return $MatchedCount / @ReferenceNodes unless @$CurrentSharedNodes;
	  next unless @$CurrentSharedNodes;
	    
	
	  #Find all Nodes that have the same path as $CurrentNode but also have the same Parent
	  my @CurrentNodeTwins = $ReferenceLookUpTable->FindSiblingNodesWithSameFullState($ParentNode, $CurrentNode);
	
	  my @ReferenceNodesWithSamePathandParent = ($CurrentNode, @CurrentNodeTwins);
	
	  my $MinAvailableSpots = @ReferenceNodesWithSamePathandParent;
	
	  foreach my $MatchingNodes (@$CurrentSharedNodes) { $MinAvailableSpots = @$MatchingNodes if $MinAvailableSpots > @$MatchingNodes }
		
		my $SortedKeysAndValues = GenerateBestNodeMatches(\@ReferenceNodesWithSamePathandParent, $CurrentSharedNodes, $NameToNodeHash, $MinAvailableSpots);

    return 0 unless $MinAvailableSpots == @ReferenceNodesWithSamePathandParent;
 
		for(my $i = 0; $i < $MinAvailableSpots; $i++) {
	 
	    my $ReferenceNodeName =shift @{$SortedKeysAndValues->[$i]};
			
			my @AssociatedNodes = @{$SortedKeysAndValues->[$i]};
			my $Node = $NameToNodeHash->{$ReferenceNodeName};
														
			$_->setUsed($Node) foreach ($Node, @AssociatedNodes);
					
			push(@NotProcessedYet, map { [$Node, $_] }@{$Node->getChildren});
			
			$Node->getAtom->setMatchNum($MatchCount);
			
			#$ConsensusStateLookUpByAtomName->{$Node->getAtom->getMatchNum} = getConsensusAtomInformation([$Node, @AssociatedNodes]);
			
			my $ConsensusHash = getConsensusAtomInformation([$Node, @AssociatedNodes]);
			
			if(!($Node->getNoRing == 1 && $AssociatedNodes[0]->getRingAtom == 1) || !($Node->getRingAtom == 1 && $AssociatedNodes[0]->getNoRing == 1)) { 
			
			  $Score += length($Node->Stringify) >= length($AssociatedNodes[0]->Stringify) ? length($AssociatedNodes[0]->Stringify) : length($Node->Stringify);
			
			}
			
			else {
			
			  my %NodeRings; 
			
		  	foreach my $NodeRing (@{$Node->getRings}) { $NodeRings{$NodeRing}++ }
		  	foreach my $AssociatedRing (@{$AssociatedNodes[0]->getRings}) { $NodeRings{$AssociatedRing}++ }
			
			  my %ConsensusRings; 
			
			  foreach my $ConsensusRing (@{$ConsensusHash->{'Rings'}}) { $ConsensusRings{$ConsensusRing}++ }
			
		  	while(my ($RingValue, $RingCount) = each %NodeRings) {
				
				  if(exists $ConsensusRings{$RingValue}) {
					
				    my $Diff = $ConsensusRings{$RingValue}*2 - $RingCount;
				
				    $Score -= $Diff;
					
			  	}
			 	
			  	else {
					
				  	$Score -= $RingCount;
				  
					
				  }
				
		  	}
		  }
			
			$MatchCount++;
    
    }

    #@NotProcessedYet = sort { $a->[1]->getLevel <=> $b->[1]->getLevel } @NotProcessedYet;	

  }
	  
	
	#return 0 unless keys %$ConsensusStateLookUpByAtomName;

	#return $ReferenceLookUpTable->Copy($ConsensusStateLookUpByAtomName);
		
			
  #return $MatchedCount*2 / (@ReferenceNodes + @{$CompareToLookUpTables->[0]->getNodes});

  return $Score;

}


=head2 GenerateBestNodeMatches

Usage: LookUpTable->GenerateBestNodeMatches($ReferenceNodes, $CurrentSharedNodes, $NameToNodeHash, MinAvailableSpots);

Arguments:

Synopsis:

=cut

sub GenerateBestNodeMatches {
	
	my ($ReferenceNodes, $CurrentSharedNodes, $NameToNodeHash, $MinAvailableSpots) = @_;
	
	my %TrackCorrectMatches = map { $_->getAtomUniqueString => [] } @$ReferenceNodes;
	
	my $NextLevel = $ReferenceNodes->[0]->getLevel + 1;
		
	my $Redo = 0; my @SortedKeysAndValues; my $Finished = 0;
		
	while(!$Finished) {
					
	  $Redo = 0;
    
    foreach my $MatchingNodes (@$CurrentSharedNodes) {
																											
	    my @Results = DetermineBestMatchUps($NextLevel, $ReferenceNodes, $MatchingNodes);
			
			push @{$TrackCorrectMatches{$_->[0]}}, $_->[1] foreach (@Results);
			
	  }
		
    @SortedKeysAndValues = sort { @$b <=> @$a } map { [$_, @{$TrackCorrectMatches{$_}}] } keys(%TrackCorrectMatches);

    foreach (0 .. $MinAvailableSpots - 1) { $Redo = 1 if @{$SortedKeysAndValues[$_]} != @$CurrentSharedNodes + 1  }
			
	  if(!$Redo) { $Finished = 1; last }
			
	  $ReferenceNodes = [ map { $NameToNodeHash->{$SortedKeysAndValues[$_]->[0]} } (0 .. $MinAvailableSpots - 1) ];

    $TrackCorrectMatches{$_} = [ ] foreach (keys %TrackCorrectMatches);

  }

  return \@SortedKeysAndValues;
	
	
}

=head2 FindMatchingNodesWithFullMatch

Usage: LookUpTable->FindMatchingNodesWithFullMatch($ParentNode, $CurrentNode, $CompareToLookUpTables);

Arguments:

Synopsis:

=cut

sub FindMatchingNodesWithFullMatch($$$) {
	
	my ($ParentNode, $CurrentNode, $CompareToLookUpTables) = @_;
	
  #In every LookUpTable see if there is are Nodes that has the same Path to the Head Node as $CurrentNode, collect all that do, If there is a LookUpTable
  #  that does not have any Nodes that match $CurrentNode than set $NotMatched to 1 and exit out, this means that it is not possible to include $CurrentNode
	#  in the consensus since it is a LookUpTable is missing it
	
	my @CurrentSharedNodes;  
	
	foreach my $CompareToLookUpTable (@$CompareToLookUpTables) {
		
	  my @CompareToCurrentLevelNodes = $CompareToLookUpTable->FindCorrespondingChildNodes($ParentNode);
	
	  #print "USED: " . join(" ", map { $_->getAtomName . " " . $_->Stringify  } grep { $_->getUsed != 0} @{$CompareToLookUpTable->getNodes}) . "\n";
	
	  #print "NODES: " . join(" ", map {$_->getAtomName . " " . $_->Stringify } @CompareToCurrentLevelNodes) . "\n";
		  
		my @MatchingNodes = grep { $CurrentNode->FullMatchesTo($_) } @CompareToCurrentLevelNodes;
		  
		#print "MATCHING: " . join(" ", map {$_->getAtomName . " " . $_->Stringify } @MatchingNodes) . "\n";
	  
		
		if(@MatchingNodes == 0) { return undef }
		
		push @CurrentSharedNodes, \@MatchingNodes;
	}
		
	return \@CurrentSharedNodes;
	
	
}

=head2 FindMatchingNodesWithStringMatch

Usage: LookUpTable->FindMatchingNodesWithStringMatch($ParentNode, $CurrentNode, $CompareToLookUpTables);

Arguments:

Synopsis:

=cut

sub FindMatchingNodesWithStringMatch($$$) {
	
	my ($ParentNode, $CurrentNode, $CompareToLookUpTables) = @_;
	
  #In every LookUpTable see if there is are Nodes that has the same Path to the Head Node as $CurrentNode, collect all that do, If there is a LookUpTable
  #  that does not have any Nodes that match $CurrentNode than set $NotMatched to 1 and exit out, this means that it is not possible to include $CurrentNode
	#  in the consensus since it is a LookUpTable is missing it
	
	my @CurrentSharedNodes;  
	
	foreach my $CompareToLookUpTable (@$CompareToLookUpTables) {
		
	  my @CompareToCurrentLevelNodes = $CompareToLookUpTable->FindCorrespondingChildNodes($ParentNode);
			  
		my @MatchingNodes = grep { $CurrentNode->StringMatchesTo($_) } @CompareToCurrentLevelNodes;
		  
	 # print "NODES: " . join(" ", map {$_->getAtomName . " " . $_->Stringify . " " . $_->getAtom->getMolecule->getName . " " } @CompareToCurrentLevelNodes) . "\n";
	  
		
		if(@MatchingNodes == 0) { return undef }
		
		push @CurrentSharedNodes, \@MatchingNodes;
	}
		
	return \@CurrentSharedNodes;
	
	
}

=head2 FindMatchingNodesWithFullMatch

Usage: LookUpTable->FindMatchingNodesWithFullMatch($ParentNode, $CurrentNode, $CompareToLookUpTables);

Arguments:

Synopsis:

=cut

sub FindMatchingNodesWithHybidMatch($$$) {
	
	my ($ParentNode, $CurrentNode, $CompareToLookUpTables) = @_;
	
  #In every LookUpTable see if there is are Nodes that has the same Path to the Head Node as $CurrentNode, collect all that do, If there is a LookUpTable
  #  that does not have any Nodes that match $CurrentNode than set $NotMatched to 1 and exit out, this means that it is not possible to include $CurrentNode
	#  in the consensus since it is a LookUpTable is missing it
	
	my @CurrentSharedNodes;  
	
	foreach my $CompareToLookUpTable (@$CompareToLookUpTables) {
		
	  my @CompareToCurrentLevelNodes = $CompareToLookUpTable->FindCorrespondingChildNodes($ParentNode);
	
	 # print "USED: " . join(" ", map { $_->getAtomName . " " . $_->Stringify  } grep { $_->getUsed != 0} @{$CompareToLookUpTable->getNodes}) . "\n";
	
	  #print "NODES: " . join(" ", map {$_->getAtomName . " " . $_->Stringify } @CompareToCurrentLevelNodes) . "\n";
		  
		my @MatchingNodes = grep { $CurrentNode->FullMatchesTo($_) } @CompareToCurrentLevelNodes;
		
		if(@MatchingNodes == 0 ) {
			
			 @MatchingNodes = grep { $CurrentNode->StringMatchesTo($_) } @CompareToCurrentLevelNodes;
			
		} 
		  
		#print "MATCHING: " . join(" ", map {$_->getAtomName . " " . $_->Stringify } @MatchingNodes) . "\n";
	  
		
		if(@MatchingNodes == 0) { return undef }
		
		push @CurrentSharedNodes, \@MatchingNodes;
	}
		
	return \@CurrentSharedNodes;
	
	
}

=head2 DetermineBestMatchUps

Usage: LookUpTable->DetermineBestMatchUps($NextLevel, $ReferenceNodes, $MatchingNodes);

Arguments:

Synopsis:

=cut

sub DetermineBestMatchUps { 
	
  my ($NextLevel, $ReferenceNodes, $MatchingNodes) = @_;

  my $SizeOfMatchingNodes = @$MatchingNodes;
	my $SizeOfReferenceNodes = @$ReferenceNodes;

	my $Min = $SizeOfMatchingNodes <= $SizeOfReferenceNodes ? $SizeOfMatchingNodes : $SizeOfReferenceNodes;
	
	$_->UpdateChildScoreGrid foreach (@$ReferenceNodes, @$MatchingNodes);
	
	my @MatchUps;
	
	foreach my $ReferenceNode (@$ReferenceNodes) {
	  push (@MatchUps, map { [$ReferenceNode, $_, 0] } @$MatchingNodes);
	}	
	
	my $Recurse = sub { 
		
		my $List = shift;
		my $Allowed = shift;
		
		my @NextAllowed = grep {$_ % $SizeOfMatchingNodes != $List->[-1] % $SizeOfMatchingNodes && int($_/$SizeOfMatchingNodes) != int($List->[-1]/$SizeOfMatchingNodes) } @$Allowed;
	 
	  return [$List, \@NextAllowed] if ! @NextAllowed;
	
	  my @NewLists = map { [ @$List, $_] } @NextAllowed;
	
	  return map { [$_, \@NextAllowed] } @NewLists;
		
	};
	
	my %ResultsHash;

	foreach my $Start (0 .. $#MatchUps) {

		my @Results = &$Recurse([$Start], [0 .. $#MatchUps]);

		my $Count = 1;

		while($Count < $Min) {
			@Results = map { &$Recurse(@{$_}) } @Results;
			$Count++;
		}

		foreach (@Results) {
			
			my $StringedResult = join(" ", sort @{$_->[0]});
			
		  $ResultsHash{$StringedResult} = $_->[0] if ! exists $ResultsHash{$StringedResult};
	  }
	
	}
	
	my @Permutations = map { [[map { \$_ } @MatchUps[@{$_}]], 0] } values %ResultsHash;
	
	my $PermutationNumber = 0;
	my $BestScore = 0;
	my $Score = 0;
	
	$PermutationNumber = $#Permutations;
	
  
  while( $PermutationNumber > 0) {	
		
		$_->[2] = ScoreChildrenUsingStateGrid($_->[0], $_->[1], $NextLevel) foreach (@MatchUps);
		  		
		foreach my $Permutation (@Permutations) {
			$Score = 0;
			$Score += $$_->[2] foreach (@{$Permutation->[0]});
			$BestScore = $Score if $BestScore < $Score;
			$Permutation->[1] = $Score;
		}
		
		my $PermutationPointer = 0;
		
		while ( $PermutationPointer <= $PermutationNumber) {
			
		  if($Permutations[$PermutationPointer]->[1] < $BestScore) { $PermutationNumber--; splice(@Permutations, $PermutationPointer, 1) }
		  else																						         { $PermutationPointer++ }
			
		}
		
		$NextLevel++;
				
		$PermutationNumber = 0 if $NextLevel >= 10 && $PermutationNumber != 0;
				 			
	}
		
	return map { [$$_->[0]->getAtomUniqueString, $$_->[1]] } @{$Permutations[0]->[0]}
	
}

=head2 getConsensusAtomInformation

Usage: LookUpTable->getConsensusAtomInformation($Nodes)

Arguments:
  $Nodes: Nodes that matched to each other

Synopsis:
  Forms a consensus from the Atoms that $Nodes store, this is for the state of the atom, if there is rings or nonring
 
=cut

sub getConsensusAtomInformation {
	
	my $Nodes = shift;
	
	my @NodesWithRings = grep { @{$_->getAtomRings} != 0 } @$Nodes;
		
	return { NoRing   => @NodesWithRings == 0 ? $Nodes->[0]->getAtomElement !~ /CL|BR|H|F|I/ : 0, 
		       Rings    => [], 
		       RingAtom => 0,
		       String   => $Nodes->[0]->getString } if @NodesWithRings != @$Nodes;
	
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
		
	return { NoRing   => 0, 
		       Rings    => \@ConsensusRings,
		       RingAtom => 1,
		       String   => $Nodes->[0]->getString };
}

=head2 ScoreChildrenUsingStateGrid 

Usage: LookUpTable->ScoreChildrenUsingStateGrid($Node1, $Node2, $Level);

Arguments:

Synopsis:
  Scores how related 2 nodes are that have the same state

=cut

sub ScoreChildrenUsingStateGrid($$$) {
	
	my ($Node1, $Node2, $Level) = @_;
	
	my $StateGridForNode1 = $Node1->getChildScoreGrid;
  my $StateGridForNode2 = $Node2->getChildScoreGrid;
	my $Score = 0;
	
	while(my ($State, $GridElement) = each %$StateGridForNode1) {
	
	  my $StateLevel = substr($State, -2);
	
	  next if $StateLevel > $Level || ! exists $StateGridForNode2->{$State};	  
			
    $Score += $GridElement->getCount <= $StateGridForNode2->{$State}->getCount ? $GridElement->getCount : $StateGridForNode2->{$State}->getCount;	
	
	}
	
	return $Score;
	
}


1;
__END__
