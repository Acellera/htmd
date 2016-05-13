package LookUpTable;

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
use BaseObject;
use LookUpTable::Node;
use LookUpTable::GridElement;
use LookUpTable::Compare;

require Exporter;

our @ISA = qw(BaseObject Exporter);

our @EXPORT = qw(CopyAtomsInNodes);

our $VERSION = '0.01';

#Exported Variables


#Non Exported Variables
our $HARD_DEPTH_LIMIT = 10;


# Preloaded methods go here.

=head2 New

Usage: LookUpTable->New;

Arguments:
  $Class: should be 'LookUpTable'

Synopsis:
  creates a new LookUpTable object

=cut

sub New {
	
	my $Class = shift;
		
	my $Self = { 
	  _Nodes              => undef,
	  _LengthOfStringify	=> 0, 
	  _Name								=> 0
	};
	
	bless $Self, $Class;
  return $Self;
	
}

sub AddNode {
	
	my ($Self, $NewNode) = @_;
	
	my $Nodes = $Self->getNodes;
	
	#Already included in LookUpTable
	return if grep { $_->getAtomName eq $NewNode->getAtomName} @$Nodes;
		
	my @Parents;
	
	foreach my $BondedTo (@{$NewNode->getAtom->getBondedAtoms}) { 
		
		foreach (@$Nodes) { push @Parents, $_ if $BondedTo->getName eq $_->getAtomName && $_->getLevel + 1 == $NewNode->getLevel }
		
	}

  confess("Could not find any parent nodes") if ! @Parents;

	foreach (@Parents) { $NewNode->AddParent($_); $_->AddChild($NewNode) }
	
	my $ChildScoreGridEntry = $NewNode->getAtomState . ":" . sprintf("%02s", $NewNode->getLevel);
	
	my $GridElement = LookUpTable::GridElement->New;
	
	$GridElement->AddNodeReference(\$NewNode);
	
	$NewNode->setChildScoreGrid( { $ChildScoreGridEntry => $GridElement });
	
	my $ReachedLevelZero = 0;
	
	my %PreviousParentHash = map { $_->getAtomName => $ChildScoreGridEntry } @Parents;
	
	my @NextLevel;
		
	while($ReachedLevelZero) {
		
		@NextLevel = ();
		
		foreach my $Parent (@Parents) {
			
			my $String = $Parent->getAtomState . "-" . $PreviousParentHash{$Parent->getAtomName};
			
			if(exists $Parent->getChildScoreGrid->{$String}) { $Parent->getChildScoreGrid->{$String}->Merge($GridElement) }
			else																	 		  	   { $Parent->getChildScoreGrid->{$String} = $GridElement->Copy }
			
	    $Parent->getChildScoreGrid->{$String}++;
			
			if($Parent->Level == 0) { $ReachedLevelZero = 0; last}
			
			my @NextLevelParents = @{$Parent->getParents};
			
			foreach my $NextLevelParent (@NextLevelParents) { $PreviousParentHash{$NextLevelParent->getAtomName} = $String }
			
			push(@NextLevel, @NextLevelParents);
		}
		
		my %SeenParentName;
		my @Uniq = grep { ! $SeenParentName{ $_->getAtomName }++ } @NextLevel;
		
		@Parents = @Uniq;
	}
	
  push @$Nodes, $NewNode;
	
}

=head2 AddAtom

Usage: $LookUpTableInstance->AddAtom($Atom);

Arguments:
  $Atom: the Atom object you wish to add to the LookUpTable, must be part of the same molecule!

Synopsis:
  Adds Atom to LookUpTable in the form of a new Node, must be part of the same molecule that the LookUpTable originates from

=cut

sub AddAtom  {
	
	my ($Self, $Atom) = @_;
	
	my $AtomNode = LookUpTable::Node->New($Atom);
  
  $Self->AddNode($AtomNode);	
	
}

=head2 BuildNodesFromAtoms

Usage: $LookUpTableInstance->BuildNodesFromAtoms($HeadAtom);

Arguments:
  $HeadAtom: Starting Atom, to build LookUpTable for

Synopsis:
  Builds Nodes for LookUpTable Object

=cut
sub BuildNodesFromAtoms { 
	
	my ($Self, $HeadAtom) = @_;
	
  my %VisitedAtoms;
  my @ToVisitAtoms;
  my %AtomToNodeArrayLookUp;

  $HeadAtom->setLevel(0);
  push @ToVisitAtoms, $HeadAtom;

  $AtomToNodeArrayLookUp{$ToVisitAtoms[0]->getUniqueString} = LookUpTable::Node->New($ToVisitAtoms[0]);

  my ($CurrentAtom, $CurrentNode);

  #Breath first search, visits alls Atom Objects within Range of HARD_DEPTH_LIMIT
  while( @ToVisitAtoms > 0) {
	
	  $CurrentAtom = shift @ToVisitAtoms;
	
	  next if $CurrentAtom->getLevel > $HARD_DEPTH_LIMIT;
	
	  my $AtomsBondedToCurrent = $CurrentAtom->getBondedAtoms;
	
	  $CurrentNode = $AtomToNodeArrayLookUp{$CurrentAtom->getUniqueString} if @$AtomsBondedToCurrent;
	
	  my @SortedAtomsBondedtoCurrent = sort { $b->getState cmp $a->getState } @$AtomsBondedToCurrent;
	
	  foreach my $BondedAtom (@SortedAtomsBondedtoCurrent) {
		
		  my $InToVisitAtoms = scalar(grep {$_ eq $BondedAtom } @ToVisitAtoms);
		
		  if(!$InToVisitAtoms && ! exists $VisitedAtoms{$BondedAtom->getName} ) {
			
			  $BondedAtom->setLevel($CurrentAtom->getLevel + 1);
									
			  my $ChildNode = LookUpTable::Node->New($BondedAtom);
			
			  $ChildNode->AddParent($CurrentNode); $CurrentNode->AddChild($ChildNode);
			
			  $AtomToNodeArrayLookUp{$BondedAtom->getUniqueString} = $ChildNode;
			
			  push @ToVisitAtoms, $BondedAtom;
			
		  }
		
		  elsif(! exists $VisitedAtoms{$BondedAtom->getName} ) {
			
			  my $VisitedNode = $AtomToNodeArrayLookUp{$BondedAtom->getUniqueString};
						
			  if($VisitedNode->getLevel == $CurrentNode->getLevel + 1) { $VisitedNode->AddParent($CurrentNode); $CurrentNode->AddChild($VisitedNode) }
			
			  elsif($VisitedNode->getLevel + 1 == $CurrentNode->getLevel) { $VisitedNode->AddChild($CurrentNode); $CurrentNode->AddParent($VisitedNode) }
			
			  elsif($VisitedNode->getLevel == $CurrentNode->getLevel) {
				
				  $VisitedNode->AddRelated($CurrentNode); $CurrentNode->AddRelated($VisitedNode);
				
			  }
			
			  #Should I do related seems worthless?
			
		  }
				
	  }
	
	  $VisitedAtoms{$CurrentAtom->getName} = $CurrentAtom;
	
  }

  my @Nodes = values %AtomToNodeArrayLookUp;

  $_->setLevel(0) foreach (values %VisitedAtoms);

  $Self->setNodes(\@Nodes);
	
}
=head2 CopyAtomsInNodes

Usage: LookUpTable->CopyAtomsInNodes

Arguments:

Synopsis:
  Copies Atoms from selected Nodes perserving bond connectivity

=cut

sub CopyAtomsInNodes {
	
	my $Nodes = shift;
	
	#Uses the name of Node Atoms to check to see if they have been copied yet
	my %CopiedAtoms;
	
	my $ReturnCopiedAtomOrCreateACopy = sub { 
		my $Node = shift;
				
	  return $CopiedAtoms{$Node->getAtom->getMatchNum} if exists $CopiedAtoms{$Node->getAtom->getMatchNum};
			
		my $CopiedAtom = $Node->getAtom->Copy;
		
	  $CopiedAtom->setLevel($Node->getLevel);
		
	  $CopiedAtoms{$Node->getAtom->getMatchNum} = $CopiedAtom;
		
		return $CopiedAtom;
		
	};
	
  foreach my $Node (@$Nodes) {
	
	  my $NodeAtomCopy = &$ReturnCopiedAtomOrCreateACopy($Node);
	
	  my @RelatedNodes;
		
	  foreach my $RelatedNode (@{$Node->getParents}, @{$Node->getChildren}) { 
		
		  next unless grep { $RelatedNode eq $_ } @$Nodes;
		
		  my $RelatedNodeAtomCopy =  &$ReturnCopiedAtomOrCreateACopy($RelatedNode);
		
		  Bond->New($NodeAtomCopy, $RelatedNodeAtomCopy, 1) if ! $NodeAtomCopy->IsBondedTo($RelatedNodeAtomCopy);
		  
		}
		
  }

  return values %CopiedAtoms;
}


=head2 Copy

Usage: $LookUpTableInstance->Copy

Arguments:
  $ConsensusStateLookUpByAtomName = An optional arguement given by DetermineConsensusLookUpTable that replaces Node information with 
the corresponding consensus information

Synopsis:
  Copies Atoms from selected Nodes perserving bond connectivity

=cut

sub Copy {
	
	my ($Self, $ConsensusStateLookUpByAtomName) = @_;
	
	my $Nodes = $Self->getNodes;
	my @NodesToCopy;
	 	
	if(defined $ConsensusStateLookUpByAtomName) {	@NodesToCopy = grep { $_->getUsed } @$Nodes }
	
	else {
				
		$ConsensusStateLookUpByAtomName = { };
		
		my $Count = 1;		
		
		foreach my $Node (@$Nodes) {
			
			$Node->getAtom->setMatchNum($Count++);
			
			$ConsensusStateLookUpByAtomName->{$Node->getAtom->getMatchNum} = {String => $Node->getString, NoRing => $Node->getNoRing, Rings => $Node->getRings, RingAtom => $Node->getRingAtom };
			
		}
		
		@NodesToCopy =  @$Nodes; 
	
	}
		
	#Collect all the copied atoms
	my @CopiedAtoms = CopyAtomsInNodes(\@NodesToCopy);
	
	return undef unless @CopiedAtoms;  
	
	foreach my $CopiedAtom (@CopiedAtoms) {    
		
		next unless exists $ConsensusStateLookUpByAtomName->{$CopiedAtom->getMatchNum};

		my $SavedConsensusInfo = $ConsensusStateLookUpByAtomName->{$CopiedAtom->getMatchNum};
				
		$CopiedAtom->setState($SavedConsensusInfo->{'String'}); $CopiedAtom->setNoRing($SavedConsensusInfo->{'NoRing'}); $CopiedAtom->setRings($SavedConsensusInfo->{'Rings'});		
		$CopiedAtom->setRingAtom($SavedConsensusInfo->{'RingAtom'});
		
	}	
	
	my $HeadAtom;
	
	foreach (@CopiedAtoms) { $HeadAtom = $_ if $_->getLevel == 0}
	  
  confess("There is no HeadAtom in this Set of Atoms") unless defined $HeadAtom;

  my $CopiedLookUpTable = LookUpTable->New;
	$CopiedLookUpTable->Initiate($HeadAtom);
		
	return $CopiedLookUpTable;
	
}

=head2 DetermineConsensusLookUpTable

Usage: $LookUpTableInstance->DetermineConsensusLookUpTable(@LookUpTableInstances);

Arguments:
  @LookUpTableInstances: All LookUpTables that need to be compared to $LookUpTableInstance

Synopsis:
  Create a new LookUpTable that is the consensus of all LookUpTables included in Arguements

=cut

sub DetermineConsensusLookUpTable {
	
	my ($Self, $CompareToLookUpTables) = @_;
	
	unshift @$CompareToLookUpTables, $Self;
		
  return CompareLookUpTables($CompareToLookUpTables);
}

=head2 AreAllNodesShared

Usage: $LookUpTableInstance->DetermineConsensusLookUpTable(@LookUpTableInstances);

Arguments:
  @LookUpTableInstances: All LookUpTables that need to be compared to $LookUpTableInstance

Synopsis:
  Create a new LookUpTable that is the consensus of all LookUpTables included in Arguements

=cut

sub AreAllNodesSharedBetween {
	
	my ($Self, $CompareToLookUpTables) = @_;
	
	unshift @$CompareToLookUpTables, $Self;
	
	return AreAllNodesShared($CompareToLookUpTables, undef , 1);
	
}


sub CompareTypeLookUpTables {
	
	my ($Self, $CompareToLookUpTables) = @_;
	
	unshift @$CompareToLookUpTables, $Self;
	
	return CompareTypes($CompareToLookUpTables, undef , 1);
	
}

=head2 FindCorrespondingChildNodes

Usage: $LookUpTableInstance->FindCorrespondingChildNodes($ParentNode);

Arguments:
  $ParentNode: The Parent Node from the Reference Table that one is trying to match up its children too

Synopsis:
  Finds which Node in this LookUpTable corresponds to $ParentNode from the Reference LookUpTable during DetermineConsensusLookUpTable
return its Children.
 
=cut

sub FindCorrespondingChildNodes { 
	
	my ($Self, $ParentNode) = @_;

	my $Nodes = $Self->getNodes;
	
	#No Parent Node this must be Level 0
	if(! defined $ParentNode){
		
		foreach (@$Nodes) { return $_ if $_->getLevel == 0 }
		
	}
	
	my $SelfPairedParentNode;
	
	foreach (@$Nodes) {
    if($_->getUsed != 0 && $_->getUsed eq $ParentNode) { $SelfPairedParentNode = $_; last }
	}
	
	confess("No Corresponding Parent Found") if ! defined $SelfPairedParentNode;	
			
	return grep { ! $_->getUsed } (@{$SelfPairedParentNode->getChildren}, @{$SelfPairedParentNode->getRelated});

}

=head2 FindSiblingNodesWithSameState

Usage: $LookUpTableInstance->FindSiblingNodesWithSameState($ParentNode, $Node);

Arguments:
  $ParentNode: The Parent Node from the Reference Table
  $Node: The Node that one is trying to find Siblings of

Synopsis:
  Find all Nodes with the Same state and same Parent
 
=cut

sub FindSiblingNodesWithSameState {
	my ($Self, $Parent, $Node) = @_;
	
	#$Node is the HeadNode
	return () if ! defined $Parent;
	
	my $Children = $Parent->getChildren;
	return grep { $Node->StringMatchesTo($_) && $Node ne $_ && $_->getUsed == 0} @$Children;
}

=head2 FindSiblingNodesWithSameFullState

Usage: $LookUpTableInstance->FindSiblingNodesWithSameFullState($ParentNode, $Node);

Arguments:
  $ParentNode: The Parent Node from the Reference Table
  $Node: The Node that one is trying to find Siblings of

Synopsis:
  Find all Nodes with the Same state and same Parent
 
=cut

sub FindSiblingNodesWithSameFullState {
	my ($Self, $Parent, $Node) = @_;
	
	#$Node is the HeadNode
	return () if ! defined $Parent;
	
	my $Children = $Parent->getChildren;
	return grep { $Node->FullMatchesTo($_) && $Node ne $_ && $_->getUsed == 0} @$Children;
}

=head2 Initiate

Usage: $LookUpTableInstance->Initiate($HeadAtom);

Arguments:
  $HeadAtom: Starting Atom, to build LookUpTable for

Synopsis:
  Setups LookUpTableInstance, must be called!

=cut

sub Initiate {
	
  my ($Self, $HeadAtom) = @_;

  #No Atoms!
  if(!$HeadAtom) { return undef}

	#Builds the LookUpTable using an array of Atoms which are assumed to represent a Type or Molecule
	$Self->BuildNodesFromAtoms($HeadAtom);
	
  $Self->SetupNodeConnectivity;

  $Self->Reset;
  
  $Self->setLengthOfStringify(length $Self->Stringify); 

}

=head2 RemoveNode

Usage: $LookUpTableInstance->RemoveNode($Node);

Arguments:
  $Node: the Node object you wish to remove to the LookUpTable

Synopsis:
  Removes a Node from $LookUpTableInstance
 
=cut

sub RemoveNode {
	
	my ($Self, $DeletedNode) = @_;
	
	my $Nodes = $Self->getNodes;
	
	my @Parents = @{$DeletedNode->getParents};
	
	my $Children = $DeletedNode->getChildren;
	
	confess("Can only remove nodes that have no children!") if @$Children;
	
	my $DeletedNodeGridEntry = $DeletedNode->getString . ":" . sprintf("%02s", $DeletedNode->getLevel); 
	
	foreach (@Parents) {
		
		$_->RemoveChild($DeletedNode); $DeletedNode->RemoveParent($_);
		
	}
	
	foreach (@{$DeletedNode->getRelated}) {
		
		$_->RemoveRelated($DeletedNode); $DeletedNode->RemoveRelated($_);
		
	}
	
	foreach (0 .. @$Nodes-1) { 
		
		 if($Nodes->[$_] eq $DeletedNode) { splice(@$Nodes, $_, 1); last }
		
  }

	my %PreviousParentHash = map { $_->getAtomName => $DeletedNodeGridEntry } @Parents;
	
	my @NextLevel;
	
	my $ReachedLevelZero = 0;
	
	my %SeenParent;
	
	my $ToMany = 0;
			
	while(!$ReachedLevelZero) {
		
		@NextLevel = ();
		
		foreach my $Parent (@Parents) {
			
		#	$SeenParent { $Parent->getAtomName } = 1;
			
			my $String = $Parent->getAtomState . "-" . $PreviousParentHash{$Parent->getAtomName};
			
			$Parent->getChildScoreGrid->{$String}->RemoveNodeReference(\$DeletedNode);
			
			if($Parent->getLevel == 0) { $ReachedLevelZero = 1; last}
			
			my @NextLevelParents = grep { ! $SeenParent{$_->getAtomName}  } @{$Parent->getParents};
			
			foreach my $NextLevelParent (@NextLevelParents) { $PreviousParentHash{$NextLevelParent->getAtomName} = $String }
			
			push(@NextLevel, @NextLevelParents);
		}
		
		my %SeenParentName;
		my @Uniq = grep { ! $SeenParentName{ $_->getAtomName }++ } @NextLevel;
						
		@Parents = @Uniq;
				
		$ToMany++;
		
	#	exit 1 if $ToMany > 100;
								
	}
						
}

=head2 Reset

Usage: $LookUpTableInstance->Reset;

Arguments:

Synopsis:
  Resets all Used variables for Nodes, Used is important for doing LookUpTable comparisons, must be reset before each comparison

=cut

sub Reset {
	
	my $Self = shift;
	
	my $Nodes = $Self->getNodes;
  
  $_->setUsed(0) foreach (@$Nodes);

  $_->getAtom->setMatchNum(-1) foreach (@$Nodes);


}

=head2 SetupNodeConnectivity

Usage: $LookUpTableInstance->SetAtomLevelsUsingNodeLevels;

Arguments:

Synopsis:
  Builds ChildNodeGrids, which allows evalution of whether node have similar children

=cut

sub SetAtomLevelsUsingNodeLevels {
	
	my $Self = shift;
	
	my $Nodes = $Self->getNodes;
	
	$_->getAtom->setLevel($_->getLevel) foreach (@$Nodes);
	
}

=head2 SetupNodeConnectivity

Usage: $LookUpTableInstance->SetupNodeConnectivity;

Arguments:

Synopsis:
  Builds ChildNodeGrids, which allows evalution of whether node have similar children

=cut

sub SetupNodeConnectivity {
	
	my $Self = shift;
	
	my $Nodes = $Self->getNodes;
	
	#Sort Node objects by level from highest to lowest
  my @SortedNodes = sort {$b->getLevel <=> $a->getLevel } @$Nodes;

  foreach my $Node (@SortedNodes) {
		
    my $Children = $Node->getChildren;

    my $NodeChildScoreGrid = { };

		#This portion of code is to create child comparison hashes, this is completely for computational time saving, what these hashes are string representations
	  # of the path to each of a nodes children. "C.3-C.4-N.3:02", would be a node at level 0 with state C.3 which has a child of C.4 which has a child of N.3 
	  # 02 on the end is a quick way to compare what level this string is this is in place so when comparing these hashes to each other we only
	  # compare elements at a certain level. A finished hash would look something like this, this would be for a ethane molecule starting at one of the carbons.
	  # key 					count
	  # C.4:00				1
	  # C.4-C.4:01		1
		# C.4-H.1:01		3
		# C.4-C.4-H:02  3
		
		#Known Issue is that some atoms get double counted in rings, See MoleculeFileHandler/TestPDBFiles/1byj_ligand_correct.pdb, C17's SizeOfChildScoreGrid is 34 
		#instead of 27. I am not sure if this is a problem or not

    foreach my $Child (@$Children) {
	
	    my $ChildNodeChildScoreGrid = $Child->getChildScoreGrid;
	
	    while(my ($State, $ChildGridElement) = each %$ChildNodeChildScoreGrid ) {
				
			  if(exists $NodeChildScoreGrid->{$State}) { $NodeChildScoreGrid->{$State}->Merge($ChildGridElement) }
				else																	   { $NodeChildScoreGrid->{$State} = $ChildGridElement->Copy }
				
		  }
		
    }

    $NodeChildScoreGrid->{$Node->getString . "-" . $_} = delete $NodeChildScoreGrid->{$_} foreach( keys %$NodeChildScoreGrid);
		
		if($Node->getString !~ /\|/) {

	    my $GridElement = LookUpTable::GridElement->New;
	    $GridElement->AddNodeReference(\$Node);
			
			$NodeChildScoreGrid->{$Node->getString. ":" . sprintf("%02s", $Node->getLevel)} = $GridElement;
			
		}
		
		#Handle compound elements such as F|CL|BR|I.1
		else {
			
			my $NumOfBonds = substr($Node->getString, -1);
		  
		  my @Elements;
		
		  if($NumOfBonds !~ /\d+/) {
			
			  @Elements = split /\|/, $Node->getString;
			
			  $NumOfBonds = undef;
			
		  }
		
		  else {
			
		    @Elements = split /\|/, substr($Node->getString, 0, -2);
		
		  }
							
      foreach my $Element (@Elements) {
	
	      my $GridElement = LookUpTable::GridElement->New;
		    $GridElement->AddNodeReference(\$Node);
		
		    if(defined $NumOfBonds) {
		
		      $NodeChildScoreGrid->{"$Element.$NumOfBonds:" . sprintf("%02s", $Node->getLevel)} = $GridElement;
		
	      }
	
	      else {
		
		      $NodeChildScoreGrid->{"$Element:" . sprintf("%02s", $Node->getLevel)} = $GridElement;
		
	      }
	
			}
			
		}
				
    my $Total = 0;
	  $Total += $_->getCount foreach (values %$NodeChildScoreGrid);

	  $Node->setChildScoreGrid($NodeChildScoreGrid);
	  $Node->setSizeOfChildScoreGrid($Total);
 
  } 

}

=head2 Stringify 

Usage: $LookUpTableInstance->Stringify;

Arguments:

Synopsis:
  prints out a MATCH string for LookUpTable

=cut

sub Stringify {
	
	my $Self = shift;
	
	my $Nodes = $Self->getNodes;
	
	$Self->Reset;
	
	my $LowestLevel = 0;
	
	foreach (@$Nodes) { $LowestLevel = $_->getLevel if $LowestLevel < $_->getLevel }
	
	my $FinalString;
	my %NodeToStringLookUp;
	
	while ( $LowestLevel >= 0) {
		
		foreach my $Node (@$Nodes) {
			
			next if $Node->getLevel != $LowestLevel;
			
			my $ChildNodes = $Node->getChildren;
			my @ChildStrings;
		
		  foreach (@$ChildNodes) {
			  
			  next if $_->getUsed;
			
				push @ChildStrings, defined $NodeToStringLookUp{$_} ? $NodeToStringLookUp{$_} : $_->getAtom->Stringify;
				
				$_->setUsed($_);
			
		  }
			
			my @SortedChildStrings = sort { $b cmp $a } @ChildStrings;
			@ChildStrings = map { '(' . $_ . ')' } @SortedChildStrings;
			
			$NodeToStringLookUp{$Node} = $Node->getAtom->Stringify . join("", @ChildStrings); 
			
			$FinalString = $NodeToStringLookUp{$Node} if $Node->getLevel == 0;
			
		}
		
		$LowestLevel--;
				
	}
	
	$Self->Reset;
	
	return $FinalString;
	
}


sub Update {
	my $Self = shift;
	
	my @Nodes = @{$Self->getNodes};
	
	foreach my $Node (@Nodes) { $Node->Update }
}

1;
__END__
