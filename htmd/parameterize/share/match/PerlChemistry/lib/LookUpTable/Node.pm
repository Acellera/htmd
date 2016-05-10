package LookUpTable::Node;

=head1 NAME


=head1 SYNOPSIS


=head1 DESCRIPTION


=head2 EXPORT


=head1 AUTHOR

Joseph Yesselman, E<lt>jyesselm@umich.eduE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2010 by Joseph Yesselman

This library is free software; you can redistribute it and/or modifycd
it under the same terms as Perl itSelf, either Perl version 5.10.0 or,
at your option, any later version of Perl 5 you may have available.

=head1 FUNCTIONS

=cut

#use 5.010000;
use strict;
use warnings;
use Carp;

require Exporter;

our @ISA = qw(Exporter);

# Items to export into callers namespace by default. Note: do not export
# names by default without a very good reason. Use EXPORT_OK instead.
# Do not simply export all your public functions/methods/constants.

# If you do not need this, moving things directly into @EXPORT or @EXPORT_OK
# will save memory.
our %EXPORT_TAGS = ( 'all'  => [ qw( )],
											 
										 'vars' => [ qw () ],
										
										 'func' => [ qw ()]);
							

our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

our @EXPORT = qw(
	
);

our $VERSION = '0.01';

#Exported Variables


#Non Exported Variables
our $AUTOLOAD;


# Preloaded methods go here.

=head2 New

Usage: LookUpTable::Node->New($Atom);

Arguments:
  $Class: should be 'LookUpTable::Node'
  $Atom: the Atom object this node will represent

Synopsis:
  creates a new Node object

=cut

sub New {
	
	my ($Class, $Atom) = @_;
	
	croak("Usage For Node->New(\$Atom), \$Atom must be an Atom Object") if ! defined $Atom;
	croak("Usage For Node->New(\$Atom), \$Atom must be an Atom Object") if ref($Atom) ne 'Atom';
		
  my $Self = { 
	  _Atom                 => $Atom,
	  _BestScore						=> 0,
	  _Children             => [],
	  _DoubleBondToParent   => 0,
	  _HighScore            => undef, 
	  _Level                => $Atom->getLevel,
	  _NoRing				        => $Atom->getNoRing,
		_Parents              => [],
		_String               => $Atom->getState,
		_Related			        => [],
		_Rings					      => $Atom->getRings,
		_RingAtom							=> $Atom->getRingAtom,
		_MissingCharge        => $Atom->getMissingCharge,
		_Used                 => 0,
		_ChildScoreGrid       => undef,
		_SizeOfChildScoreGrid => 0
	  
	};
	
	bless $Self, $Class;
	
	$Self->setBestScore($Self->StringCompare($Self));
	
  return $Self;
}

=head2 AddParent

Usage: $NodeInstance->AddParent;

Arguments:
  $ParentNode: a Parent of $NodeInstance

Synopsis:
  Adds $ParentNode to the internal Parent array for $NodeInstance

=cut

sub AddParent {
	my ($Self, $ParentNode) = @_;
	
	confess("Cannot add this parent node to \$NodeInstance because its not a Node Object") if ref($ParentNode) ne 'LookUpTable::Node';
		
	my $Parents = $Self->getParents;
	
	push @$Parents, $ParentNode;

}

=head2 AddChild

Usage: $NodeInstance->AddChild;

Arguments:
  $ChildNode: A Child of $NodeInstance

Synopsis:
  Adds $ChildNode to the internal Parent array for $NodeInstance

=cut

sub AddChild {
	my ($Self, $ChildNode) = @_;
	
	confess("Cannot add this child node to \$NodeInstance because its not a Node Object") if ref($ChildNode) ne 'LookUpTable::Node';
		
	my $Children = $Self->getChildren;
	
	push @$Children, $ChildNode;

}


sub AddRelated {
	my ($Self, $RelatedNode) = @_;
	
	confess("Cannot add this related node to \$NodeInstance because its not a Node Object") if ref($RelatedNode) ne 'LookUpTable::Node';
		
	my $Related = $Self->getRelated;
	
	push @$Related, $RelatedNode;

}
=head2 RemoveChild

Usage: $NodeInstance->RemoveChild($ChildNode);

Arguments:
  $ChildNode: a Child of $NodeInstance

Synopsis:
  Removes $ChildNode to the internal Child array for $NodeInstance

=cut

sub RemoveChild {
	
	my ($Self, $ChildNode) = @_;
	
	my $Children = $Self->getChildren;
	
	foreach (0 .. @$Children-1) {
		
		if($Children->[$_] eq $ChildNode) { splice(@$Children, $_, 1); last; }
		
	}

}

=head2 RemoveParent

Usage: $NodeInstance->RemoveParent($ParentNode);

Arguments:
  $ParentNode: a Parent of $NodeInstance

Synopsis:
  Removes $ParentNode to the internal Parent array for $NodeInstance

=cut

sub RemoveParent {
	
	my ($Self, $ParentNode) = @_;
	
	my $Parents = $Self->getParents;
	
	foreach (0 .. @$Parents-1) {
		
		if($Parents->[$_] eq $ParentNode) { splice(@$Parents, $_, 1); last; }
		
	}

}

sub RemoveRelated {
	
	my ($Self, $RelatedNode) = @_;
	
	my $Related = $Self->getRelated;
	
	foreach (0 .. @$Related-1) {
		
		if($Related->[$_] eq $RelatedNode) { splice(@$Related, $_, 1); last; }
		
	}

}


=head2 StringCompare

Usage: $NodeInstance->StringCompare($OtherNodeInstance);

Arguments:
  $OtherNodeInstance: Another Node object which to compare the String variable to

Synopsis:
  Scores how related the String variables are, these are infact the Atom State variables, I have no idea why I am calling them String

=cut

{
	our %StringCompareCache;

  sub StringCompare {
	
	  my ($Self, $CompareTo) = @_;
	  
	  my $SelfString = $Self->getString;
	  my $CompareToString = $CompareTo->getString;
	
	  return $StringCompareCache{$SelfString . " " . $CompareToString } if exists $StringCompareCache{$SelfString . " " . $CompareToString};
	
  	return 1 if $SelfString eq "*" || $CompareToString eq "*";

    my $First;
    my $Second;

    if(length $SelfString <= length $CompareToString){
	
	    $First = $SelfString; $Second = $CompareToString;
	
    }

    else {
	
	    $First = $CompareToString; $Second = $SelfString;
	    
	
    }
 

    my @SelfCharacterArray = split /\./, $First;
    my @CompareToCharacterArray =  split /\./, $Second;

    my $count = 0;
    my $min = $#SelfCharacterArray > $#CompareToCharacterArray ? $#CompareToCharacterArray : $#SelfCharacterArray;

    foreach (0 .. $min) {
		  if($SelfCharacterArray[$_] =~ $CompareToCharacterArray[$_]) { $count++ }
		  else                                                        { last;    }
    }

    $StringCompareCache{$SelfString . " " . $CompareToString} = $count;

    return $count;
  }
}

=head2 StringCompare

Usage: $NodeInstance->StringCompare($OtherNodeInstance);

Arguments:
  $OtherNodeInstance: Another Node object which to compare the String variable to

Synopsis:
  Scores how related the String variables are, these are infact the Atom State variables, I have no idea why I am calling them String

=cut

{
	my %FullCompareCache;
 	
	sub FullStateCompare {
		
		my ($Self, $OtherNode) = @_;
		
		my ($SelfStringify, $OtherNodeStringify) = ($Self->Stringify, $OtherNode->Stringify);
		
		#if($Self->getDoubleBond)
						
    return $FullCompareCache{$SelfStringify. " " . $OtherNodeStringify} if exists $FullCompareCache{$SelfStringify . " " . $OtherNodeStringify};
				
		my $SelfRings = $Self->getRings;
		my $OtherNodeRings = $OtherNode->getRings;
		
		#print $Self->Stringify . " " . $OtherNode->Stringify . " " . $Self->getMissingCharge . " " . $OtherNode->getMissingCharge . "\n";
		
		if($Self->getDoubleBondToParent == 1 && $OtherNode->getDoubleBondToParent == 0) {
						
		  $FullCompareCache{$SelfStringify . " " . $OtherNodeStringify} = 0; return 0;
			
		}
				
		if($Self->getNoRing == 1 && $OtherNode->getRingAtom == 1 ||  $OtherNode->getNoRing == 1 && $Self->getRingAtom == 1 || 
			$OtherNode->getRingAtom == 0 && $Self->getRingAtom == 1 )  { 
			
			$FullCompareCache{$SelfStringify . " " . $OtherNodeStringify} = 0; return 0;
			
		}		
			
		if($Self->getMissingCharge != 0 && 
			($Self->getMissingCharge > 0 && $OtherNode->getMissingCharge <= 0 || $Self->getMissingCharge < 0 && $OtherNode->getMissingCharge >= 0)) {
			
			$FullCompareCache{$SelfStringify . " " . $OtherNodeStringify} = 0; return 0;
			
		}
		
 		my %SeenOtherNodeRing;
    my %SeenSelfRings;

    my @Pairs;
    my @Matches;

    foreach my $i (0 .. @$SelfRings-1) {

	    foreach my $j (0 .. @$OtherNodeRings-1) {

 		    my $Ring = $OtherNodeRings->[$j];
								
				if($SelfRings->[$i] =~ /\*/ || $OtherNodeRings->[$j] =~ /\*/) {
					
					next if length($SelfRings->[$i]) == 2 && substr($SelfRings->[$i],-1) ne substr($OtherNodeRings->[$j],-1);
					
				  push @Pairs, [$i, $j];
					push @Matches, { $i . "i" => 1, $j . "j" => 1};
					next;
					
				}
				
				if($SelfRings->[$i] =~ qr/$Ring/ || length($SelfRings->[$i]) == 1 && $SelfRings->[$i] eq substr($Ring,0,1) ) { 
					#print $SelfRings->[$i] . " " . $OtherNodeRings->[$j] . "\n";
					
					push @Pairs, [$i, $j];
					push @Matches, { $i . "i" => 1, $j . "j" => 1};
					
				}

		  }

	  }
	
	  my %SeenMatches;
	  my $Done  = 0;
	  my $MatchedAllRings = 0;
	  my @NewMatches;
	  my $Count = 0;
	
	  while( !$Done ) {
		
	    foreach my $Match (@Matches) {
		
		    if(keys %$Match == @$SelfRings*2) {
			
			    $MatchedAllRings = 1; $Done = 1; last;
			
		    }
		
		    foreach my $Pair (@Pairs) {
			
			    next if exists $Match->{$Pair->[0] . "i"} ||  exists $Match->{$Pair->[1] . "j"};
			
			    my %NewMatch = %$Match;
			
			    $NewMatch{$Pair->[0] . "i"} = 1; $NewMatch{$Pair->[1] . "j"} = 1;
			
			    next if exists $SeenMatches{ join(" ", sort keys %NewMatch)};
			
			    $SeenMatches{ join(" ", sort keys %NewMatch)} = 1;
					
			    push @NewMatches, \%NewMatch;  
			
		    }
		
	    }
	
	    $Count++ if @Matches == @NewMatches;
	
	    @Matches = @NewMatches;	  
	
	    $Done = 1 if $Count > 2;
	
	  }
	 			
		
	#	print "\n" . @$OtherNodeRings . " " .  (keys %SeenOtherNodeRing) . " " . @$SelfRings . "\n";
		
		if(!$MatchedAllRings && @$SelfRings != 0 ) { 
			
			$FullCompareCache{$SelfStringify . " " . $OtherNodeStringify} = 0; return 0;
			
		}
		
		my $StringCompare = $Self->StringCompare($OtherNode);
		
		if ($StringCompare < $Self->getBestScore && $StringCompare < $OtherNode->getBestScore) {
			
			$FullCompareCache{$SelfStringify . " " . $OtherNodeStringify} = 0; return 0;
			
		}
		
		
		
		$FullCompareCache{$SelfStringify . " " . $OtherNodeStringify} = $StringCompare + @$SelfRings;
		
		return $FullCompareCache{$SelfStringify . " " . $OtherNodeStringify};
		
	}
}


sub FindBondTypeWithParent {
	
	my $Self = shift;
	
	my $Count = 0;
	
	foreach my $Parent (@{$Self->getParents}) {
	
	  $Count++;
	
	  next if ! defined $Parent->getUsed;
	
	  my $Atom = $Self->getAtom;  
	
	  my $ParentAtom = $Parent->getAtom;
	
	  foreach my $Bond (@{$Atom->getBonds}) {
		
		  next if $Bond->getBondingPartner($Atom) ne $ParentAtom;
	 	
	    return $Bond->getType; 
	
	  }
	
  } 

  croak ("Could not find bond type with parent node") if $Count > 0;
	
}

=head2 FullMatchesTo

Usage: $NodeInstance->FullMatchesTo($OtherNodeInstance);

Arguments:
  $OtherNodeInstance: Another Node object which to compare the String variable to

Synopsis:
  Scores how related the String variables are, these are infact the Atom State variables, I have no idea why I am calling them String

=cut

sub FullMatchesTo {
	
	my ($Self, $OtherNode) = @_;
	
	#return 0 if $Self->getLevel != $OtherNode->getLevel;
		
	my $StringScore = $Self->FullStateCompare($OtherNode);
	return 0 if $StringScore < $Self->getBestScore && $StringScore < $OtherNode->getBestScore;
	
	return 1;
		
}

=head2 StringMatchesTo

Usage: $NodeInstance->StringMatchesTo($OtherNodeInstance);

Arguments:
  $OtherNodeInstance: Another Node object which to compare the String variable to

Synopsis:
  Scores how related the String variables are, these are infact the Atom State variables, I have no idea why I am calling them String

=cut

sub StringMatchesTo {
	
	my ($Self, $OtherNode) = @_;
	
	#return 0 if $Self->getLevel != $OtherNode->getLevel;
		
	my $StringScore = $Self->StringCompare($OtherNode);
	return 0 if $StringScore < $Self->getBestScore && $StringScore < $OtherNode->getBestScore;
	
	return 1;
		
}

=head2 UpdateChildScoreGrid

Usage: $NodeInstance->UpdateChildScoreGrid;

Arguments:

Synopsis:
  Required during LookUpTable::DeterminedBestMatchUps, updates whether a node has been used or not in ChildScoreGrid

=cut

sub UpdateChildScoreGrid {

	my $Self = shift;
	
	$_->Update foreach (values %{$Self->getChildScoreGrid});
	
}

=head2 Update

Usage: $NodeInstance->Update;

Arguments:

Synopsis:
  Required during LookUpTable::DeterminedBestMatchUps, updates whether a node has been used or not in ChildScoreGrid

=cut

sub Update {
	
	my $Self = shift;
	
	my $SelfAtom = $Self->getAtom;
	
	$Self->setNoRing($SelfAtom->getNoRing);
	$Self->setRings($SelfAtom->getRings);
	$Self->setRingAtom($SelfAtom->getRingAtom);
	$Self->setString($SelfAtom->getState);
	$Self->setMissingCharge($SelfAtom->getMissingCharge);
	
	$Self->setDoubleBondToParent(1) if !$Self->getDoubleBondToParent && $Self->FindBondTypeWithParent() == 2;
	
}


sub Stringify {
	
	my $Self = shift;
			
	my $SelfState = $Self->getString;
	my $SelfNoRing = $Self->getNoRing;
	my @SelfRings = @{$Self->getRings};
				
	my $String = "";
	
	$String .= "=" if($Self->getDoubleBondToParent);
	
	$String .= "!" if($SelfNoRing);
	
	$String .= $SelfState;
	
	$String .= "%" if $Self->getRingAtom;
	
	$String .= join(",", @SelfRings) if scalar(@SelfRings) > 0;
	
	$String .= "+" if $Self->getMissingCharge > 0;
  $String .= "-" if $Self->getMissingCharge < 0;
	
	return $String;
	
}




sub AUTOLOAD {
  my ($Self, $newvalue) = @_;
  my $type = ref($Self) 
             or print "$Self is not an object\n";


  if($AUTOLOAD =~ /DESTROY/ ) { return; } 

  my ($operation, $attribute) = ($AUTOLOAD =~ /(get|set)(\w+)$/);
  
  confess("Method name $AUTOLOAD in " . ref($Self) . " is not in the recognized form (get|set)_attribute\n") unless ($operation && $attribute);
	
	my $AtomAttribute;
	if(! exists $Self->{"_" . $attribute}) {
		$AtomAttribute = $1 if $attribute =~ /Atom(\w+)/;
		confess("Parameter $AtomAttribute does not exist in Atom") unless (exists $Self->{"_Atom"}->{"_" . $AtomAttribute});
	  
	}
	else {
	  confess("Parameter $attribute does not exist in " . ref($Self)) unless (exists $Self->{"_" . $attribute});
	}
	
  no strict "refs";
 
  if($operation eq 'set')  { 
	  $Self->{"_" . $attribute} = $newvalue;

	  *{$AUTOLOAD} = sub { 
		  						   my ($Self, $newvalue) = @_;
		  							 $Self->{"_" . $attribute} = $newvalue; 	
								   };
  }

  elsif($operation eq 'get')    { 
	
	  if($AtomAttribute) {
	    *{$AUTOLOAD} = sub { 
			                 my ($Self) = @_;
					  					 $Self->{"_Atom"}->{"_" . $AtomAttribute};
					           };
					
		  use strict "refs";

			return $Self->{"_Atom"}->{"_" . $AtomAttribute};	
			
		}
	
	  *{$AUTOLOAD} = sub { 
		                 my ($Self) = @_;
									   $Self->{"_" . $attribute};
									 }; 
	}

  use strict "refs";

  #confess("Cannot retrieve parameter: " . $attribute . " because it is undefined, add its value to your .par file\n") if(! defined $Self->{"_" . $attribute});

  return $Self->{"_" . $attribute};
}




1;
__END__
