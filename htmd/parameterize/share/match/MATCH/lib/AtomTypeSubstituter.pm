package AtomTypeSubstituter;

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
use MATCHBaseObject;
use BaseObject ':vars';
use Type;

require Exporter;

our @ISA = qw(MATCHBaseObject Exporter);

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


sub New {
	
	my $Class = shift;
	
	my $Self = {
    _TypeStringAssociationHash     => { },	  
    _DiffNumBondsHash							 => { },	
    _IncrementTrainingData				 => [ ] 

	};
	
	bless $Self, $Class;	
	return $Self;
	
}

sub Initiate { 
	
  my $Self = shift;

  my $RelationMatrixFilePath = $Parameters->getRelationMatrixFilePath;
  my $IncrementTrainingFilePath = $Parameters->getIncrementTrainingFilePath;

  return unless defined $RelationMatrixFilePath;

  open(FILE, $RelationMatrixFilePath);

	my @TypeRelationMatrix = map { my @spl = split(/\s+/, $_); [@spl] } <FILE>;

	close(FILE);
	
	$Self->InitiateIncrementChargeMatrix;
	
	my $TypeStringAssociationHash = { };
	my $DiffNumBondsHash = { };

	foreach my $line (@TypeRelationMatrix) {
		
		my $type = shift(@$line);
		
		for (my $i = 0; $i < scalar(@$line); $i+=2) {
			
			$TypeStringAssociationHash->{$type}{$line->[$i]} = $line->[$i+1];
			
		}
	}
		
	$Self->setTypeStringAssociationHash($TypeStringAssociationHash);
  	
}

sub OldFindIncrementSubstitutionForBond {
	my ($Self, $Bond, @Increments) = @_;
	
	my @PossibleIncrementSubstituions;
	
	foreach my $Increment (@Increments) {
		
		push @PossibleIncrementSubstituions, $Increment if rindex($Increment->[0], $Bond->getPrimaryElement, 0) == 0 && rindex($Increment->[1], $Bond->getSecondaryElement) == 0;
		
		push @PossibleIncrementSubstituions, [$Increment->[1], $Increment->[0], $Increment->[3], $Increment->[2]] if rindex($Increment->[1], $Bond->getPrimaryElement, 0) == 0 && rindex($Increment->[0], $Bond->getSecondaryElement) == 0;
		
	}
	
	my ($BestMatch, $BestMatchingIncrement) = (0,0);
	
	foreach my $Increment (@PossibleIncrementSubstituions) {
		
		if($Bond->getPrimaryType eq $Bond->getSecondaryType && $Increment->[0] ne $Increment->[1] ) { next;}

    if(!$Self->NewBondSharesMissingChargeDifferenceAsIncrement($Bond, $Increment)) { next; }

    my $FowardSubstitutionScore = $Self->ScoreIncrementSubstitutionForBond($Bond, $Increment);
    $Bond->SwapPrimaryAtom;

    my $ReversedSubstitutionScore  = $Self->ScoreIncrementSubstitutionForBond($Bond, $Increment);
    $Bond->SwapPrimaryAtom; 

    if($BestMatch < $FowardSubstitutionScore)   { $BestMatch = $FowardSubstitutionScore; $BestMatchingIncrement = $Increment}
    if($BestMatch < $ReversedSubstitutionScore) { $BestMatch = $ReversedSubstitutionScore; $BestMatchingIncrement = [$Increment->[1], $Increment->[0], $Increment->[3], $Increment->[2]]}
		
	}
	
	return $BestMatchingIncrement;
	
}

sub FindIncrementSubstitutionForBond {
	my ($Self, $Bond, @Increments) = @_;
	
	my $OrginalIncrement =[ map { $_->getType } $Bond->ToAtomArray ];
	
	my @PossibleIncrementSubstituions;
	
	foreach my $Increment (@Increments) {
		
		#next if $OrginalIncrement eq $Increment;
				
		push @PossibleIncrementSubstituions, $Increment if substr($Increment->[0], 0, 1) eq substr($OrginalIncrement->[0], 0, 1) && substr($Increment->[1], 0, 1) eq substr($OrginalIncrement->[1], 0, 1);
		
		push @PossibleIncrementSubstituions, [$Increment->[1], $Increment->[0], $Increment->[3], $Increment->[2]] if substr($Increment->[0], 0, 1) eq substr($OrginalIncrement->[1], 0, 1) && substr($Increment->[0], 0, 1) eq substr($OrginalIncrement->[1], 0, 1)
	
	}
		
	my ($BestMatch, $BestMatchingIncrement) = (0,0);
	
	foreach my $Increment (@PossibleIncrementSubstituions) {
		
		next if $OrginalIncrement->[0] eq $OrginalIncrement->[1] && $Increment->[0] ne $Increment->[1];
		
    next unless $Self->NewBondSharesMissingChargeDifferenceAsIncrement($Bond, $Increment);

    my $FowardSubstitutionScore = $Self->NewScoreTypeRelation($OrginalIncrement->[0], $Increment->[0]) + $Self->NewScoreTypeRelation($OrginalIncrement->[1], $Increment->[1]);

    my $ReversedSubstitutionScore  = $Self->NewScoreTypeRelation($OrginalIncrement->[1], $Increment->[0]) + $Self->NewScoreTypeRelation($OrginalIncrement->[0], $Increment->[1]);
		
    if($BestMatch < $FowardSubstitutionScore)   { $BestMatch = $FowardSubstitutionScore; $BestMatchingIncrement = $Increment}
    if($BestMatch < $ReversedSubstitutionScore) { $BestMatch = $ReversedSubstitutionScore; $BestMatchingIncrement = [$Increment->[1], $Increment->[0], $Increment->[3], $Increment->[2]]}
				
	}
	
  return $BestMatchingIncrement;
	
}


sub ScoreIncrementSubstitutionForBond {
	
	my ($Self, $Bond, $Increment, $SkipDiffNumBondCheck) = @_;
	
	my $Score = $Self->ScoreTypeRelation($Bond->getPrimaryType, $Increment->[0]) + $Self->ScoreTypeRelation($Bond->getSecondaryType, $Increment->[1]);
	
	return $Score;
	
}


sub ScoreArrayOfTypes {
	my ($Self, $ArrayOfTypes1, $ArrayOfTypes2, $Length, $SkipDiffNumBondCheck) = @_;
	
	my $ForwardScore = 0;
	my $ReverseScore = 0;
	
	for(my $i = 0; $i < $Length; $i++) {
		$ForwardScore += $Self->ScoreTypeRelation($$ArrayOfTypes1[$i], $$ArrayOfTypes2[$i], $SkipDiffNumBondCheck);
	}
	
	my $j = $Length - 1;
	
	for(my $i = 0; $i < $Length; $i++) {
		$ReverseScore += $Self->ScoreTypeRelation($$ArrayOfTypes1[$i], $$ArrayOfTypes2[$j], $SkipDiffNumBondCheck);
		$j--;
	}
	
	return $ForwardScore > $ReverseScore ? $ForwardScore : $ReverseScore;
}


sub NewScoreArrayOfTypes {
	
	my ($Self, $Array1, $Array2) = @_;
	
	my ($ForwardScore, $ReverseScore) = (0,0);
	
	foreach (0 .. @$Array1-1) {
		
		$ForwardScore += $Self->NewScoreTypeRelation($Array1->[$_], $Array2->[$_]);
		
	}
	
  my $j = @$Array1 - 1;
	
	foreach (0 .. @$Array1-1) {
		
		$ReverseScore += $Self->NewScoreTypeRelation($Array1->[$_], $Array2->[$j]);
		
		$j--;
		
	}
		
	return $ForwardScore > $ReverseScore ? $ForwardScore : $ReverseScore;
	
}

sub NewScoreTypeRelation {
	
	my ($Self, $Type1, $Type2) = @_;

  return 1 if $Type2 eq "X";
	
	my $TypeStringAssociationHash = $Self->getTypeStringAssociationHash;
		
	return exists $TypeStringAssociationHash->{$Type1}{$Type2} ? $TypeStringAssociationHash->{$Type1}{$Type2} : 0;
  
}


sub ScoreTypeRelation {
	my ($Self, $FirstType, $SecondType, $SkipDiffNumBondCheck) = @_;
	
	if(defined $SkipDiffNumBondCheck) {
	
	  my $DiffNumBondsHash = $Self->getDiffNumBondsHash;
	
  	#if(! exists $DiffNumBondsHash->{$FirstType}{$SecondType}) { return 0 }
	
	  #if($DiffNumBondsHash->{$FirstType}{$SecondType} == 0) { return 0} 
	
  }
	
  my $TypeStringAssociationHash = $Self->getTypeStringAssociationHash;

  if(! exists $TypeStringAssociationHash->{$FirstType}{$SecondType}) { return 0 }
	else 																															 { return $TypeStringAssociationHash->{$FirstType}{$SecondType} }
}

sub BondSharesMissingChargeDifferenceAsIncrement {
	my ($Self, $Bond, $Increment) = @_;
	my ($FirstType, $SecondType) = ($Increment->[0], $Increment->[1]);
	my (@FirstTypeAllowedMissingCharge, @SecondTypeAllowedMissingCharge);
	
	my @IncrementTrainingData = @{$Self->getIncrementTrainingData};
	
	foreach my $TrainedOnIncrement (@IncrementTrainingData) {
		if(($FirstType eq $TrainedOnIncrement->[0] && $SecondType eq $TrainedOnIncrement->[1]) || 
		   ($FirstType eq $TrainedOnIncrement->[1] && $SecondType eq $TrainedOnIncrement->[0])) { next; }
		
		if($FirstType eq $TrainedOnIncrement->[0])  { push(@FirstTypeAllowedMissingCharge, $TrainedOnIncrement->[4]) }
		if($FirstType eq $TrainedOnIncrement->[1])  { push(@FirstTypeAllowedMissingCharge, $TrainedOnIncrement->[5]) }
		if($SecondType eq $TrainedOnIncrement->[0]) { push(@SecondTypeAllowedMissingCharge, $TrainedOnIncrement->[4]) }
		if($SecondType eq $TrainedOnIncrement->[1]) { push(@SecondTypeAllowedMissingCharge, $TrainedOnIncrement->[5]) }	
	}
	
	my %SeenMissingChargeState;
	
	my @FirstTypeUniqueMissingChargeStates = grep { ! $SeenMissingChargeState{$_}++ } @FirstTypeAllowedMissingCharge;
	%SeenMissingChargeState = ( );
	
	my @SecondTypeUniqueMissingChargeStates = grep { ! $SeenMissingChargeState{$_}++ } @SecondTypeAllowedMissingCharge;
	
	my $MatchedFound = 0;
	
	foreach my $UniqueMissingChargeState ( @FirstTypeUniqueMissingChargeStates) { 
	  if($Bond->getPrimaryMissingCharge == $UniqueMissingChargeState) { $MatchedFound = 1; last }	
	}
	
	if(!$MatchedFound) { return 0 }
	$MatchedFound = 0;
	
	foreach my $UniqueMissingChargeState ( @SecondTypeUniqueMissingChargeStates) { 
	  if($Bond->getSecondaryMissingCharge == $UniqueMissingChargeState) { $MatchedFound = 1; last }	
	}
	
	if(!$MatchedFound) { return 0 }
	
	return 1;
	
}

sub InitiateIncrementChargeMatrix {
	
	my ($Self, $Increments) = @_;
	
	my $IncrementFilePath = $Parameters->getIncrementFilePath;
	
	open(FILE, $IncrementFilePath);
	
	my @Increments = map { my @spl = split /\s+/, $_; [@spl] } grep { $_ !~ /^#/ } <FILE>;
	
	close(FILE);
	
	my $ChargeRelationshipLookUp;
	
	foreach my $Increment (@Increments) {
		
		my $Foward = $Increment->[2] + $Increment->[3];
		my $Backward = $Increment->[3] + $Increment->[2];
		
		push @{$ChargeRelationshipLookUp->{$Increment->[0]}}, $Foward;
		push @{$ChargeRelationshipLookUp->{$Increment->[1]}}, $Backward;
		
	}
	
	while (my ($Type, $ChargeRelationships) = each %$ChargeRelationshipLookUp) {
		
		my @UniqueRelationships = keys % { { map { $_ => 1 } @$ChargeRelationships } };
		
		$ChargeRelationshipLookUp->{$Type} = \@UniqueRelationships;
				
	}
	
	$Self->{_ChargeRelationshipLookUp} = $ChargeRelationshipLookUp;	
}


sub NewBondSharesMissingChargeDifferenceAsIncrement {
	
	my ($Self, $Bond, $Increment) = @_;
			
	#my $ChargeRelationshipLookUp = $Self->getChargeRelationshipLookUp;
		
	#my @ChargeRelationshipsForIncrement = map { $ChargeRelationshipLookUp->{ $_ } } ($Increment->[0], $Increment->[1]);
	
	my @MissingChargeForBonds = map { $_->getMissingCharge } $Bond->ToAtomArray;
	
	my $Foward = $Increment->[2] + $Increment->[3];
	my $Backward = $Increment->[3] + $Increment->[2];
	
	if(($Foward == $Bond->getPrimaryMissingCharge && $Backward == $Bond->getSecondaryMissingCharge) || 
	   ($Backward == $Bond->getPrimaryMissingCharge && $Foward == $Bond->getSecondaryMissingCharge))  {
		
		return 1;
		
	}
	
	return 0;
				
}


