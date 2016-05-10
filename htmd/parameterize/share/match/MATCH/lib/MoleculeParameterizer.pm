package MoleculeParameterizer;

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
use MATCHFunctions ':all';
use BaseObject ':vars';

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


sub New {
	
	my $Class = shift;

  my $Self = {
	  _Parameters          => undef,
	  _AtomTypeSubstituter => undef,
	  _BondParameters      => [ ],
	  _AngleParameters     => [ ],
	  _DihedralParameters  => [ ],
	  _ImproperParameters  => [ ],
	  _NonbondedParameters => [ ]
  };

  bless $Self, $Class; 

  return $Self;	
	
	
}

sub Initiate {
	
	my ($Self, $AtomTypeSubstituter) = @_;
	
	my $ParameterFilePath = $Parameters->getParameterFilePath;
	
	my $ParameterLists = LoadAtomicParametersFromFile($ParameterFilePath);
	
	$Self->setAtomTypeSubstituter($AtomTypeSubstituter);
	$Self->setParameters($Parameters);
	$Self->setBondParameters($ParameterLists->{'Bond'});
	$Self->setAngleParameters($ParameterLists->{'Angle'});
	$Self->setDihedralParameters($ParameterLists->{'Dihedral'});
	$Self->setImproperParameters($ParameterLists->{'Improper'});
	$Self->setNonbondedParameters($ParameterLists->{'Nonbond'});
	
}

sub GenerateAllPossibleParametersForMolecule {
	
  my ($Self, $Molecule) = @_;
	
	my $MoleculeBonds = $Molecule->getBonds;
	
	my @PossibleBondParameters = map { [ $_->ToAtomArray ] }  @$MoleculeBonds;
	
	my %SeenNonbonds; 
	my @PossibleNonBondedParameters = grep { ! $SeenNonbonds{$_} ++ } map { ($_->[0], $_->[1]) } @PossibleBondParameters; 

  my @PossibleAngleParameters;
  my %SeenAngles;

  foreach my $Bond (@PossibleBondParameters) {
	
	  my @PossibleAngles = LengthenAtomArray(@$Bond);
	
		foreach my $Angle (@PossibleAngles) {
			
		  if(! exists $SeenAngles{join(" ", @$Angle)} || ! exists $SeenAngles{join(" ", reverse( @$Angle))} ) {
			 
			  $SeenAngles{join(" ", @$Angle)} = 1; push @PossibleAngleParameters, $Angle;
			
			} 	
		
		}
	
	}
		
	my @PossibleDihedralParameters;
	my %SeenDihedrals;
	
	foreach my $Angle (@PossibleAngleParameters) {
		
	  my @PossibleDihedrals = LengthenAtomArray(@$Angle);
	
		foreach my $Dihedral (@PossibleDihedrals) {
			
		  if(! exists $SeenDihedrals{join(" ", @{$Dihedral})} || ! exists $SeenDihedrals{join(" ", reverse( @$Dihedral))} ) {
			
			  $SeenDihedrals{join(" ", @$Dihedral)} = 1; push @PossibleDihedralParameters, $Dihedral;
				
			} 
				
		}
		
	}
	
	return (\@PossibleNonBondedParameters, \@PossibleBondParameters, \@PossibleAngleParameters, \@PossibleDihedralParameters);
	
}

sub LengthenAtomArray {
	
	my @AtomArray = @_;
	
	my @LengthedAtomArrays;
	my $Bonded = $AtomArray[0]->getBondedAtoms;
	
	foreach my $BondedAtom (@$Bonded) {
		
    push @LengthedAtomArrays, [$BondedAtom, @AtomArray] if! (grep { $_ eq $BondedAtom } @AtomArray);
	
	}
	
	$Bonded = $AtomArray[-1]->getBondedAtoms;

  foreach my $BondedAtom (@$Bonded) {
	
		push @LengthedAtomArrays, [@AtomArray, $BondedAtom] if ! (grep { $_ eq $BondedAtom } @AtomArray);
	
	}
	
	return @LengthedAtomArrays;
}

sub NewParameterizeMolecule {
	
	my ($Self, $Molecule) = @_;
	
	my $MoleculeFileName = defined $Molecule->getFileName ? $Molecule->getFileName : "UNK";
	
	my $AtomTypeSubstituter = $Self->getAtomTypeSubstituter; 
		
  my @AllPossibleParameters = $Self->GenerateAllPossibleParametersForMolecule($Molecule);

  my $RequiredNonBondedParameters = ExtractUniqueParameters( [map { [$_] }  @{$AllPossibleParameters[0]}] );  
  my $RequiredBondParameters      = ExtractUniqueParameters($AllPossibleParameters[1]);
  my $RequiredAngleParameters     = ExtractUniqueParameters($AllPossibleParameters[2]); 
  my $RequiredDihedralParameters  = ExtractUniqueParameters($AllPossibleParameters[3]);
  
  my $RequiredImproperParameters  = ExtractUniqueParameters($Molecule->getImpropers);  

  my (@FinalNonBondedParameters, @FinalBondParameters, @FinalAngleParameters, @FinalDihedralParameters, @FinalImproperParameters);

  my @KnownNonbondParameters   = map { $_->{'Types'} = [ $_->{'Type'} ]; $_ } @{$Self->getNonbondedParameters};
  my $KnownBondParameters      = $Self->getBondParameters;
  my $KnownAngleParameters     = $Self->getAngleParameters;
  my $KnownDihedralParameters  = $Self->getDihedralParameters;
  my $KnownImproperParameters  = $Self->getImproperParameters;

  #Handle Impropers
  foreach my $RequiredImproperParameter (@$RequiredImproperParameters) {
	
	  my $BestFitParameter = FindBestFitParameter($AtomTypeSubstituter, $RequiredImproperParameter, $KnownImproperParameters);
	
    my $BestScore = $AtomTypeSubstituter->NewScoreArrayOfTypes($RequiredImproperParameter, $BestFitParameter->{'Types'});
	
	  my @BestParameters = ($BestFitParameter);
	
    foreach my $KnownImproperParameter (@$KnownImproperParameters) {
	
	    my $Score = $AtomTypeSubstituter->NewScoreArrayOfTypes($RequiredImproperParameter, $KnownImproperParameter->{'Types'});
	
	    next if $Score != $BestScore;
	
	    push @BestParameters, $KnownImproperParameter;
	
    }

    my %NumOfXTypes;

    foreach my $BestParameter (@BestParameters) {
	
	    my $NumOfXs = 0;
	
	    foreach my $Type (@{$BestParameter->{'OriginalTypes'}}) {
		
		    $NumOfXs++ if $Type eq "X";
		
	    }
	
	    $NumOfXTypes{$BestParameter} = $NumOfXs;
	
    }

    my @SortedParameters = sort { $NumOfXTypes{$a} <=> $NumOfXTypes{$b} } @BestParameters;

    push @FinalImproperParameters, $SortedParameters[0];
	
  }

  my $BestFitNonBondedParameters = FindParameters($AtomTypeSubstituter, $RequiredNonBondedParameters, \@KnownNonbondParameters);
  
  my $BestFitBondParameters      = FindParameters($AtomTypeSubstituter, $RequiredBondParameters, $KnownBondParameters);  
  my $BestFitAngleParameters     = FindParameters($AtomTypeSubstituter, $RequiredAngleParameters, $KnownAngleParameters);
  my $BestFitDihedralParameters  = FindParameters($AtomTypeSubstituter, $RequiredDihedralParameters, $KnownDihedralParameters);
  
  #Handle multiple dihedrals with different n values
  foreach my $BestFitDihedralParameter (@$BestFitDihedralParameters) {
	
	  my @SameParameter;
	  my $StringedParameter = join " ", @{$BestFitDihedralParameter->{'OriginalTypes'}};
	
	  foreach my $KnownDihedralParameter (@$KnownDihedralParameters) {
		
		  my $StringedKnownParameter = join " ", @{$KnownDihedralParameter->{'Types'}};
		
      if($StringedParameter eq $StringedKnownParameter && $BestFitDihedralParameter->{'n'} != $KnownDihedralParameter->{'n'}) {
	
	      my $CopiedParameter = CopyParameter($KnownDihedralParameter);
	
	      $CopiedParameter->{'Types'} = $BestFitDihedralParameter->{'Types'};
	
	      push @FinalDihedralParameters, $CopiedParameter;
	
      }
		
	  }
	
	  push @FinalDihedralParameters, $BestFitDihedralParameter;
	
  }

  push @FinalBondParameters,      @$BestFitBondParameters;
  push @FinalAngleParameters,     @$BestFitAngleParameters;
  push @FinalNonBondedParameters, @$BestFitNonBondedParameters; 
  
  #Handle Appending to an existing parameter file that might have duplicate parameters
  my @FinalArray = (\@FinalBondParameters, \@FinalAngleParameters, \@FinalDihedralParameters, \@FinalImproperParameters, \@FinalNonBondedParameters);

  if($Parameters->getStandAloneParameterFile == 0) {
  	 
	  #my @KnownParameters = ($KnownBondParameters, $KnownAngleParameters, $KnownDihedralParameters, $KnownImproperParameters, \@KnownNonbondParameters);
	
	  #my $NewUniqueArray = RemoveDuplicateParameters(\@FinalArray, \@KnownParameters);
	  
	  #@FinalArray = @{$NewUniqueArray};
	
  }


  if(defined $Parameters->getAppendingParameterFilePath) {
						
		my $AppendingParameters = LoadAtomicParametersFromFile($Parameters->getAppendingParameterFilePath);
		
	  my $AppendingImpropers = $AppendingParameters->{'Improper'};
	
	  $AppendingImpropers = [] unless defined $AppendingImpropers;

    my @AppendingArray = ($AppendingParameters->{'Bond'}, $AppendingParameters->{'Angle'}, $AppendingParameters->{'Dihedral'}, $AppendingImpropers, [map { $_->{'Types'} = [ $_->{'Type'} ]; $_ } @{$AppendingParameters->{'Nonbond'}}] );

	  my $NewUniqueArray = RemoveDuplicateParameters(\@FinalArray, \@AppendingArray);
	
	  foreach my $i (0 .. 4) {
						
		  push @{$AppendingArray[$i]}, @{$NewUniqueArray->[$i]};
		
	  }
		
	  $Self->WriteParameterFile($Parameters->getAppendingParameterFilePath, @AppendingArray);
		
	  return;
	
	}
	
  $Self->WriteParameterFile("$MoleculeFileName.prm", @FinalArray);

}

sub FindParameters {
	
	my ($AtomTypeSubstituter, $RequiredParameters, $KnownParameters, $ExactMatchOnly) = @_;
	
	my @BestFitParameters;

  foreach my $RequiredParameter (@$RequiredParameters) {

	  my $BestFitParameter = FindBestFitParameter($AtomTypeSubstituter, $RequiredParameter, $KnownParameters, $ExactMatchOnly);

	  croak "Could not find a Parameter for " . join(" ", @$RequiredParameter) . "\n" if ! defined $BestFitParameter && $Parameters->getExitifNotParameterized;  
	
    next if  ! defined $BestFitParameter && !$Parameters->getExitifNotParameterized;
	
	  push @BestFitParameters, $BestFitParameter;

  }

  return \@BestFitParameters;
 
	
}

sub NewDoesParameterExist {
	
  my ($NewParameter, $AppendingParameter) = @_;

  my $NewTypes = [map { ShortenType($_) } @{$NewParameter->{'Types'}}];

  my $AppendingTypes = $AppendingParameter->{'Types'};

  my @list;
  my $Count = 0;

  foreach my $i (@$NewTypes) { push @list, $Count++ }

  my $InCommon = grep { $NewTypes->[$_] eq $AppendingTypes->[$_] } (@list);
	
	return 1 if $InCommon == @$NewTypes;
	
	@$NewTypes = reverse @$NewTypes;
	
	$InCommon = grep { $NewTypes->[$_] eq $AppendingTypes->[$_] } (@list);
	
	return 1 if $InCommon == @$NewTypes;
	
	return 0;

}

sub FindBestFitParameter {
	
	my ($AtomTypeSubstituter, $NeededParameter, $KnownParameters, $ExactMatchOnly ) = @_;
	
	my $BestFitParameter; my $BestScore = 0;
	
	my $Score;
		
	foreach my $KnownParameter (@$KnownParameters) {
		
		$Score = $AtomTypeSubstituter->NewScoreArrayOfTypes($NeededParameter, $KnownParameter->{'Types'});
		 		
		if($Score > $BestScore) {
			
			next if $ExactMatchOnly && $Score != @$NeededParameter;
			
		  $BestScore = $Score; $BestFitParameter = $KnownParameter;
			
			
		}
		
	}
	
	if($BestFitParameter) {
		
    my $CopiedParameter = CopyParameter($BestFitParameter);
		
		$CopiedParameter->{'Types'} = $NeededParameter;
		
		return $CopiedParameter;
		
	}
	
	return $BestFitParameter;

}

sub CopyParameter {
	
	my $Parameter = shift;
	
  my $CopiedParameter;

  while (my ($key, $value) = each %$Parameter) {
	
	  $CopiedParameter->{$key} = $value;
	
  }	

  @{$CopiedParameter->{'OriginalTypes'}} = @{$CopiedParameter->{'Types'}};
	
	return $CopiedParameter;
	
}

sub ExtractUniqueParameters {
	
	my $Parameters = shift;
	
	my %Seen; 
	
	my @UniqueParameters;
	
	foreach my $Parameter (@$Parameters) {
		
		@$Parameter =   map { $_->getType } @$Parameter;
		
		push @UniqueParameters, $Parameter if ! $Seen{join(" ", @$Parameter)} && ! $Seen{join(" ", reverse @$Parameter)};
		
		$Seen{join(" ", @$Parameter)} = 1;
		
	}
		
	return \@UniqueParameters;
	
}


sub WriteParameterFileWithSubs {
	
	my ($Self, $FileName, $Bonds, $Angles, $Dihedrals, $Impropers, $Nonbonds) = @_;
	
	my @PrintableBondParameters      =  map { [@{$_->{'Types'}}, $_->{'Kb'}, $_->{'b0'}] } @$Bonds;
	my @PrintableAngleParameters     =  map { [@{$_->{'Types'}}, $_->{'Ktheta'}, $_->{'Theta0'}] } @$Angles;
	my @PrintableDihedralParameters  =  map { [@{$_->{'Types'}}, $_->{'Kchi'}, $_->{'n'}, $_->{'delta'}] } @$Dihedrals;
	my @PrintableImproperParameters  =  map { [@{$_->{'Types'}}, $_->{'Kchi'}, $_->{'n'}, $_->{'delta'}] } @$Impropers;
	my @PrintableNonBondParameters   =  map { [@{$_->{'Types'}}, "0.0000" , $_->{'Epsilon'}, $_->{'Rmin'} ] } @$Nonbonds;

	open(FILE, ">$FileName");
	
	print FILE "* prm file built by MATCH\n*\n\n";
	
	print FILE "BONDS\n";
	
	print FILE $Self->PrintParameterForParameterFile($PrintableBondParameters[$_]) . "! Original: " . join(" ", @{$Bonds->[$_]->{'OriginalTypes'}}) . "\n" foreach (0 .. @PrintableBondParameters-1);
	
	print FILE "\nANGLES\n";
	
	print FILE $Self->PrintParameterForParameterFile($PrintableAngleParameters[$_]) . "! Original: " . join(" ", @{$Angles->[$_]->{'OriginalTypes'}})  . "\n" foreach (0 .. @PrintableAngleParameters-1);

  print FILE "\nDIHEDRALS\n";
  
  print FILE $Self->PrintParameterForParameterFile($PrintableDihedralParameters[$_]) . "! Original: " . join(" ", @{$Dihedrals->[$_]->{'OriginalTypes'}})  . "\n" foreach (0 .. @PrintableDihedralParameters -1);
	
  print FILE "\nIMPROPER\n";

  print FILE $Self->PrintParameterForParameterFile($PrintableImproperParameters[$_]) . "! Original: " . join(" ", @{$Impropers->[$_]->{'OriginalTypes'}}) . "\n" foreach (0 .. @PrintableImproperParameters-1);
  
  print FILE "\nNONBONDED nbxmod  5 atom cdiel shift vatom vdistance vswitch -\ncutnb 14.0 ctofnb 12.0 ctonnb 10.0 eps 1.0 e14fac 1.0 wmin 1.5\n";
	
	print FILE $Self->PrintParameterForParameterFile($PrintableNonBondParameters[$_]) . "! Original: " . join(" ", @{$Nonbonds->[$_]->{'OriginalTypes'}})  . "\n" foreach (0 .. @PrintableNonBondParameters-1);
	
	print FILE "\n\n";
 
	close(FILE);
	
	
}


sub WriteParameterFile {
	
	my ($Self, $FileName, $Bonds, $Angles, $Dihedrals, $Impropers, $Nonbonds) = @_;
	
	my @PrintableBondParameters      =  map { [@{$_->{'Types'}}, $_->{'Kb'}, $_->{'b0'}] } @$Bonds;
	my @PrintableAngleParameters     =  map { [@{$_->{'Types'}}, $_->{'Ktheta'}, $_->{'Theta0'}] } @$Angles;
	my @PrintableDihedralParameters  =  map { [@{$_->{'Types'}}, $_->{'Kchi'}, $_->{'n'}, $_->{'delta'}] } @$Dihedrals;
	my @PrintableImproperParameters  =  map { [@{$_->{'Types'}}, $_->{'Kchi'}, $_->{'n'}, $_->{'delta'}] } @$Impropers;
	my @PrintableNonBondParameters   =  map { [@{$_->{'Types'}}, "0.0000" , $_->{'Epsilon'}, $_->{'Rmin'} ] } @$Nonbonds;

	open(FILE, ">$FileName");
	
	print FILE "* prm file built by MATCH\n*\n\n";
	
	print FILE "BONDS\n";
	
	print FILE $Self->PrintParameterForParameterFile($_) . "\n" foreach (@PrintableBondParameters);
	
	print FILE "\nANGLES\n";
	
	print FILE $Self->PrintParameterForParameterFile($_) . "\n" foreach (@PrintableAngleParameters);

  print FILE "\nDIHEDRALS\n";
  
  print FILE $Self->PrintParameterForParameterFile($_) . "\n" foreach (@PrintableDihedralParameters);
	
  print FILE "\nIMPROPER\n";

  print FILE $Self->PrintParameterForParameterFile($_) . "\n" foreach (@PrintableImproperParameters);
  
  print FILE "\nNONBONDED nbxmod  5 atom cdiel shift vatom vdistance vswitch -\ncutnb 14.0 ctofnb 12.0 ctonnb 10.0 eps 1.0 e14fac 1.0 wmin 1.5\n";
	
	print FILE $Self->PrintParameterForParameterFile($_) . "\n" foreach (@PrintableNonBondParameters);
	
	print FILE "\n\n";
 
	close(FILE);
	
}


sub RemoveDuplicateParameters {
	
	my ($FinalArray, $AppendingArray) = @_;
	
	my @NewUniqueArray;

	foreach my $i (0 .. 4) {
		
		$NewUniqueArray[$i] = [];
		
  	foreach my $FinalParameter (@{$FinalArray->[$i]}) {
			
			my $Flag = 0;
			
			foreach my $AppendingParameter (@{$AppendingArray->[$i]}) {
				
				if(NewDoesParameterExist($FinalParameter, $AppendingParameter)) {
										
					$Flag = 1; last;
				
				}
				
			}
			
			if(!$Flag) {
								
				push @{$NewUniqueArray[$i]}, $FinalParameter; 
								
			}
			
		}
		
	}
	
	return \@NewUniqueArray;
	
}

sub DoesParameterExist {
	
	my ($ExistingParameters, $NewParameter, $Length) = @_;
		
	my $NewTypes = [map { ShortenType($_) } @$NewParameter[0 .. $Length-1]];
	
	my $Exists = 0;
	
	foreach my $ExistingParameter (@$ExistingParameters) {
		
	  my $Types = $ExistingParameter->{'Types'};
	
	  my $InCommon = grep { $NewTypes->[$_] eq $Types->[$_] } (0 .. @$Types-1);
		
		if($InCommon == @$Types) {
			
			$Exists = 1; last
			
		}
		
		my $ReverseTypes = [reverse (@$Types)];
		
		$InCommon = grep { $NewTypes->[$_] eq $ReverseTypes->[$_] } (0 .. @$Types-1);
	  
		if($InCommon == @$Types) {
			
			$Exists = 1; last
			
		}
		
	}
	
	return $Exists;
	
}

sub PrintParameterForParameterFile {
	
	my ($Self, $Parameter) = @_; 
	
	my @ArrayedParameter = @{$Parameter};
		
	my $FinalString;
	
	foreach my $Element (@ArrayedParameter) {
		
		if($Element =~ /^!/)             { last; }
		elsif($Element =~ /^[A-Z0-9]+$/) { $FinalString .= sprintf("%-6s", ShortenType($Element)) . " " }
		elsif($Element =~ /^\d$/) 			 { $FinalString .= sprintf("%-3s", $Element) . " " }
		elsif($Element =~ /^-?\d+\.\d+/) { $FinalString .= sprintf("%-10s", $Element) . " " }
		
	}
	
	return $FinalString;
	
}



1;
__END__



