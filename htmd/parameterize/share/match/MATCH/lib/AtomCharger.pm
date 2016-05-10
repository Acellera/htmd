package AtomCharger;

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
	
    _AtomTypeSubstituter			     => undef,	  
    _Increments										 => [ ],
    _IncrementLookUpTable          => { },
    _RefiningIncrementsLookUp      => { },

	};
	
	bless $Self, $Class;	
	return $Self;
	
}


sub Initiate {
	
	my ($Self, $AtomTypeSubstituter) = @_;
	
	my $IncrementFilePath = $Parameters->getIncrementFilePath;
	my $RefiningIncrementsFilePath = $Parameters->getRefiningIncrementsFilePath;
	
	croak("IncrementFilePath is not set, cannot read .incr file for this Forcefield") unless defined $IncrementFilePath;
	
	open(FILE, $IncrementFilePath);
	
	my @Increments = map { my @spl = split /\s+/, $_; [@spl] } grep { $_ !~ /^#/ } <FILE>;
	
	close(FILE);
	
	$Self->setIncrements(\@Increments);
	
	my %IncrementLookUpTable;
	
	foreach my $Incr (@Increments) {
		
		if($Incr->[0] le $Incr->[1]) { 
			
		  $IncrementLookUpTable{$Incr->[0] . " " . $Incr->[1]} = [$Incr->[2], $Incr->[3]];
						
		}
		
		else {
		
		  $IncrementLookUpTable{$Incr->[1] . " " . $Incr->[0]} = [$Incr->[3], $Incr->[2]];				
			
		}
		
	}
		
	$Self->setIncrementLookUpTable(\%IncrementLookUpTable);
	
	$Self->setAtomTypeSubstituter($AtomTypeSubstituter);
	
	$Self->SetupRefinementMATCHStrings($RefiningIncrementsFilePath) if defined $RefiningIncrementsFilePath;
	
}

sub SetupRefinementMATCHStrings {
	
	my ($Self, $RefiningIncrementsFilePath) = @_;
	
	open(FILE, $RefiningIncrementsFilePath);
	my @RefiningIncrementFileContents =  <FILE>;
	close(FILE);
	
	my @RefiningIncrements;
	my @NameShortCuts = grep { $_ =~ /\s*\S+\s+=\s+\S+\s+/ } @RefiningIncrementFileContents;
	
	my @Refinements = grep { $_ !~ /\s*\S+\s+=\s+\S+\s+/ } @RefiningIncrementFileContents;
	
  my $Substituted = ApplyStringNamingShortCuts(\@NameShortCuts, \@Refinements);	
	
	foreach my $RefiningIncrementStringRepresentation (@$Substituted) {
			
		next if $RefiningIncrementStringRepresentation =~ /^(?:\#| )/;
		
		my @spl = split(/\s+/, $RefiningIncrementStringRepresentation);
		
		next if @spl != 6;
		
		my $RefiningType = Type->New($spl[0], $spl[4], $spl[5]);
		
		$RefiningType->Initiate;
		
		push @RefiningIncrements, [$spl[0], $spl[1], $spl[2], $spl[3], $RefiningType] ;

	}
	
	my $RefiningIncrementsLookUp = { };
	
	foreach my $RefiningIncrement (@RefiningIncrements) {
		
		if(exists $RefiningIncrementsLookUp->{join " ", sort ($RefiningIncrement->[0], $RefiningIncrement->[1])}) {
			
			push @{$RefiningIncrementsLookUp->{join " ", sort ($RefiningIncrement->[0], $RefiningIncrement->[1]) }}, [$RefiningIncrement->[2], $RefiningIncrement->[3], $RefiningIncrement->[4]];
		
		}
		
		else {
			
			$RefiningIncrementsLookUp->{join " ", sort ($RefiningIncrement->[0], $RefiningIncrement->[1]) } = [[ $RefiningIncrement->[2], $RefiningIncrement->[3], $RefiningIncrement->[4] ]]; 
		
		}
	}
	
	while( my ($Bond, $Increments) = each(%$RefiningIncrementsLookUp)) {
		
		my @SortedIncrements = sort { length $b->[2]->getString <=> length $a->[2]->getString } @$Increments;
		
		$RefiningIncrementsLookUp->{$Bond} = [@SortedIncrements]; 
	
	}
	
	$Self->setRefiningIncrementsLookUp($RefiningIncrementsLookUp);	
	
}


sub ChargeAtomsInChain {
	
	my ($Self, $Chain) = @_;
	
	my @Bonds = map { @{$_->getBonds} } @{$Chain->getMolecules};
	
	#If atoms are already charged for some reason

# MJH
#	foreach my $Bond (@Bonds) {
#		$Bond->getPrimary->setCharge(0); $Bond->getSecondary->setCharge(0);
#	}
	
	my $Result =  $Self->ChargeAtoms(\@Bonds);
	
	return 0 unless $Result == 1;
	
	foreach my $Molecule (@{$Chain->getMolecules}) {
		
		my $TotalCharge = 0;
		
		my $Atoms = $Molecule->getAtoms;
		
		foreach my $Atom (@$Atoms) { 
			
			$TotalCharge += $Atom->getCharge if $Atom->getMolecule eq $Molecule;
		
		}
		
		if(sprintf("%.0f", $TotalCharge) != $TotalCharge) { 
			
		 # return 0 if $Parameters->getExitifNotCharged;
		
		  #croak("Charge is not interger!");
			
		}

		$Molecule->setCharge($TotalCharge);
		
	}
	
	return $Result;
	
	
}


sub ChargeAtomsInMolecule {
	
	my ($Self, $Molecule) = @_;
	
	my $MoleculeBonds = $Molecule->getBonds;
	
	my $Result = $Self->ChargeAtoms($MoleculeBonds);
	
	return 0 unless $Result == 1;
	
  my $MoleculeAtoms = $Molecule->getAtoms;
	my $TotalCharge = 0;
	
	foreach my $Atom (@$MoleculeAtoms) { $TotalCharge += $Atom->getCharge }
	
	if(sprintf("%.0f", $TotalCharge) != $TotalCharge) { 
		
	 # return 0 if $Parameters->getExitifNotCharged;
	
	 # croak("Charge is not interger!");
		
	}
	
	$Molecule->setCharge($TotalCharge);
	
	return $Result;

	
}


sub ChargeBond {
	
	my ($Self, $Bond, $Print) = @_;
	
	return 1 if $Parameters->getDoNotCharge == 1;
		
	my $Increment;
		
	#Handle Refinement Increments
  $Increment = $Self->FindRefinementIncrement($Bond) if $Parameters->getUsingRefiningIncrements;

  #Found a Refinement Increment = Sucess
  return SuccessfulDeterminationofIncrement($Bond, $Increment, $Print) if defined $Increment && scalar @$Increment > 0;
		
  #Look for normal increment, as no refinement existed or no refinement was desired
	my $IncrementLookUpTable = $Self->getIncrementLookUpTable;
	
  if(exists $IncrementLookUpTable->{SortedTypeString($Bond)}) { 
		
	  $Increment = $IncrementLookUpTable->{SortedTypeString($Bond)}; 
		
		return SuccessfulDeterminationofIncrement($Bond, $Increment, $Print);
			
	}
	
	print SortedTypeString($Bond) . "\n";
	
	#Try and find substitute if we have not suceeded yet
	return -2 if !$Parameters->getSubstituteIncrements && !$Parameters->getExitifNotCharged;
	
	croak "the increment between atoms " . $Bond->getPrimaryName . " " . $Bond->getPrimary->Stringify . " (" . $Bond->getPrimaryType . ") and " . $Bond->getSecondaryName . " " . $Bond->getSecondary->Stringify . " (" . $Bond->getSecondaryType . ") does not exist 
	and Increment Substitution is turned off, consider turning it on -SubstituteIncrements 1\n"  if !$Parameters->getSubstituteIncrements;
	
	my $AtomTypeSubstituter = $Self->getAtomTypeSubstituter;
	
	croak "SubstituteIncrements is turned on yet no AtomTypeSubstituter object is defined!" unless defined $AtomTypeSubstituter;
	
	my $SubstituteIncrement = $AtomTypeSubstituter->FindIncrementSubstitutionForBond($Bond, @{$Self->getIncrements});
	
	return SuccessfulDeterminationofIncrement($Bond, $SubstituteIncrement, $Print) if defined $SubstituteIncrement;
	
	#No substitute was found
	
	#Theory is that if they are truely the same type there should be no exchange of charge upon bonding since both have the same electron 
	#cloud density, this is only imposed if no other substitue can be found
	return SuccessfulDeterminationofIncrement($Bond, [0,0], $Print) if $Bond->getPrimaryType eq $Bond->getSecondaryType;
	
	#Failure but ExitifNotCharged is flaged thus only an error is sent and does not crash out
  return -2 if !$Parameters->getExitifNotCharged;

  croak "the increment between atoms " . $Bond->getPrimaryName . " " . $Bond->getPrimary->Stringify . " (" . $Bond->getPrimaryType . ") and " . $Bond->getSecondaryName . " " . $Bond->getSecondary->Stringify . " (" . $Bond->getSecondaryType . ") does not exist 
	and Increment Substitution failed to yeild an exceptable result\n";

}

sub SuccessfulDeterminationofIncrement {
	
	my ($Bond, $Increment, $Print) = @_;
	
	$Bond->ApplyIncrement(@$Increment);
		
	print sprintf("%-7s", $Bond->getPrimaryName) . " " . sprintf("%-7s", $Bond->getPrimaryType)  . " "  . sprintf("%-7s", $Bond->getSecondaryName) . " " . sprintf("%-7s", $Bond->getSecondaryType) . " " . sprintf("%-7s", $Bond->getType)  . " " . sprintf("%-7s", $Bond->getPrimary->getChain->getName) . " " . sprintf("%+9s", $Increment->[0]) . " - " . sprintf("%+9s", $Increment->[1]) . "\n" if defined $Print && $Print == 1;
  
  #Sucess!!
  return 1;
	
}

=head2 FindRefinementIncrement

Usage: $AtomChargerInstance->FindRefinementIncrement($Bond);

Arguments:
  $Bond: A bond object 

Synopsis:
  Locates if any refinement increment matches the chemical space of either atom in the bond object

=cut

sub FindRefinementIncrement {
	
	my ($Self, $Bond) = @_;
	
	my $RefiningIncrementsLookUp = $Self->getRefiningIncrementsLookUp;
		
	my $ChargeToApplyToBond;
	
	my $PotentialRefinements = $RefiningIncrementsLookUp->{SortedTypeString($Bond)} if  $RefiningIncrementsLookUp->{SortedTypeString($Bond)};
		 
  foreach my $Refinement (@$PotentialRefinements) {
	
	  if($Refinement->[2]->getLookUpTable->AreAllNodesSharedBetween([$Bond->getPrimaryLookUpTable])) {
												
		  $ChargeToApplyToBond = [$Refinement->[0], $Refinement->[1]]; last;
			
		}
		
		elsif($Refinement->[2]->getLookUpTable->AreAllNodesSharedBetween([$Bond->getSecondaryLookUpTable])) {
			
		  $ChargeToApplyToBond = [$Refinement->[1], $Refinement->[0]]; last;
			
		}
	
  }

  return $ChargeToApplyToBond;
	
}





sub ChargeAtoms {
	
	my ($Self, $Bonds) = @_;
	
	my $ChargeComplete = 1;
	
	my @NoAvailableIncrements;
	
	return 1 if $Parameters->getDoNotCharge == 1;
	
	my $AtomTypeSubstituter = $Self->getAtomTypeSubstituter;

	my $UsingRefiningIncrements = $Parameters->getUsingRefiningIncrements;

	my $SubstituteIncrements = $Parameters->getSubstituteIncrements;
	
	my $ExitifNotCharged = $Parameters->getExitifNotCharged;
	
  croak("SubstituteIncrements is turned on yet no AtomTypeSubstituter object is defined!") if $SubstituteIncrements && ! defined $AtomTypeSubstituter;
	
	my $RefiningIncrementsLookUp = $Self->getRefiningIncrementsLookUp;
	
	my $Increments = $Self->getIncrements;
		
	my %IncrementLookUpTable = map { $_->[0] . "-" . $_->[1] => [$_->[2], $_->[3]] } @$Increments;
	
	foreach my $Bond (@$Bonds) {
				
		my $PotentialRefinements = [];
		
		if($UsingRefiningIncrements) {
						
		  $PotentialRefinements = $RefiningIncrementsLookUp->{SortedTypeString($Bond)} if exists $RefiningIncrementsLookUp->{SortedTypeString($Bond)}; 
				
		}
		
		my @ChargeToApplyToBond = ( );
				
		if(@$PotentialRefinements) {
								
			foreach my $RefiningIncrement (@$PotentialRefinements) {
												
			  if($RefiningIncrement->[2]->getLookUpTable->AreAllNodesSharedBetween([$Bond->getPrimaryLookUpTable])) {
				
			#	  print "Refined: " . $Bond->getPrimaryName . " " . $Bond->getPrimaryType . " " . $Bond->getSecondaryName . " " . $Bond->getSecondaryType . " ".  $RefiningIncrement->[2]->getString . "\n";
								
				  @ChargeToApplyToBond = ($RefiningIncrement->[0], $RefiningIncrement->[1]); last;
			
			  }	
			
			  elsif($RefiningIncrement->[2]->getLookUpTable->AreAllNodesSharedBetween([$Bond->getSecondaryLookUpTable])) {
			
			#	  print "Refined: " . $Bond->getPrimaryName . " " . $Bond->getPrimaryType . " " . $Bond->getSecondaryName . " " . $Bond->getSecondaryType . " " . $RefiningIncrement->[2]->getString . "\n";
				  
				 
			    @ChargeToApplyToBond = ($RefiningIncrement->[1], $RefiningIncrement->[0]); last;
			  }
			}
		}
		
	
		
		if($#ChargeToApplyToBond == -1) {
			if(exists $IncrementLookUpTable{$Bond->getPrimaryType . "-" . $Bond->getSecondaryType})     { @ChargeToApplyToBond = @{$IncrementLookUpTable{$Bond->getPrimaryType . "-" . $Bond->getSecondaryType}} }
			elsif(exists $IncrementLookUpTable{$Bond->getSecondaryType . "-" . $Bond->getPrimaryType})  { @ChargeToApplyToBond = @{$IncrementLookUpTable{$Bond->getSecondaryType . "-" . $Bond->getPrimaryType}}; @ChargeToApplyToBond = ($ChargeToApplyToBond[1], $ChargeToApplyToBond[0]) }
		}
		
		if($#ChargeToApplyToBond == -1 && $SubstituteIncrements) { 
			my $SubstituteIncrement = $AtomTypeSubstituter->FindIncrementSubstitutionForBond($Bond, @$Increments);
			
			if(!$SubstituteIncrement && $Bond->getPrimaryType eq $Bond->getSecondaryType) {
				
				$SubstituteIncrement = [0,0,0,0];
				
			}
			
			elsif(!$SubstituteIncrement) { 
			  
			  $ChargeComplete = 0;
			
			  push @NoAvailableIncrements, $Bond;
					
			  next;
				
			}
			
		  @ChargeToApplyToBond = ($SubstituteIncrement->[2], $SubstituteIncrement->[3]);
		}
		
		#No Substitution 
		if($#ChargeToApplyToBond == -1) { 
			
			return -2 if !$ExitifNotCharged;
			
		  print "FATAL ERROR!: EXITING \n";
			print "the increment between atoms " . $Bond->getPrimaryName . " " . $Bond->getPrimary->Stringify . " (" . $Bond->getPrimaryType . ") and " . $Bond->getSecondaryName . " " . $Bond->getSecondary->Stringify . " (" . $Bond->getSecondaryType . ") does not exist \n";
			print "If you would like the closest match please use the -SubCharge argument\n";
		  exit 1;	
		
		}
		
	  print sprintf("%-7s", $Bond->getPrimaryName) . " " . sprintf("%-7s", $Bond->getPrimaryType)  . " "  . sprintf("%-7s", $Bond->getSecondaryName) . " " . sprintf("%-7s", $Bond->getSecondaryType) . " " . sprintf("%-7s", $Bond->getType)  . sprintf("%+9s", $ChargeToApplyToBond[0]) . " - " . sprintf("%+9s", $ChargeToApplyToBond[1]) . "\n";
    

		$Bond->ApplyIncrement(@ChargeToApplyToBond);
		
	}


  if(@NoAvailableIncrements) {
	 
	  return -2 if !$ExitifNotCharged;
	
	  print "The Following Bonds Could Not be Charged!\n";
	 
	  print sprintf("%-7s", "Name") . " " . sprintf("%-7s", "Type")  . " "  . sprintf("%-7s", "Name") . " " . sprintf("%-7s", "Type") . " " . sprintf("%-7s", "Bond Type")  . "\n";
	  
	
	  foreach my $Bond (@NoAvailableIncrements) {
				
		  print sprintf("%-7s", $Bond->getPrimaryName) . " " . sprintf("%-7s", $Bond->getPrimaryType)  . " "  . sprintf("%-7s", $Bond->getSecondaryName) . " " . sprintf("%-7s", $Bond->getSecondaryType) . " " . sprintf("%-7s", $Bond->getType)  . "\n";
		  
	
	  }

    croak("Could not Finish Charging of Molecule!");
	
  }

  my %Seen;

  my @Atoms = grep { ! $Seen{$_} ++ } map { $_->ToAtomArray }@$Bonds; 


  #Smooth Charge!
  my %CanBeRounded;

  foreach my $Atom (@Atoms) {
	
	  my $Diff = ($Atom->getCharge - sprintf("%.3f", $Atom->getCharge));
		
		if(-0.000001 < $Diff && $Diff > 0.000001) {
			
			$CanBeRounded{$Atom->getName} = [$Atom, $Diff];
			
		}
		
  }

  foreach my $key (keys %CanBeRounded) {
	
	  last;
	
	  my $value = $CanBeRounded{$key};
	
	  my $Atom = $CanBeRounded{$key}->[0];
	
	  my $MissingOppositeCharge = undef;
	
	  foreach my $Bonded (@{$Atom->getBondedAtoms}) {
		
		  my $BondedDiff = $CanBeRounded{$Bonded->getName};
		  
		  next unless defined $BondedDiff;
		
		  if($value->[1] < 0 && $BondedDiff->[1] > 0 || $value->[1] > 0 && $BondedDiff->[1] < 0) {
			
			  print "made it!\n";
			
		  }
		 
		}

# MJH
#    if(! defined $MissingOppositeCharge) {
#	    $Atom->setCharge( sprintf("%.3f", $Atom->getCharge));
#    }
	
	}
  	
	return 1;
		
}

sub SortedTypeString {
	
	my $Bond = shift;
	
	my @SortedTypes = sort map { $_->getType} $Bond->ToAtomArray;
	
  $Bond->SwapPrimaryAtom if $SortedTypes[0] ne $Bond->getPrimaryType;

  return join(" ", @SortedTypes);
	
}
