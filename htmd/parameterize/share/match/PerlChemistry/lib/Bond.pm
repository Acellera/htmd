package Bond;

=head1 NAME

Bond - describes a chemical bond between two atoms

=head1 SYNOPSIS

use Bond;
@ISA = qw(Bond);

A cool thing to note is that once a Bond is initiated it can directly access information from either of the two atoms it contains:

$BondInstance->getPrimaryCharge is equivalent to $BondInstance->getPrimary->getCharge if one would like to get the partial charge of the primary atom  

=head1 DESCRIPTION

Bond objects represent a chemical bond between two atoms, for PerlChemistry they are basically wraper objects that allow access to both Atoms in the bond
for quick manipulations.

=head2 EXPORT

Nothing

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
use Atom;

require Exporter;

our @ISA = qw(BaseObject Exporter);

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

Usage: Bond->New($PrimaryAtom, $SecondaryAtom, $Type);

Arguments:
  $Class: should be 'Bond'
  $PrimaryAtom: The first atom object that will take part in this bond
  $SecondaryAtom: The second atom object that will take part in this bond
  $Type: single (1), double (2), or triple (3)

Synopsis:
  Creates Bond object 

=cut

sub New {
	
	my ($Class, $PrimaryAtom, $SecondaryAtom, $Type) = @_;

  my $Self = { 
	  _Increment	 => [],									#Contains the Bond Increment after it is determined
	  _Primary     => $PrimaryAtom,       #The first atom object
	  _Secondary   => $SecondaryAtom,     #The second atom object
	  _Type        => $Type,              #Whether the bond is a single(1), double(2) or triple(3)
	  _Molecule    => undef,              #What molecule object the bond belongs too
	  _Applied     => 0,                  #A boolean on whether an Increment has been applied to this bond doing atom charging
	  
	};
	
  bless $Self, $Class;

  #Assumes that this bond will be part of the same molecule that the primary atom is apart this may pose a problem at some point
  $Self->setMolecule($PrimaryAtom->getMolecule); 

  #Adds itself to both atoms, so they are aware of it exists and can update their interal variables
  $PrimaryAtom->AddBond($Self);
  $SecondaryAtom->AddBond($Self);

  return $Self;
}

=head2 ApplyIncrement

Usage: $BondInstance->ApplyIncrement(@ChargeToApply);

Arguments:
  @ChargeToApply: An array of length 2, the first number will be added to the atomic charge of the first atom and the second number will be added to the second

Synopsis:
  Adds the charge from the Set increment to Atoms involved in Bond object

=cut

sub ApplyIncrement { 
  
 	my ($Self, @ChargeToApply) = @_;

  croak("\$ChargeToApply not defined cannot ApplyIncrement") if ! @ChargeToApply;
	
	$Self->getPrimary->setCharge($Self->getPrimaryCharge + $ChargeToApply[0]);
	$Self->getSecondary->setCharge($Self->getSecondaryCharge + $ChargeToApply[1]);
	
	my $ChangeInNum = $Self->IsIncrementSolved ? -1 : 1;
		
	$Self->getPrimary->setNumOfNotSolvedBonds($Self->getPrimaryNumOfNotSolvedBonds + $ChangeInNum );
	$Self->getSecondary->setNumOfNotSolvedBonds($Self->getSecondaryNumOfNotSolvedBonds + $ChangeInNum);
	
	#$Self->setApplied(1);

}

=head2 getBondingPartner

Usage: $BondInstance->getBondingPartner($Atom);

Arguments:
  $Atom: on of the atoms taking part in $BondInstance, returns other Atom

Synopsis:
  returns the Atom object bonded to $Atom

=cut

sub getBondingPartner {
	
  my ($Self, $Atom) = @_;
	
	if($Self->getPrimary eq $Atom)      { return $Self->getSecondary }
	elsif($Self->getSecondary eq $Atom) { return $Self->getPrimary   }

  croak("\$Atom is not contained in this Bond");

}

=head2 IsIncrementSolved

Usage: $BondInstance->IsIncrementSolved;

Arguments:

Synopsis:
  return 1 if Bond has been solved and 0 if it has not been solved, this is part of MATCH

=cut

sub IsIncrementSolved {
	
	my $Self = shift;
  
  my $Increment = $Self->getIncrement;

  return @$Increment ? 1 : 0;

}

=head2 MakeThisAtomPrimary

Usage: $BondInstance->MakeThisAtomPrimary($Atom);

Arguments:
  $Atom: on of the atoms taking part in $BondInstance

Synopsis:
  moves $Atom into the primary position, useful for sorting and manipulations in MATCH

=cut

sub MakeThisAtomPrimary {

  my ($Self, $Atom) = @_;

  if($Self->getSecondary eq $Atom) {  $Self->SwapPrimaryAtom; return }

  croak("This Atom is not in this Bond object!\n") if $Self->getPrimary ne $Atom;
	
}

=head2 SwapPrimaryAtom

Usage: $BondInstance->SwapPrimaryAtom;

Arguments:

Synopsis:
  moves secondary atom into primary position, used for organization during MATCH, also moves the increment positions if they exist

=cut

sub SwapPrimaryAtom {
	
	my $Self = shift;
	
	my $TempAtom = $Self->getPrimary;
	$Self->setPrimary($Self->getSecondary);
	$Self->setSecondary($TempAtom);
	
	return if !$Self->IsIncrementSolved;
	
	my $Increment = $Self->getIncrement;
	
	($$Increment[0], $$Increment[1]) = ($$Increment[1], $$Increment[0]);
	
}

=head2 ToAtomArray

Usage: $BondInstance->ToAtomArray;

Arguments:

Synopsis:
  returns an array of both Atom objects in Bond object, useful for doing sorting

=cut

sub ToAtomArray {
	
	my $Self = shift;
	
	return ($Self->getPrimary, $Self->getSecondary);
		
}

=head2 UpdateType

Usage: $BondInstance->UpdateType($NewType);

Arguments:
  $NewType: the bond type that one wishes to have this bond be

Synopsis:
  Changes the bond type of this bond, i.e. if one wishes to make it a double bond $BondInstance->UpdateType(2).

=cut

sub UpdateType {
	
	my ($Self, $NewType) = @_;
	
	my $OldType = $Self->getType;
	
  $Self->setType($NewType);

  $Self->getPrimary->setSumofBondTypes($Self->getPrimarySumofBondTypes + ($NewType - $OldType));
  $Self->getSecondary->setSumofBondTypes($Self->getSecondarySumofBondTypes + ($NewType - $OldType));

  $Self->getPrimary->setProtonationState($Self->getPrimary->calculateProtonationState);
  $Self->getSecondary->setProtonationState($Self->getSecondary->calculateProtonationState);

}

=head2 UnBondAtoms

Usage: $BondInstance->UnBondAtoms;

Arguments:

Synopsis:
  Removes this bond between the Atom objects involved in it

=cut

sub UnBondAtoms {
  
  my $Self = shift;

  $Self->getPrimary->RemoveBond($Self);
  $Self->getSecondary->RemoveBond($Self);

  $Self->UnSetIncrement; #Cleanup!
	
	$Self->DESTROY;

}

=head2 UnSetIncrement

Usage: $BondInstance->UnSetIncrement;

Arguments:

Synopsis:
  removes the _Increment value and remove Charge from Atoms

=cut

sub UnSetIncrement { 
  
  my $Self = shift;
	
	return if !$Self->IsIncrementSolved;
	
	my @Increment = @{$Self->getIncrement};
	
	$Self->setIncrement([]);
	
	$Self->ApplyIncrement(@Increment);
	
	#Correct the addition of not solved bonds from $Self->AppyIncrement
	$Self->getPrimary->setNumOfNotSolvedBonds($Self->getPrimaryNumOfNotSolvedBonds - 1 );
	$Self->getSecondary->setNumOfNotSolvedBonds($Self->getSecondaryNumOfNotSolvedBonds - 1);
	

} 


sub AUTOLOAD {
  my ($Self, $newvalue) = @_;
  my $type = ref($Self) 
             or print "$Self is not an object\n";


  if($AUTOLOAD =~ /DESTROY/ ) { return; } 

  my ($operation, $attribute) = ($AUTOLOAD =~ /(get|set)(\w+)$/);
  
  #print $operation . " " . $attribute . "\n";

	confess("Method name $AUTOLOAD in " . ref($Self) . " is not in the recognized form (get|set)_attribute\n") unless ($operation && $attribute);
	
	my ($AtomAttribute, $AtomPosition);

	if (! exists $Self->{"_" . $attribute}) {
		if($attribute =~ /(Primary|Secondary)(\w+)/ ) { 
			$AtomPosition = $1;
		  $AtomAttribute = $2;	
		}
		else { confess("Parameter $attribute does not in exist in " . ref($Self) . "!\n") unless (exists $Self->{"_" . $attribute}) }	
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
					  					 $Self->{"_" . $AtomPosition }->{"_" . $AtomAttribute};
					           };
					
		  use strict "refs";

			return $Self->{"_" . $AtomPosition}->{"_" . $AtomAttribute};	
			
		}
	
	  else {
	    *{$AUTOLOAD} = sub { 
		                   my ($Self) = @_;
				  					   $Self->{"_" . $attribute};
					  				 }; 
	  }
	}


  use strict "refs";

  return $Self->{"_" . $attribute};
}



1;
__END__
