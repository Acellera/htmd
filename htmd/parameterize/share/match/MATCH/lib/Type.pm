package Type;

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
use Bond;
use LookUpTable;

require Exporter;

our @ISA = qw(MATCHBaseObject Exporter);

# Items to export into callers namespace by default. Note: do not export
# names by default without a very good reason. Use EXPORT_OK instead.
# Do not simply export all your public functions/methods/constants.

# This allows declaration	use BaseObject ':all';
# If you do not need this, moving things directly into @EXPORT or @EXPORT_OK
# will save memory.
our %EXPORT_TAGS = ( 'all'  => [ qw(GetAtomInformationFromString GetRingsFromString)],
											 
										 'vars' => [ qw () ],
										
										 'func' => [ qw ()],
										
										 'test' => [ qw (GetAtomInformationFromString GetRingsFromString)]);
							

our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

our @EXPORT = qw(
	
);

our $VERSION = '0.01';

#Exported Variables


#Non Exported Variables


# Preloaded methods go here.

=head2 New

Usage: Type->New;

Arguments:
  $Class: should be 'Type'
  $Name: Name of Type 
  $BondNumber: the number of bonds that first
  $String: MATCH string that represents Type

Synopsis:
  creates a new Type object

=cut

sub New {
	
	my ($Class, $Name, $BondNumber, $String) = @_;
	
	confess("New Types are created Type->(TypeName, BondNumber, TypeString)") unless defined $Name && defined $BondNumber && defined $String;
	
	my $Self = {
		
		_Atoms        => [], #Treated like an molecule, these dummy atoms will be used to store the chemical environment of this atomtype
		_BondNumber   => $BondNumber,
		_LookUpTable  => undef,
	  _Name         => $Name,
	  _String       => $String,
		
	};
	
	bless $Self, $Class;
	
	return $Self;
	
	
}

=head2 Initiate

Usage: $TypeInstance->Initiate;

Arguments:

Synopsis:
  Setups new Type object needs to be called!

=cut

sub Initiate {
	
	my $Self = shift;
	
	my $String = $Self->getString;
  
  #Break up Type string into an array, each element should represent 1 atom
  my @ArrayedTypeString = split(/([\( \{ \) \}])/, $String);

  #Initiate the first atom which should be contained in the first element of the ArrayedTypeString.
  my $HeadAtom = Atom->New;
  my $HeadAtomString = $ArrayedTypeString[0];

  GetAtomInformationFromString($HeadAtom, $HeadAtomString);

	my $HeadAtomName = $HeadAtom->getElement . 0;
  
  #Set up the internal variables of the Headatom
  $HeadAtom->setLevel(0);
  $HeadAtom->setName($HeadAtomName);
  $HeadAtom->setNum(0);

  $Self->AddAtom($HeadAtom);

  $Self->FindChildren($HeadAtom, @ArrayedTypeString);

  my $Atoms = $Self->getAtoms;

  #Fix State and num of bonds due to how $Atom->AddBond works
  foreach my $Atom (@$Atoms) {
	
	  my $TypeStringNumOfBonds = $Atom->getNumOfBonds - @{$Atom->getBondedAtoms};

    $Atom->setNumOfBonds($TypeStringNumOfBonds);
    $Atom->setSumofBondTypes($TypeStringNumOfBonds);

    if($Atom->getNoNumofBonds == 1) {
	
	    $Atom->setState($Atom->getElement);
	
	  }
	
	  else {

      $Atom->setState($Atom->getElement . "." . $TypeStringNumOfBonds);

    }

    $Atom->setState($Atom->getExclude ? "[^" . $Atom->getElement . "]" : $Atom->getState);
	
  }

 foreach my $Atom (@$Atoms) { $Atom->setUniqueString($Atom->getName) }

  my $SelfLookUpTable = LookUpTable->New;

  $SelfLookUpTable->Initiate($HeadAtom);

  $SelfLookUpTable->Update();

  $Self->setLookUpTable($SelfLookUpTable);
		
}

=head2 AddAtom

Usage: $TypeInstance->AddAtom($Atom);

Arguments:
  $Atom: Type Atom that needs to be added to internal array of Type Object

Synopsis:
  Adds Type Atom to internal array of Type Object

=cut

sub AddAtom {
	
  my ($Self, $Atom) = @_;

  my $Atoms = $Self->getAtoms;

  push @$Atoms, $Atom;	
	
}

=head2 FindChildren

Usage: $TypeInstance->FindChildren($Atom, @ArrayedTypeString);

Arguments:
  $Atom: Type Atom that needs to be added to internal array of Type Object
  @ArrayedTypeString:

Synopsis:
  Adds Type Atom to internal array of Type Object

=cut

sub FindChildren {
	
	my ($Self, $CurrentAtom, @ArrayedTypeString) = @_;
	
	#Level describes how deep into the string we are by counting "(" and ")"s, |: is the location in the array
	#EXAMPLES: |C.3(C.3) = 0
	#           C.3(|C.3) = 1
	#           C.3(C.3)| = 0
	
	my $Level = 0;
	
	foreach my $i (0 .. $#ArrayedTypeString) {
		
		if($ArrayedTypeString[$i] eq "(" || $ArrayedTypeString[$i] eq "{" )   { $Level++ }
		elsif($ArrayedTypeString[$i] eq ")" || $ArrayedTypeString[$i] eq ")") { $Level-- }
		
	  #If we are in level 1 meaning that this is now a child string of level 0
		#C.3(N.2), N.2 is a child of C.3
		if($Level == 1 && $ArrayedTypeString[$i] =~ /\w+/) {
			
		  my $ChildAtom = Atom->New;
			my $ChildAtomLevel = $CurrentAtom->getLevel + 1;
			my $ChildAtomNum = @{$Self->getAtoms} + 1;
			my $ChildAtomString = $ArrayedTypeString[$i];
			
			my $BondOrder = GetAtomInformationFromString($ChildAtom, $ChildAtomString);
			
			#Name is not really important just helps distingiush between atoms
			my $ChildAtomName = $ChildAtom->getElement . $ChildAtomNum;

			#ChildArrayTypeString is going to be the Array used to determine the Children of the Child
			# C.3(N.2(C.1)), C.1 is a child of this N.2 which is a child of C.3
			my @ChildArrayedTypeString;			
			my $ChildLevel = 1;
			
			push @ChildArrayedTypeString, $ArrayedTypeString[$i];
			
			for(my $j = $i+1; $j < @ArrayedTypeString; $j++) {
				
			  if($ArrayedTypeString[$j] eq "(" || $ArrayedTypeString[$j] eq "{" )   { $ChildLevel++ }
				elsif($ArrayedTypeString[$j] eq ")" || $ArrayedTypeString[$j] eq ")") { $ChildLevel-- }	
				
				push @ChildArrayedTypeString, $ArrayedTypeString[$j];
				
				last if $ChildLevel == 0;			
				
			}
			
			#Stores ChildAtoms information into its internal variables.
			$ChildAtom->setLevel($ChildAtomLevel);
			$ChildAtom->setName($ChildAtomName);
			$ChildAtom->setNum($ChildAtomNum);
						
			#Create a Bond object between CurrentAtom and its ChildAtom
			my $ParentChildBond = Bond->New($CurrentAtom, $ChildAtom, $BondOrder);

			$Self->AddAtom($ChildAtom);

			#Recursive, finds Children in Child ArrayedTypeString
			$Self->FindChildren($ChildAtom, @ChildArrayedTypeString);
			
		}
				
	}
	
}

=head2 GetAtomInformationFromString

Usage: Type->GetAtomInformationFromString;

Arguments:
  $Atom: the atom object which corresponds to the string
  $String: the portion of the MATCH string that represents this atom

Synopsis:
  Extracts all of the information about an atom from its MATCH string representation

=cut

sub GetAtomInformationFromString {
	
	my ($Atom, $String) = @_;
	
	my @AtomRings = GetRingsFromString($String);
  my $AtomState;
  my $AtomNoRing = 0;
  my $Charge = 0;
  my $Exclude = 0;
  my $DoubleBondToParent = 0;
  

  $Exclude = 1 if $String =~ /\^/;

  $AtomNoRing = 1 if $String =~ /\!/;

  $DoubleBondToParent = 1 if $String =~ /\=/;

  $AtomState = $1 if $String =~ /=?!?(\S+\.\d|\w+)/;

  $Charge = 1 if $String =~ /\+/;

  $Charge = -1 if $String =~ /\-/;

	#If the Atom state cannot be determined something is wrong with the Type String format, exit out!
  croak("MATCH string not in correct format! : Element . NumofBonds") if !$AtomState;
		
  my ($AtomElement, $NumofBonds) = split(/\./, $AtomState);
  my $NoNumofBonds = 0; 

  #Handles (C) instead (C.1)
  if (! $NumofBonds) {
	
	  $NumofBonds = 1;
	  $NoNumofBonds = 1;	
	
  }

  croak("Incorrect MATCH String format for $String") unless defined $NumofBonds && $AtomElement;

  croak("NoRing and Rings were both present for this Atom cannot be") if $AtomNoRing == 1 && @AtomRings; 

  #Set up the interal variables of the Atom object for the atom.
  $Atom->setElement($AtomElement);
  $Atom->setNoRing($AtomNoRing);
  $Atom->setNumOfBonds($NumofBonds);
  $Atom->setRings(\@AtomRings);
  $Atom->setMissingCharge($Charge);

  if($String =~ /%/) {
	
	  $Atom->setRingAtom(1);
	
  }

  $Atom->setState($AtomState);
  $Atom->{_NoNumofBonds} = $NoNumofBonds;
  $Atom->{_Exclude} = $Exclude;

  return 1 + $DoubleBondToParent;

}

=head2 GetRingsFromString

Usage: Type->GetRingsFromString;

Arguments:
  $Atom: the atom object which corresponds to the string
  $String: the portion of the MATCH string that represents this atom

Synopsis:
  Extracts the rings from a MATCH string (i.e. C.3%6, would return (6). and N.3%5,6 would return (5 6))

=cut

sub GetRingsFromString {
	
	my $String = shift;
	
	my @ArrayedString = split(/\%/, $String);

	return ( ) if $#ArrayedString == 0;
	
	my @Rings = split(/\,/, $ArrayedString[1]);
		
	return @Rings;
	
}

=head2 DoesLookUpTableMatch

Usage: $TypeInstance->DoesLookUpTableMatch($LookUpTable);

Arguments:
  $LookUpTable: A LookUpTable object to compare to $TypeInstances LookUpTabke

Synopsis:
  Compare $LookUpTable to $TypeInstances LookUpTable instance, if it returns 1 then $LookUpTable accurately represents $TypeInstances's
Chemical Space

=cut

sub DoesLookUpTableMatch {
	
	my ($Self, $LookUpTable) = @_;
	
	confess("\$LookUpTable is not a LookUpTable object") if ref $LookUpTable ne 'LookUpTable';
	
	my @TableArray = ( $LookUpTable );
	
	return $Self->getLookUpTable->AreAllNodesSharedBetween(\@TableArray);
	
}


1;
__END__
