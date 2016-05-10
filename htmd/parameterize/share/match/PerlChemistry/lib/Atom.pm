package Atom;

use Data::Dumper;
=head1 NAME

Atom - Representation of an Atom, contains all atomic variables such as charge, bonds, etc

=head1 SYNOPSIS

use Atom;
@ISA = qw(Atom);

=head1 DESCRIPTION

The Atom object represents an Atom as you would think of it in chemistry, it is the "Atom" of PerlChemisty, HAHA, 
anyway it basically stores all relevant information about the Atom it represents such as Charge, what other Atoms its bonded too etc

=head2 EXPORT

Nothing

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

use BaseObject ':all';
use Utilities ':func';

require Exporter;

our @ISA = qw(BaseObject Exporter);

# Items to export into callers namespace by default. Note: do not export
# names by default without a very good reason. Use EXPORT_OK instead.
# Do not simply export all your public functions/methods/constants.

# If you do not need this, moving things directly into @EXPORT or @EXPORT_OK
# will save memory.
our $VERSION = '0.01';

#Exported Variables


#Non Exported Variables

#The acceptable number of bonds a element contains, a double bond counts as 2!
our %BondNumHash = ( D => [1], H => [1], I => [1], F => [1], C => [4], N => [3], O => [2], S => [2,4,6], P => [3,5], CL => [1], BR => [1], AL => [5], FE => [4], SI => [4]);

#This is a knowledge based hash of protonation states of elements that can exists, 0 is neutral
our %AllowedProtonationStates = ( C4 => [0], C3 => [0], C2 => [0],
																	O2 => [0], O1 => [-1, 0], 
																	N4 => [1], N3 => [0, 1], N2 => [0], N1 => [0], 
																	F1 => [0], 
																	CL1 => [0], 
																	BR1 => [0], 
																	I1 => [0], 
																	FE4 => [0],
																	H1 => [0], 
																	D1 => [0],
																	S1 => [-1, 0], S2 => [0], S3 => [0,1], S4 => [0], S6 => [0], 
																	SI4 => [0],
																	P4 => [0], 
																	AL4 => [-1,0] );

our %ElementWeight = ( D =>  "1.008000",  H => "1.008000",  C => "12.01100", N => "14.00700", O => "15.99900", F => "18.99800", P => "30.974000", S => "32.06000", CL => "35.45000", BR => "79.90400", I => "126.90447"  );

# Preloaded methods go here.

=head2 New

Usage: Atom->New;

Arguments:
  $Class: should be 'Atom'
  $Fields: a Hash of parameters that need to be loaded into this Atom object, i.e. [Charge => 0, Name => 'N1']

Synopsis:
  Creates Atom object 

=cut

sub New {
	
	my ($Class, $Fields) = @_;
		
	my $Self = { 
	  _Bonds        						  => [],          #Array of Bond Objects
	  _BondedAtoms   						  => [],          #Array of Atom Objects that are bound to this Atom
	  _Chain	      						  => undef,       #The Chain object that the Atom belongs to
	  _Charge       						  => 0,           #Atomic Partial Charge
	  _Dummy								      => 0,           #Is this atom a RTF dummy in a Topology Patch, i.e. it was not declared but it is referred too see XXX
	  _Element      						  => undef,       #The String element of the Atom if its Carbons it would be "C"
	  _ElementWeight              => 0,           #The Numerical Weight of the element 1.008000 for H etc
	  _Fragment     						  => undef,       
	  _Level         						  => 100,         #Used in Molecular Graphs how many bonds away is this atom from the Head node
	  _LookUpTable  						  => undef,       #The LookUpTable object that contains the Molecular Graph for the chemical space of this Atom
	  _MatchNum										=> -1,          #Index used in Matching Molecular Graphs
	  _MissingCharge 						  => 0,           #The amount of Partial Charge that is required in each bond increment rule to return a charge atom to neutral 
	  _Molecule	    						  => undef,       #The Molecule object that the Atom belongs to 
	  _Name                       => undef,       #The name of the atom in string form
	  _NameChanged								=> 0,           #A flag used to keep track of renaming atoms
	  _NoRing                     => 0,           #When building Molecular Graphs keeps track of weather a group of Atoms are ever in a ring
	  _Num          						  => undef,       #The number of the atom, such as in PDB files 
	  _NumOfBonds   						  => 0,           #The Count of Bond Objects found in _Bonds
	  _NumOfNotSolvedBonds        => 0,           #The Count of Bonds Objects that do not have their Increments solved doing the solving process
	  _Occupancy                  => undef,       #The PDB Occupancy
	  _TempFactor                 => undef,       #The PDB TempFactor
	  _PlanarRingAtom							=> 0,           #Boolean Flag for whether an Atom is within a Planar ring
	  _ProtonationState           => 0,           #The Numerical value of the atomic protonation state, i.e. +1 for N in NH4+
	  _Resonance                  => 0,           #Is this atom within a resonant charge structure where its increment will be averaged
	  _Rings                      => [],          #Array of Rings, each ring is just the size of the ring, i.e. 6 for a 6-membered ring
	  _RingAtom										=> 0,           #Boolean Flag for whether an Atom is within any Ring
	  _State                      => undef,       #The State of an atom which is what is used in TypeString declarations which is defined as Atomic Element "." Number of Bonds
	  _SumofBondTypes             => 0,           #The Sum of Bond Types, if an Atom is within 2 singles bonds and 1 double it would be 4
    _Type                       => undef,       #The CHARMM atom type in string form
	  _UniqueString               => undef,	      #The Unique Indenifier to distinguish atoms
	  _X            						  => 0,           #the X coordinate in cartesian coordinates
	  _Y                          => 0,           #the Y coordinate in cartesian coordinates
	  _Z                          => 0,           #the Z coordinate in cartesian coordinates
	};

  if(defined $Fields) {
	  
	  while( my ($AtomField, $Value) = each(%$Fields)) {
		  
		  confess("$AtomField is not a Field of the Atom Object:\n") if ! exists $Self->{"_" . $AtomField};
		
		  $Self->{"_" . $AtomField} = $Value;
				
	  }
	
  }

  bless $Self, $Class;
  return $Self;
	
}

=head2 AmountofChargeToAllowedProtonationState

Usage: $AtomInstance->AmountofChargeToAllowedProtonationState($AtomProtonationState);

Arguments:
  $AtomProtonationState: the protonation state of atom if not specified will use the internal variable _ProtonationState

Synopsis:
  Returns the difference in charge from an allowed Protonation State for $AtomInstance, will only return a value if the 
  Element + NumberOfBonds is possible, i.e. C with one bond is not possible thus will flag a error, These errors can be turned off
  by setting ExistIfNotInitiated to 0 in command line

=cut

{ 
	
	#This hash keeps track of previously calculated values based on the atomic element, number of bonds and protonation state
	our $AllowedProtonationStateHash;

	sub AmountofChargeToAllowedProtonationState {

    my ($Self, $AtomProtonationState) = @_;

    $AtomProtonationState = $Self->getProtonationState unless defined $AtomProtonationState;

    return $AllowedProtonationStateHash->{$Self->getElement . $Self->getNumOfBonds . " " . $AtomProtonationState} if exists  $AllowedProtonationStateHash->{$Self->getElement . $Self->getNumOfBonds . " " . $AtomProtonationState};

    my $ElementAllowedProtonationStates = $AllowedProtonationStates{$Self->getElement . $Self->getNumOfBonds};
    my $BestScore = 10;


#    return -9 if ! defined @$ElementAllowedProtonationStates && $Parameters->getExitifNotInitiated == 0;
#    croak("Allowed Prontation Not defined for " . $Self->getElement . $Self->getNumOfBonds . " for " . $Self->getName . " in " . $Self->getMolecule->getName  ) if  ! defined @$ElementAllowedProtonationStates && $Parameters->getCheckAtomProtState == 1;   
#
    return -9 if ! @$ElementAllowedProtonationStates && $Parameters->getExitifNotInitiated == 0;
    croak("Allowed Prontation Not defined for " . $Self->getElement . $Self->getNumOfBonds . " for " . $Self->getName . " in " . $Self->getMolecule->getName  ) if  ! @$ElementAllowedProtonationStates && $Parameters->getCheckAtomProtState == 1;   

    foreach (@$ElementAllowedProtonationStates) { $BestScore = $AtomProtonationState - $_ if abs($AtomProtonationState - $_) < abs($BestScore) }
	  
		$AllowedProtonationStateHash->{$Self->getElement . $Self->getNumOfBonds . " " . $AtomProtonationState} = $BestScore;
		
		return $BestScore;
		
	}
	
	
}

=head2 calculateProtonationState

Usage: $AtomInstance->calculateProtonationState($SumofBondTypes);

Arguments:
  $SumodBondTypes: the sum of bond types i.e. 4 for Carbon is neutral and would yield a result of 0, if nothing is specified the in the
  interal variable _SumofBondTypes will be used

Synopsis:
  Calculates the Protonation state of Atom, if one wants to update the internal _ProtonationState variable:

  $AtomInstance->setProtonationState($AtomInstance->calculateProtonationState);

=cut

{

  #Store results of calculateProtonationState in this hash for reuse
  our $ProtonationStateHash;

  sub calculateProtonationState {
	
	  my ($Self, $SumofBondTypes) = @_;
	
	  $SumofBondTypes = $Self->getSumofBondTypes unless defined $SumofBondTypes;
	
	  my $Element = $Self->getElement;
	
	  if (exists $ProtonationStateHash->{$SumofBondTypes . $Element}) {
		  
		  return $ProtonationStateHash->{$SumofBondTypes . $Element};
		
	  }
	
	  my $ElementBondNumHash = $BondNumHash{$Element};
	  my $BestScore = 10;
	
	  #Find difference between acceptable protonation and current protonation state
	  foreach (@$ElementBondNumHash) { $BestScore = $SumofBondTypes - $_ if abs($SumofBondTypes - $_) < abs($BestScore) }
	
	  $ProtonationStateHash->{$SumofBondTypes . $Element} = $BestScore;
	
	  return $BestScore;
	
  }

}


=head2 Copy

Usage: $AtomInstance->Copy;

Arguments:

Synopsis:
  Copies an Atom object, all fields are copied except Bonds, BondedAtoms, SumofBondedTypes, NumOfBonds, and NumOfNotSolvedBonds
  In addition also generates a new UniqueString indentifier for the atom

=cut

sub Copy {

  my $Self =shift;	
	
	my $CopiedAtom = Atom->New;
	
	my %DoNotCopy = ( _Bonds => 1, _BondedAtoms => 1, _SumofBondTypes => 1, _NumOfBonds => 1, _NumOfNotSolvedBonds => 1, _Rings => 1 ); 
	
	while ( my ($Field, $Value) = each(%$Self)) {
		
		next if exists $DoNotCopy{$Field};
		
		$CopiedAtom->{$Field} = $Value;
		
	}
	
	push @{$CopiedAtom->getRings}, $_ foreach (@{$Self->getRings});
		
	$CopiedAtom->setUniqueString(GenerateRandomString());
	
	return $CopiedAtom;
}

=head2 AddBond

Usage: $AtomInstance->AddBond($Bond);

Arguments:
  $Bond: The bond object that involves $AtomInstance

Synopsis:
  Adds a Bond object to internal Bond array, only bonds that contain this atom object should be added, (i.e. Bonds that describe this atoms interaction with another atom) 

=cut

sub AddBond {
	
	my($Self, $Bond) = @_;
	
	croak("\$Bond is not a Bond Object in Atom Object:\n" . $Self->DebugInfo) if  ref $Bond ne 'Bond';
	
	croak("This Bond does not contain the Atom Instance that is being added too in Atom Object:\n" . $Self->DebugInfo) if $Bond->getPrimary ne $Self && $Bond->getSecondary ne $Self;
	
	my $SelfBonds = $Self->getBonds;
	
	push(@$SelfBonds, $Bond);
	
	my $SelfBondedAtoms = $Self->getBondedAtoms;
	
	$Bond->getPrimary eq $Self ? push(@$SelfBondedAtoms, $Bond->getSecondary) : push(@$SelfBondedAtoms, $Bond->getPrimary);
		
  #Update Internal Variables for Atom object
	$Self->setNumOfBonds($Self->getNumOfBonds + 1);
	$Self->setNumOfNotSolvedBonds($Self->getNumOfNotSolvedBonds + 1) if ! $Bond->IsIncrementSolved;
	$Self->setSumofBondTypes($Self->getSumofBondTypes + $Bond->getType);
	$Self->setState($Self->getElement . "." . $Self->getNumOfBonds);	
	$Self->setProtonationState($Self->calculateProtonationState);
	
}

=head2 getCartesianCoordinates

Usage: $AtomInstance->getCartesianCoordinates;

Arguments:

Synopsis:
  A wraper function to return [X,Y,Z] coordinates of $AtomInstance

=cut

sub getCartesianCoordinates {
	
	my $Self = shift;
	
	return [$Self->getX, $Self->getY, $Self->getZ];
	
}

=head2 IsBondedTo

Usage: $AtomInstance->IsBondedTo($Atom);

Arguments:
  $Atom: An Atom object that one wants to see if it is bonded to $AtomInstance

Synopsis:
  Checks the _BondedAtoms interal array of $AtomInstance to see if $Atom is contained in it

=cut

sub IsBondedTo {
	
  my ($Self, $Atom) = @_;

  croak("\$Atom is not a Atom object!\n" . $Self->DebugInfo) if ref $Atom ne 'Atom';

  my $SelfBondedAtoms = $Self->getBondedAtoms;

  foreach (@$SelfBondedAtoms) { return $_ if $_ eq $Atom }

  return 0;

}

=head2 IsBondedTo

Usage: $AtomInstance->IsBondedTo($Atom);

Arguments:
  $Atom: An Atom object that one wants to see if it is bonded to $AtomInstance

Synopsis:
  Checks the _BondedAtoms interal array of $AtomInstance to see if $Atom is contained in it

=cut

sub getBondwithAtom {

  my ($Self, $Atom) = @_;

  return 0 unless $Self->IsBondedTo($Atom);  	
	
	my $Bonds = $Self->getBonds;
	
	foreach my $Bond (@$Bonds) {
		
		if($Bond->getPrimary eq $Atom || $Bond->getSecondary eq $Atom) {
			
			return $Bond;
			
		}
		
	}
	
	croak ("Could Not find the bond between " . $Self->getName . " and" . " " . $Atom->getName);

}


=head2 DoesLookUpTableMatch

Usage: $AtomInstance->DoesLookUpTableMatch($LookUpTable);

Arguments:
  $LookUpTable: A LookUpTable object to compare to $AtomInstances LookUpTabke

Synopsis:
  Compare $LookUpTable to $AtomInstances LookUpTable instance, if it returns 1 then $LookUpTable accurately represents $AtomInstance's
Chemical Space

=cut

sub DoesLookUpTableMatch {
	
	my ($Self, $LookUpTable) = @_;
	
	croak("\$LookUpTable is not a LookUpTable object") if ref $LookUpTable ne 'LookUpTable';
	
	my @TableArray = ( $LookUpTable );
	
	return $Self->getLookUpTable->AreAllNodesSharedBetween(\@TableArray);
	
}

=head2 Initiate

Usage: $AtomInstance->Initiate;

Arguments:

Synopsis:
  Setups of Atom object, every Atom should call this function!

=cut

sub Initiate {
	
	my $Self = shift;
	
	if(! defined $Self->getElement) {
		
		croak("Neither Type nor Name is set cannot set up Initiate this Atom\n" . $Self->DebugInfo) if !$Self->getType and !$Self->getName;
		
		my $ElementString;
		
		if($Self->getType && $Self->getName) {
			
			my @TypeCharArray = split "", $Self->getType;
			my @NameCharArray = split "", $Self->getName;
			
			my $Min = @TypeCharArray >= @NameCharArray ? @NameCharArray : @TypeCharArray;
			
			foreach (0 .. $Min-1) {
				
				if($TypeCharArray[$_] eq $NameCharArray[$_] ) { $ElementString .= $TypeCharArray[$_] }
				else { last;}
				
			}
						
		  $ElementString = $Self->getType if !$ElementString;
			
		}
		
		else {
			
			$ElementString = !$Self->getType ? uc($Self->getName) : uc($Self->getType);
			
		}	
		
		#Catch names that start with digits	
		if($ElementString =~ /^\d+(\w+)/) {
			
			$ElementString = $1;
			
		}	
			
	  my $SelfElement = $ElementString !~ /^(?:AL|BR|CL|SI|FE)/ ? substr($ElementString, 0, 1) : substr($ElementString, 0, 2);
				
	  $Self->setElement($SelfElement);
  }

  else {
		
	  $Self->setElement(uc($Self->getElement));  
	
  }

  if(! defined $Self->getName && defined $Self->getElement && defined $Self->getNum) {
		
	  $Self->setName($Self->getElement . $Self->getNum);
	  
  }

  if(! defined $Self->getState) {
	
	  my $SelfState = $Self->getElement . "." . $Self->getNumOfBonds;
	
	  $Self->setState($SelfState);
		
	}		
	
	croak "Something is horribly wrong in Atom Initiation as the Atomic Element or State is not defined!" if ! defined $Self->getElement || ! defined $Self->getState;
	
	my $UniqueString = GenerateRandomString();

  #Check to make sure everything is okay!
   
	$Self->setElementWeight($ElementWeight{$Self->getElement});
	$Self->setUniqueString($UniqueString);
	
	return 1;
		
}


=head2 RemoveBond

Usage: $AtomInstance->RemoveBond($Bond);

Arguments:
  $Bond: A Bond object that $AtomInstance is apart of

Synopsis:
  Removes the bond between $Atom and $AtomInstance if they share a Bond object

=cut

sub RemoveBond {
	
	my ($Self,$Bond) = @_;
	
	my $SelfBonds = $Self->getBonds;
	my $SelfBondedAtoms = $Self->getBondedAtoms;
	
	foreach (0 .. scalar(@$SelfBonds)) {
		
		next if $SelfBonds->[$_] ne $Bond;
			
		$Self->setNumOfBonds($Self->getNumOfBonds - 1);		
		$Self->setNumOfNotSolvedBonds($Self->getNumOfNotSolvedBonds - 1) if ! $Bond->IsIncrementSolved;
				
		$Self->setSumofBondTypes($Self->getSumofBondTypes - $SelfBonds->[$_]->getType);
		$Self->setState($Self->getElement . "." . $Self->getNumOfBonds);
		
		$Self->setProtonationState($Self->calculateProtonationState);
			
		splice(@$SelfBonds, $_, 1);
		splice(@$SelfBondedAtoms, $_, 1);
		
		last;
		
	}	
	
	
}

=head2 setCartesianCoordinates

Usage: $AtomInstance->setCartesianCoordinates($Coordinates);

Arguments:
  $Coordinates: the X,Y,Z coordinates you wish to set the internal X,Y,Z cordinates of $AtomInstance as

Synopsis:
  A wraper function to set [X,Y,Z] coordinates of $AtomInstance

=cut

sub setCartesianCoordinates {
	
	my ($Self, $Coordinates) = @_;
	
	$Self->setX(shift @$Coordinates);
	$Self->setY(shift @$Coordinates);
	$Self->setZ(shift @$Coordinates);
	
}

=head2 Stringify

Usage: $AtomInstance->Stringify;

Arguments:

Synopsis:
  Prints out Atom information

=cut

sub Stringify {
	
	my $Self = shift;
			
	my $SelfState = $Self->getState;
	my $SelfNoRing = $Self->getNoRing;
	my @SelfRings = @{$Self->getRings};
	
	croak("State not defined cannot Stringify this Atom object:\n" . $Self->DebugInfo) if (! defined $SelfState); 
			
	my $String = "";
	
	$String .= "!" if($SelfNoRing);
	
	$String .= $SelfState;
	
	$String .= "%" if $Self->getRingAtom;
	
	$String .= join(",", @SelfRings) if scalar(@SelfRings) > 0;
	
	return $String;
	
}

=head2 GenerateRandomString

Usage: GenerateRandomString()

Arguments:

Synopsis:
  genereates a random string of length 10 including numbers and letters

=cut


sub GenerateRandomString {
	
	my $Length = 10;
	
	my $String = "";
	
	my $randnum = int(rand(127));
	
	foreach (0 .. $Length-1) {
		
		$randnum = int(rand(127));
		
		$String .= chr($randnum);
		
	}
	
	return $String;
}


=head2 getAverageBondAngle

Usage: $AtomInstance->getAverageBondAngle

Arguments:

Synopsis:
  Calculates the average angle between each pair of bonds, this can be used to determine the hybridization of an Atom 

=cut

sub getAverageBondAngle {
	
	my $Self = shift;
			
	my @DiffVectors;
	
	my $BondedAtoms = $Self->getBondedAtoms;
	
	foreach my $BondedAtom (@$BondedAtoms) {
		
		my @DiffVector = map { $Self->getCartesianCoordinates->[$_] - $BondedAtom->getCartesianCoordinates->[$_] } (0 .. 2);
				
		push @DiffVectors, \@DiffVector; 
		
	}
	
	my $Average = 0;
	my $Count = 0;
	
	foreach my $i (0 .. @DiffVectors-1) {
		
		foreach my $j ($i+1 .. @DiffVectors-1) {
								
			my $DotProduct = DotProduct($DiffVectors[$i], $DiffVectors[$j]);
			
			my $Angle = ArcCosine($DotProduct / (Magnitude($DiffVectors[$i])*Magnitude($DiffVectors[$j])));
			
			$Average += $Angle*(180/3.14);
			$Count++;
						
		}
		
	}
	
	return 0 if $Count == 0;
	
	return $Average / $Count;
	
	
}

=head2 DebugInfo (DEPREICATED)

Usage: $AtomInstance->DebugInfo;

Arguments:

Synopsis:
  Prints out all information possible for Atom upon fatal error

=cut

sub DebugInfo {
	
	my $Self = shift;
	
	my $String = "Atomic Details\n";
	
	if(defined $Self->getUniqueString) { $String .= "AtomUniqueString: " . $Self->getUniqueString }
	elsif(defined $Self->getName)		   { $String .= "AtomName: " . $Self->getName}
	elsif(defined $Self->getNum)		   { $String .= "AtomNumber: " . $Self->getNum }
	
	$String .= "\nType: " . $Self->getType if (defined $Self->getType);
	
	if(@{$Self->getBondedAtoms}) {
	  $String .= "\n Bonds: ";
	  $String .= $_->getBondingPartner($Self)->getName . " " . $_->getType . " " foreach (@$Self->getBonds);
	}
	
	return $String;
	
}

1;
__END__
