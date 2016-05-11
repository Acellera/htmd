package Molecule;

=head1 NAME

Molecule - Contains the bulk of functionality required to correctly characterized the chemical environment of individual atoms

=head1 SYNOPSIS

use Molecule;
@ISA = qw(Molecule);

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
use Bond;
use BaseObject ':all';
use LookUpTable;
use Utilities ':func';

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

our $DistanceBondingRules = { };

our $ResidueDefinitions = { 
	
	'ALA'  => 'Protein',
	'ARG'  => 'Protein',
	'ASN'  => 'Protein',
	'ASP'  => 'Protein',
	'ASPH' => 'Protein',
	'CYS'  => 'Protein',
	'GLN'  => 'Protein',
	'GLU'  => 'Protein',
	'GLY'  => 'Protein',
	'HIS'  => 'Protein',
	'HSD'  => 'Protein',
	'HSE'  => 'Protein',
	'HSP'  => 'Protein',
	'ILE'  => 'Protein',
	'LEU'  => 'Protein',
	'LYS'  => 'Protein',
	'MET'  => 'Protein',
	'MSE'  => 'Protein',
	'PHE'  => 'Protein',
	'PRO'  => 'Protein',
	'SER'  => 'Protein',
	'THR'  => 'Protein',
	'TRP'  => 'Protein',
	'TYR'  => 'Protein',
	'VAL'  => 'Protein',
	'GLUH' => 'Protein',
	'TIP3' => 'Water',
	'TP3M' => 'Water',
	'HOH'  => 'Water',
	'NA'	 => 'Ion',
 	'SOD'  => 'Ion',
	'MG'   => 'Ion',
	'POT'  => 'Ion',
	'CES'  => 'Ion',
	'CA'   => 'Ion',
	'CAL'  => 'Ion',	
	'CLA'  => 'Ion',
	'ZN2'  => 'Ion',
	'ZN'	 => 'Ion',
	'NCO'  => 'Ion',
	'CD'	 => 'Ion',
	'MN'	 => 'Ion',
	'SR'   => 'Ion',
	'GUA'  => 'Nucleic',
  'ADE'  => 'Nucleic',
  'CYT'  => 'Nucleic',
  'THY'  => 'Nucleic',
  'URA'  => 'Nucleic',
	'G'    => 'Nucleic',
  'A'    => 'Nucleic',
  'C'    => 'Nucleic',
  'T'    => 'Nucleic',
  'U'    => 'Nucleic',
  'RG'   => 'Nucleic',
  'RU'   => 'Nucleic',
  'RA'   => 'Nucleic',
  'RC'   => 'Nucleic',
	'DG'   => 'Nucleic',
  'DA'	 => 'Nucleic',
  'DC'	 => 'Nucleic',
  'DT'	 => 'Nucleic',
  'HEME' => 'Cofactor',
  'NAG'  => 'Factor',
  'MAN'  => 'Factor',
  'BMA'  => 'Factor',
  'SO4'  => 'Factor',
  'GOL'  => 'Factor',
  'IPA'  => 'Factor',

};

our $CommonMissnamedAtoms = {
	
  "OP1"		=> "O1P",
  "OP2"   => "O2P",
	"HO2'"	=> "H2'",
	"H2'"		=> "H2''",
	"2HO"   => "H2'",
	"1H5'"  => "H5''",
	"2H5'"  => "H5'",
#	"H3'"		=> "H3T",
	
};

our $CommonMissnamedResidues = {
	
	'RG'   => 'GUA',
  'RU'   => 'URA',
  'RA'   => 'ADE',
  'RC'   => 'CYT',
	'G'    => 'GUA',
  'U'    => 'URA',
  'A'    => 'ADE',
  'C'    => 'CYT',
  'DG'   => 'GUA',
  'DA'	 => 'ADE',
  'DC'	 => 'CYT',
  'DT'	 => 'THY',
  '5BU'  => 'URA',
  'ASPH' => 'ASP',
  'HIS'  => 'HSD',
  'HIE'  => 'HSD',
  'HIP'  => "HSP",
  'CYX'  => 'CYS',
  'CYM'  => 'CYS',
  'CA'   => 'CAL',
};


# Preloaded methods go here.

=head2 New

Usage: Molecule->New;

Arguments:
  $Class: should be 'Molecule'
  $StructuralHash: The raw data from a file to be loaded into the molecule object

Synopsis:
  creates a new Molecule object

=cut

sub New {
	
	my ($Class, $StructuralHash) = @_;	
	
	my $Self = {
		
		_Atoms          => [ ],       #The Atom Object of the atoms within this Molecule
		_Bonds					=> [ ],       #The Bond Objects of the bonds within this Molecule
		_Chain					=> undef,     #The Chain Object that this molecule is apart of
		_Charge					=> undef,     #The Formal Charge of the molecule 
	  _ConnectorAtom  => undef,     #Which atom forms intermolular bonds with other molecules in the chain such as N and C for Proteins
	  _ConnectTo			=> undef,     #Which atom object does the Connector Atom bond too to connect the molecules
		_IsPatch				=> 0,         #Is this molecule just a Topology Patch
		_Impropers			=> [ ],       #The list of atom names that need an improper angle to keep them within the plane in CHARMM
		_Indentification=> undef,     #Is it a Protein/NA/Ligand etc
		_Failed 				=> 0,         #Keeps track if something did not work properly and explict errors are turned off
		_FileName			  => undef,     #Where was this file loaded in from / or where should it be printed too
		_Name						=> undef,     #Residue/Molecule name
		_Num						=> undef,     #Molecule number
		_Rings					=> [ ],       #Each ring is an array of atoms that make up an Ring
		_SegmentId			=> "",        #The CHARMM PDB Segment ID
		
	};
	
	bless $Self, $Class;
	
	#Structural Hashes are built from loading in Files such as PDB/MOL2/SDF etc
  $Self->BuildFromStructuralHash($StructuralHash) if defined $StructuralHash;
	
	return $Self;
	
}

=head2 BuildFromStructuralHash

Usage: $MoleculeInstance->BuildFromStructuralHash;

Arguments:
  $StructuralHash: should be 'Molecule'

Synopsis:
  Loads data from StructuralHash built from loading in from file such as a PDB, MOL2 or RTF

=cut

sub BuildFromStructuralHash {
	
	my ($Self, $StructuralHash) = @_;
	
	croak("This Hash is not of a Molecule") if delete $StructuralHash->{'Structure'} ne ref($Self);
		
	my $Atoms;
	my $Bonds;
	my $Deletes;
	
	#Transfer information of the StructuralHash into the Molecule object, Atoms, Bonds, and Deletes need addition processing before adding them
	while(my ($Field, $Value) = each(%$StructuralHash)) {
				
		if($Field eq 'Atoms')      { $Atoms = $Value; next } 
		elsif($Field eq 'Bonds')   { $Bonds = $Value; next } 
		elsif($Field eq 'Deletes') { $Deletes = $Value; next }
						
		croak("$Field is not a Field of the Molecule Object:\n") if ! exists $Self->{"_" . $Field};

		$Self->{"_" . $Field} = $Value;
		
	}
	
	#print scalar(@$Atoms) . "\n";
	
	#Is it a protein/nucleic/ligand etc, based on residue name
	$Self->IndentifyMolecule();
	
	#Fix PDB A vs B crystalization thing if A and B are part of Atom names
  my @Atoms = $Self->RemoveDuplicateAtoms($Atoms);
		
	my @AtomObjects;
	my $AtomCount = 1;
	
	my %AtomNames;
  my %AtomPositions;
				
	foreach my $AtomDataHash (@Atoms) {
			
		my $AtomObject = Atom->New($AtomDataHash);
				
		$AtomObject->setMolecule($Self);
		
		$AtomObject->setChain($Self->getChain) if defined $Self->getChain;
		
		$AtomObject->setNum($AtomCount++) if ! defined $AtomObject->{_Num};
		
		#Catch Errors if ExitifInitiated is set to 0
		if($AtomObject->Initiate == 0) {
			
			$Self->setFailed(1); last;
			
		}
		
		if(exists $AtomNames{$AtomObject->getName}) {
			
			#If atoms are named by only element add numbers to make atom names C1 C2 C3 etc instead of C C C 			
			$AtomObject->setName($AtomObject->getElement . $AtomObject->getNum) if $AtomObject->getName eq $AtomObject->getElement && ! ( defined $Bonds || $Self->getIsPatch);
			
		}
						
		$AtomNames{$AtomObject->getName} = 1;
		
		push @AtomObjects, $AtomObject;
		
	}
	

	ReNameAtoms(\@AtomObjects) if $Parameters->getRenameAtoms;
		
	$Self->setAtoms(\@AtomObjects);
	 		
	if(defined $Bonds || $Self->getIsPatch) { 
				
		$Self->setFailed(1) if $Self->UseBondingInformationFromFile($Bonds) == 0;
	
	}
	
	elsif(@$Atoms > 1) { $Self->BondAtomsUsingIdealBondLengths }	 
	
	#This is a Patching Molecule, initiate Deletes internal variable
	$Self->{"_Deletes"} = defined $Deletes ? $Deletes : [] if $Self->getIsPatch == 1;

	
}

=head2 RemoveDuplicateAtoms

Usage: $MoleculeInstance->RemoveDuplicateAtoms($AtomDataHashes);

Arguments:
  $AtomDataHashes: the set of atoms data before they are turned into molecules

Synopsis:
  Weeds out duplicate atoms when there is an A/B atom due to crystalization possibilities in PDBss

=cut

sub RemoveDuplicateAtoms {
	
	my ($Self, $AtomDataHashes) = @_;
	
	my %NameExists;
	my %RemoveAtom;
		
	if($Self->getIndentification ne 'Ligand')  {
				
	  foreach my $Atom (@$AtomDataHashes) {
	
	    if($Atom->{'Name'} =~ /^D/) {

        my $TempName = "H" . substr($Atom->{'Name'},1);

        if($TempName =~ /^(\S+)\s+([A-Z])$/) {
	
				  my $Name = $1;
				
				  push @{$NameExists{$1}}, $Atom;
	
				}
				
				elsif(length($TempName) > 3 && $TempName =~ /^(\S+)([A-Z0-9]+)/) {
				
				  next if length($2) > 1 || $2 =~ /^\d+/;

					my $Name = $1;

					push @{$NameExists{$1}}, $Atom;
					
				}
		
	    }
	
			elsif($Atom->{'Name'} =~ /^(\S+)\s+([A-Z])$/) {

				    my $Name = $1;

			      push @{$NameExists{$1}}, $Atom;

			    }
	
	    #Still not sure about this
	    elsif(length($Atom->{'Name'}) > 3 && $Atom->{'Name'} =~ /^(\S+)([A-Z0-9]+)/) {
		
		    #print $Atom->{'Name'} . " " . $1 . " " . $2 . "\n";
				
				#Exception H and HD12 or HD1 and HD12, hopefully this solves it!
		    next if length($2) > 1 || $2 =~ /^\d+/;
				
		    my $Name = $1;

				push @{$NameExists{$1}}, $Atom;
		
	    }
	
	    else {
			
		    push @{$NameExists{$Atom->{'Name'}}}, $Atom; 
		
		  } 
	
    }

    while (my ($Name, $Atoms) = each %NameExists) {
		
	    #print $Self->getName . " ". $Self->getNum . " " . $Name . " : " .  join(" ", map { $_->{'Name'}} @$Atoms) . "\n" if @$Atoms > 1;
	
	    next if @$Atoms == 1;

      $Atoms->[0]->{'Name'} = $Name;

      shift @$Atoms;

      foreach my $Atom (@$Atoms) {
	
	      $RemoveAtom{$Atom} = 1;
	
      }
	
    }

  }

  my @Atoms = grep { ! exists $RemoveAtom{$_} } @$AtomDataHashes;
		
	return @Atoms;
	
}

=head2 IndentifyMolecule

Usage: $MoleculeInstance->IndentifyMolecule;

Arguments:

Synopsis:
  Using a simple Residue name matching approach, determines what type of molecule this is (Protein/NA/Crysalization Factor/Ion or Ligand)

=cut

sub IndentifyMolecule {
	
	my $Self = shift;
	
	$Self->setName($CommonMissnamedResidues->{$Self->getName}) if exists $CommonMissnamedResidues->{$Self->getName};
			
	my $Identification = exists $ResidueDefinitions->{$Self->getName} ? $ResidueDefinitions->{$Self->getName} : 'Ligand';
	
	$Self->setIndentification($Identification);	
	
}

=head2 ReNameAtoms

Usage: ReNameAtoms($Atoms);

Arguments:
  $Atoms: the atoms of the molecule 

Synopsis:
  If the naming is really weird, this function can rename all the atoms using the naming convention of Element + AtomNumber

=cut

sub ReNameAtoms {
	
	my $Atoms = shift;
	
	my %ElementCount;
	
	foreach my $Atom (@$Atoms) {
		
		my $Count = exists $ElementCount{$Atom->getElement} ? $ElementCount{$Atom->getElement} : 0;
		
		$ElementCount{$Atom->getElement}++;
		
		$Atom->setName($Atom->getElement . ($Count+1));
		
	}
	
}

=head2 Initiate

Usage: $MoleculeInstance->Initiate;

Arguments:

Synopsis:
  Setups Molecule Object, must be called if any of the MATCH functions are to be used!

=cut

sub Initiate {
	
	my $Self = shift;
		
	my $Atoms = $Self->getAtoms;
	
	#Loop Through all atoms and build there LookUpTables (Molecular Graphs) with that atom as the head node		
	foreach my $Atom (@$Atoms) {
		
		return 0 if $Atom->AmountofChargeToAllowedProtonationState == -9;
				
		my $AtomLookUpTable = LookUpTable->New;
		
		$AtomLookUpTable->Initiate($Atom);
	  
	  $Atom->setLookUpTable($AtomLookUpTable);	 
			
	}
	
	#Calculate Formal Charge of Molecule	
	my $TotalCharge = 0;
	
	$TotalCharge += $_->getCharge foreach (@$Atoms);
		
	$Self->FindRings;
	
	$Self->CleanEnsureCorrectBondTypes;		
	#$Self->NewEnsureCorrectBondTypes;
	
	$Self->NewReduceAtomicChargesToLowestPossible;
		
	$Self->FindResonanceStructures;
	
	$_->getLookUpTable->Update foreach (@$Atoms);
	
	return 1;
	
}


=head2 FindRings

Usage: $MoleculeInstance->FindRings;

Arguments:

Synopsis:
  Find Atom Rings in this MoleculeInstance, stores them in internal variable _Rings

=cut


sub SetCustomCharge {
	my $Self = shift;
	my $charge = shift;
	print( "Set charge $charge\n" );
	my $Atoms = $Self->getAtoms;
  my 	$i=0;
	foreach my $Atom (@$Atoms) {
			if ( $i==0 ) {
				$Atom->setCharge( $charge );
				$i=1;
			}
	}   


}

sub FindRings {
	
	my $Self = shift;
	
	#Take only atoms that have more than one bond, not possible for 1 bonded atoms to be in a ring
	my @HeavyAtoms = grep { $_->getNumOfBonds > 1 } @{$Self->getAtoms};
	
	my $NumOfAtoms = scalar(@HeavyAtoms);
	
	my %UsedAtoms;
	my $SelfRings = $Self->getRings;
	
	my %Aromatic = ( 'C.3' => 1, 'N.2' => 1, 'N.3' => 1, 'S.2' => 1, 'O.2' => 1 );
	
	foreach my $Atom (@HeavyAtoms) {
	
	  next if $UsedAtoms{$Atom};
	
	  #Use Levels from $Atom's LookUpTable, this allows more quicker generation of rings, since if $Atom's level is always 0 so to get back
	  #to $Atom the most direct path is to follow the path of lowest level.
	  $Atom->getLookUpTable->SetAtomLevelsUsingNodeLevels;
		
		my @Pathes;
		my @NextLevelPathes;
 		
		push @Pathes, [$Atom];
		
		my $Finished = 0;
		my @Ring;
		
		while(!$Finished) {
	  
	    foreach my $Path (@Pathes) {
		
		    #No path can be longer than 1/2 the size of the molecule and need to turn back to get back to $Atom
		    next if $Path->[-1]->getLevel > ($NumOfAtoms - @$Path + 1);
				
		    #Remove first element of @ArrayedPath which is $Atom, this makes it easier to process
		    my @ArrayedPath = @{$Path};
				shift(@ArrayedPath);
		
		    my $Children = $Path->[-1]->getBondedAtoms;
		    my @AcceptableChildren;
		
		    foreach my $Child (@$Children) {
					
					#If this Child has only one bond (Dead End) or is already included in the list its not worth following that path
					#Keep in mind $Atom cannot be here since it was removed already with the shift statement 
			    next if $Child->getNumOfBonds < 2 || grep { $Child eq $_ } @ArrayedPath;
			
			    if($Child eq $Atom) {
				
				    if(@$Path > 2) { @Ring = @$Path; $Finished = 1} #Success
				    else					 { next;  }											  #There cant be 2 membered rings!
				
			    }
			
			    push @AcceptableChildren, $Child unless $Finished; #This Child is allowable unless we are Finished
									
		    }
		
				my @SortedChildren = sort { $a->getLevel <=> $b->getLevel } @AcceptableChildren; #Find the lowest level out of Children

			  my @BestChildren = grep { $SortedChildren[0]->getLevel == $_->getLevel } @SortedChildren; #Follow all pathes that have the lowest level
		
		    push @NextLevelPathes, [@$Path, $_] foreach (@BestChildren); #Prepare next iteration 
	
	    }
	
	    $Finished = 1 if @NextLevelPathes == 0;
	
	    @Pathes = @NextLevelPathes;
	    @NextLevelPathes = ( );
	
	  }
	
	  if(@Ring) { 
		 
      #Attempt to assign whether a Ring is Aromatic (A) or NotAromatic (N) based on the the Element and Bond count of each atom in ring
      #This might get phased out to the ability lower in this function to find planar rings based on geometry
      my $NotAromatic = 0;
		
		  foreach my $RingAtom (@Ring) {
								
			  unless(exists $Aromatic{ $RingAtom->getState }) {
				
				  $NotAromatic = 1;
				  last;
				
			  }
			
		  }
		
		  my $Aromatic = $NotAromatic ? "N" : "A";
		
		  foreach my $RingAtom (@Ring) {
			
			  $RingAtom->setRingAtom(1);
			
			  $UsedAtoms{$RingAtom} = 1;
		
		    push @{$RingAtom->getRings}, scalar(@Ring) . $Aromatic ;
			 		
		  }
		
		  #Start newer geometry based algorithm to determine aromaticity
		  #Calculates the difference vector between each pair of atoms the ring, and then using 3 of these 
		  #at a time calculates (D1 X D2) * D3 which whould be 0 if all the Vectors are in the same plane
		  my @DiffVector;
		
		  foreach my $i (0 .. @Ring-1) {
			
			  foreach my $j ($i+1 .. @Ring-1) {
			
				  my @DV = map { $Ring[$i]->getCartesianCoordinates->[$_] - $Ring[$j]->getCartesianCoordinates->[$_] } (0 .. 2);

			  	push @DiffVector, \@DV; 
				
			  }
			
		  }
		
		  my $Average = 0;
		  my $Count = 0;
		  
		  for (my $i = 0; $i < @DiffVector-1; $i+=3) {
			
			  $Average += DotProduct(CrossProduct($DiffVector[$i], $DiffVector[$i+1]), $DiffVector[$i+2]) if $i+2 < @DiffVector;
			  
			  $Count++;
			
		  }
	
		  if(abs($Average / $Count) < 0.05) {
			
			  foreach my $Atom (@Ring) {
				
				  $Atom->setPlanarRingAtom(1);
				
			  }
			
		  }
		
		  push @$SelfRings, \@Ring;
		
		}
	
	}
	
}

=head2 CleanCorrectBondTypes

Usage: $MoleculeInstance->CleanEnsureCorrectBondTypes;

Arguments:
  $ConsiderAllGroups: set to 1 if you wish to brute force consider all solutions instead of discarding bad solutions, this will be done automatically if
this function fails to find a suitable solution on its first attempt

Synopsis:
  Fixes all the bonding throughout the molecule makes sure that single / double / triple bonds are correct

=cut


sub CleanEnsureCorrectBondTypes {
	
	my ($Self, $ConsiderAllGroups) = @_;
	
	my $Atoms = $Self->getAtoms;
	my $Bonds = $Self->getBonds;
	
	#Not really atoms with incorrect protonation states more just atoms with non zero protontation
	my @AtomsWithIncorrectProtonationStates = sort { @{$a->getBonds} <=> @{$b->getBonds} } grep { $_->getProtonationState != 0} @$Atoms;
	
	#If are no atoms with protonation states that are non zero, exit out, since its a waste of time 
	return unless @AtomsWithIncorrectProtonationStates;
	
	#Looking for potential bonds to fix, only considering bonds with atoms that can support more than one bond (by element), and contain one atom with
	#a nonzero protonation state since those are likely to change	
	my @BondsWithIncorrectProtonationStates = grep { (! CanAtomOnlySupportOneBond($_->getPrimary) && ! CanAtomOnlySupportOneBond($_->getSecondary)) &&
																									 ($_->getPrimary->getProtonationState != 0 || $_->getSecondary->getProtonationState != 0)  } @$Bonds;
	
	#Quick way of looking up bonds
	my %BondKeys = map {  $_ => join( " ", sort map { $_->getName } $_->ToAtomArray) } @$Bonds;
	
	#Keeps track of which bonds should be investigated in the upcoming algorithm	
	my %BondsWithMissingCharge = map {  $BondKeys{$_} => 1 } @BondsWithIncorrectProtonationStates;
									
  #To speed up this brute force procedure, using arrays instead of hashes to store the current values of each bond type/ atom sum of bonds
  #helps. Using the $Bonds/$Atoms arrays each respective bond or atom will be referenced by its position in its array. Thus the first bond in the $bonds
  #array will be referenced as 0. %BondCountHash keeps track of this number when crawling to a bond for the first time
  my $BondCount = 0;

  my %BondCountHash;

  my @BondValues;

  foreach my $Bond (@$Bonds) {
	
	  push @BondValues, $Bond->getType;
	  
	  $BondCountHash{ $BondKeys{$Bond} } = $BondCount++;
	
  }

  my $AtomCount = 0;

  my %AtomCountHash;

  my @AtomValues;

  foreach my $Atom (@$Atoms) {
	
	  push @AtomValues, $Atom->getSumofBondTypes;
	
	  $AtomCountHash{ $Atom->getName } = $AtomCount++;
	
  }
	
	#Groups keep track of the current bond number of each atom and the current bond type of each bond ( single=1/double=2/triple=3)
	#=> Should probably be called Solutions instead of Groups = 8.8.11 JDY
	my @Groups = ({ Atoms => \@AtomValues, BondValues => \@BondValues });
	
	my $CurrentAtom = shift @AtomsWithIncorrectProtonationStates;
	my @NewBonds = @{$CurrentAtom->getBonds};
	
	my %Visited;
	
	my $Count = 0;
	
	my $Done = 0;
	
	#This algorithm exhustively tries different combinations of 1/2/3 bond types for bonds that can support them to try and find
	#the combination that produces the lowest amount of formal charge for the molecule that is allowed 
	
	while ( !$Done ) {
		
	  #Find next bond to consider
	  my $CurrentBond;
	
	  #Look through NewBond to find the next bond to consider that meets all these conditions
	  while ( @NewBonds) {
				  
		  my $TempBond = shift @NewBonds;
	  	
	    next unless $BondKeys{$TempBond}; #joined residues
	
	    next if exists $Visited{ $BondKeys{$TempBond} }; #Have we seen this bond already?
		
		  #This removes any bond that has H/F/CL etc any atom that can only bond with a single bond, these cannot be changed so not worth considering
	    next if CanAtomOnlySupportOneBond($TempBond->getPrimary) || CanAtomOnlySupportOneBond($TempBond->getSecondary); 
							    
		  $Visited{ $BondKeys{$TempBond} } = 1;
		
		  next unless exists $BondsWithMissingCharge{  $BondKeys{$TempBond} }; #Look at only bonds that have missing charge
		
		  $CurrentBond = $TempBond;
		
		  last;
		
	  }
	
	  #All bonds have been considered exit out
    unless(defined $CurrentBond) {
	
	    #This is a final check to make sure we have looked at each bond with incorrect charge, this helps if the algorithm gets stuck at an end
	    my $LastCheck = 0;
	
	    foreach my $Bond (@BondsWithIncorrectProtonationStates) {
		
		    next if $Visited { $BondKeys {$Bond}};
		
		    $LastCheck = 1; $CurrentBond = $Bond;
		
	    }
	
	    #If we have checked each bond with incorrect protonation states we are done and can exit out     
	    if(!$LastCheck) {
	
	      $Done = 1; last
	
      }
	
    }

    #Lookup Count numbers to be used in the groups, this was save space and computation time!
    my $BondNumber = $BondCountHash{ $BondKeys{$CurrentBond} };
    my $PrimaryAtomNum = $AtomCountHash{ $CurrentBond->getPrimaryName };
    my $SecondaryAtomNum = $AtomCountHash { $CurrentBond->getSecondaryName };
		
		#Crawl to new bonds, these are the bonds we will consider next    
	  push @NewBonds, (@{$CurrentBond->getPrimary->getBonds}, @{$CurrentBond->getSecondary->getBonds});
	   
	  #Get all the different possible bond states given the current bond number of each atom in a bond, see getAllowableBondTypes for more info	   
	  foreach my $Group (@Groups) {
				     
		  $Group->{PossibleNext} = getAllowableBondType($CurrentBond, $Group->{Atoms}->[$PrimaryAtomNum], $Group->{Atoms}->[$SecondaryAtomNum]);

	  }
	
    #Update Groups
		my @NextGroups;

		foreach my $Group (@Groups ) {

	    my $PossibleNexts = $Group->{PossibleNext};

			#If there are no new bonds in this group skip, this should really only happen at the end
		  unless(defined $PossibleNexts) {

				push @NextGroups, $Group; next;

			}

			foreach my $PossibleNext (@$PossibleNexts) {

        my ($ChangeInBondType, $Bond) = @$PossibleNext;

        #If the Change in Type is 0 then don't create a new group just use existing group, this avoids generating exact copies on the current group
			  if($ChangeInBondType == 0) {

					push @NextGroups, $Group; next;

			  }
			
			  #NewGroup will have the same properties as the old group other than this one bond state
			  my $NewGroup = { Atoms => [map { $_} @{$Group->{Atoms}}], BondValues => [map { $_ } @{$Group->{BondValues}}]  };
			
			  $NewGroup->{Atoms}->[$PrimaryAtomNum]   -= $ChangeInBondType;
			  $NewGroup->{Atoms}->[$SecondaryAtomNum] -= $ChangeInBondType; 
			
			  $NewGroup->{BondValues}->[$BondNumber]  = -$ChangeInBondType + $Bond->getType;   
        
			  push @NextGroups, $NewGroup;
			
			}

		}
				 
    $Count++;

    #Finished updating groups with new bond, move them into @Groups for next round
    @Groups = @NextGroups;  

    my $GroupCount = 0;

    #ConsiderAllGroups is only defined if this function has already failed to produce a suitable solution, this is a last ditch effort to get a solution by 
    #considering all solutions/groups, this will significantly slow down the speed but is more reliable for some cases.
    next if defined $ConsiderAllGroups;

    #Because this method is brute force considering all possibilities for single/double/triple bond combos that are allowable, for very large molecule
    #the computational cost gets very large, thus every 5 rounds, the groups with the most missing charge (i.e. least number of correct bond types) are removed
    next if $Count == 0 || $Count % 5 != 0;

    #Calculate the average amount of missing charge per group 
    my $TotalCharge = 0;
		
	  foreach my $Group (@Groups) {
				
		  my $MoleculeCharge = 0;
		
		  foreach my $AtomNumber (0 .. @$Atoms-1) {
			
				$MoleculeCharge += abs($Atoms->[$AtomNumber]->calculateProtonationState($Group->{Atoms}->[$AtomNumber]));

		  }
		
		  $TotalCharge += $MoleculeCharge;
				
		  $GroupCount++;
		
	  }
		 
    my $AvgCharge = @Groups != 0 ? $TotalCharge / @Groups : 0;
		
	  @NextGroups = ();
	
	  #Go through all the groups again and remove ones that are below the current average missing charge
		foreach my $Group (@Groups) {
		 
		  my $MoleculeCharge = 0;

			foreach my $AtomNumber (0 .. @$Atoms-1) {
			
				$MoleculeCharge += abs($Atoms->[$AtomNumber]->calculateProtonationState($Group->{Atoms}->[$AtomNumber]));

		  }
		
		  push @NextGroups, $Group if $MoleculeCharge < $AvgCharge;
		
	  }   
	
	  @Groups = @NextGroups;    
   
	}
	
	#Done with generating all the possiblities now its scoring time, looking for groups that have zero missing charge (every atoms has a possible
	#protonation state) and the minimize the formal charge of the molecule. 
		
	#Final Scoring of Groups, ones with lowest charge
	my $BestScore = 100; my @BestGroups;	
		
	foreach my $Group (@Groups) {
		
	  my $TotalMissingCharge = 0;
		
	  my $MoleculeCharge = 0;

    #Calculate the Total amount of missing charge and the formal charge of molecule		
		foreach my $AtomNumber (0 .. @$Atoms-1) {
			
			$TotalMissingCharge += AllowedProtonationState($Atoms->[$AtomNumber], $Group->{Atoms}->[$AtomNumber]);
						
			$MoleculeCharge += abs($Atoms->[$AtomNumber]->calculateProtonationState($Group->{Atoms}->[$AtomNumber]));
									
		}
		
		#If there is any missing charge then this group is not correct go to the next	
		next unless $TotalMissingCharge == 0;	  		
		
		#Assuming there is no missing charge lets find groups that have the lowest amount of formal charge
		if( abs($MoleculeCharge) < abs($BestScore) ) {
							
		  $BestScore = $MoleculeCharge; 
		
		  @BestGroups = ($Group);
		
	  }
	
	  elsif( abs($MoleculeCharge) == abs($BestScore) ) {
		
		  push @BestGroups, $Group;
		
	  }
			
	}
	
	#No solution was found, try again considering all groups
	return $Self->CleanEnsureCorrectBondTypes(1) if @BestGroups == 0 && ! defined $ConsiderAllGroups;
	
	#Possible Error or for some reason went through the algorithm without needing it	
	if(@BestGroups == 0 && defined $ConsiderAllGroups) {
		
		print "Could not determine a suitable solution for bond types for " . $Self->getName . "! be wary\n";
		
		return 0;
		
	}	
		
	#Keep Charge as close to the outside of the molecule as possible
  #Set in MATCH.pl
	my $AtomicProtonationStates = $Parameters->getAtomicProtonationStates;	
	
	my $FinalGroup; my $FinalScore = -1;
	
	foreach my $Group (@BestGroups) {
		
		my $CurrentScore = 0;
		my $Viable = 1;
				
	  foreach my $AtomNumber (0 .. @$Atoms-1) {
		 
		  #A set protonation state was set by user and this solution does not contain that protonation state	
		  if(defined $AtomicProtonationStates->{$Atoms->[$AtomNumber]->getName} && $AtomicProtonationStates->{$Atoms->[$AtomNumber]->getName} !=  $Atoms->[$AtomNumber]->calculateProtonationState($Group->{Atoms}->[$AtomNumber])) {
						
			  $Viable = 0; last;
			
		  }

      #Look at only atoms that have formal charge
			next if $Atoms->[$AtomNumber]->calculateProtonationState($Group->{Atoms}->[$AtomNumber]) == 0;
			
		  if($Atoms->[$AtomNumber]->getNumOfBonds == 1) { $CurrentScore++ }
		
		  else {
			
			  foreach my $BondedAtom (@{$Atoms->[$AtomNumber]->getBondedAtoms}) {
				
				  if($BondedAtom->getNumOfBonds == 1) { $CurrentScore += 1 / $Atoms->[$AtomNumber]->getNumOfBonds }  
				
			  }
			
		  }

		}
		
		next if !$Viable;
		
    if($CurrentScore > $FinalScore ) {
	
	    $FinalScore = $CurrentScore; $FinalGroup = $Group;
	
    }
		
	}	
	
	if(! defined $FinalGroup && keys %$AtomicProtonationStates) {
		
		croak("Could not find a suitable set of bond orders that would satisfy the protonation state constraint given, please reduce or remove the constraints");
		
	}
	
  #Apply Final Solution	
	foreach my $BondNumber (0 .. @$Bonds-1) {    
						
		$Bonds->[$BondNumber]->UpdateType(  $FinalGroup->{BondValues}->[$BondNumber] );
		
	}
	
	#Success!
	return 1;
		
}

#Functions for CleanEnsureCorrectBondTypes

=head2 getAllowableBondType

Usage: getAllowableBondType($Bond, $PrimarySumOfBondNum, $SecondarySumOfBondNum);

Arguments:
  $Bond: The bond object that you wish to calulate the possible types (i.e Single/Double/Triple) given the bond number of both atoms in the bond
  $PrimarySumOfBondNum: The current bond number for the first atom in the bond, number 4 for Carbon would be neutral
  $SecondarySumOfBondNum: The current bond number for the second atom in the bond

Synopsis:
  Calculates the possible bond types that $Bond can occupy, returns an array, in the form: [[ChangeInBondType1, $Bond], [ChangeInBondType2, $Bond]], where 
the change in bond types are the difference from the current bond type to the new one. Example if the bond is current a single bond, then -1 would be a double.

=cut

{

  #Since there will be many groups in EnsureCorrectBondTypes that will require the same calculation this hash, saves the return value given a certain a combination
  #of function arguments for speedy recall
  our $AllowedBondTypeHash;

  sub getAllowableBondType {

    my ($CurrentBond, $PrimarySumOfBondNum, $SecondarySumOfBondNum) = @_;

    #If these arguements have been seen and calculated already return the results from the hash
    return $AllowedBondTypeHash->{"$CurrentBond,$PrimarySumOfBondNum,$SecondarySumOfBondNum"} if exists $AllowedBondTypeHash->{"$CurrentBond,$PrimarySumOfBondNum,$SecondarySumOfBondNum"};
		
	  my @PossibleBondTypeChanges;

    #Calculate missing charge for both atoms given the current sum of bond numbers
		my $PrimaryMissingCharge   = $CurrentBond->getPrimary->calculateProtonationState($PrimarySumOfBondNum);
		my $SecondaryMissingCharge = $CurrentBond->getSecondary->calculateProtonationState($SecondarySumOfBondNum);
		
	  my @Range;

	  my $Greater = 0;
	  my $Lesser = 0; 

	  if($PrimaryMissingCharge >= $SecondaryMissingCharge) {

	    $Greater = $PrimaryMissingCharge; $Lesser = $SecondaryMissingCharge;

	  }

	  else {

		  $Greater = $SecondaryMissingCharge; $Lesser = $PrimaryMissingCharge;

		}				
		 
		do {
									
	    push (@Range, $Greater); $Greater--;	
											
	  } while ($Greater >= $Lesser);
			
		#Here we calculate if changing the bond type will improve the amount of missing charge for this bond, if so save it	
		foreach my $i (@Range) {

			#no reason to do all these calculations if there is no change plus 0 is always added regardless
	    next if $i == 0;

			#Cannot have a Bond with Type 0
			next if $CurrentBond->getType - $i == 0;

      #New missing charge given the change in bond type, i.e. if a Carbon atom had 3 single bonds it would have a sum of bond num of 3 and giving a missing
      #charge of -1, making one of those bonds a double would give a sum of 4 removing missing charge
			my $NewPrimaryMissingCharge   = $CurrentBond->getPrimary->calculateProtonationState($PrimarySumOfBondNum - $i);
			my $NewSecondaryMissingCharge = $CurrentBond->getSecondary->calculateProtonationState($SecondarySumOfBondNum - $i);

      #Is this a bond type we want to keep
			my $AllowableBondType = 0;

			#If both atoms have reduced missing charge then this is definitely good
			$AllowableBondType = 1 if abs($PrimaryMissingCharge) > abs($NewPrimaryMissingCharge) && abs($SecondaryMissingCharge) > abs($NewSecondaryMissingCharge);  

			#If the first atom can support its protonation state and the second atom reduces its missing charge this is also possible
			$AllowableBondType = 1 if $AllowableBondType == 0 && AllowedProtonationState($CurrentBond->getPrimary, $PrimarySumOfBondNum - $i) == 0 && abs($SecondaryMissingCharge) > abs($NewSecondaryMissingCharge);

			#Same as last case except reversed
			$AllowableBondType = 1 if $AllowableBondType == 0 && abs($PrimaryMissingCharge) > abs($NewPrimaryMissingCharge) && AllowedProtonationState($CurrentBond->getSecondary, $SecondarySumOfBondNum - $i) == 0; 

      #If any of these conditions are met save the bond type change 
			push @PossibleBondTypeChanges, [$i, $CurrentBond] if $AllowableBondType;

		}

		#No change in Bond Type is always a possibility	
	  push @PossibleBondTypeChanges, [0, $CurrentBond];
	  
	  #Save results in hash for repeat case
		$AllowedBondTypeHash->{"$CurrentBond,$PrimarySumOfBondNum,$SecondarySumOfBondNum"} = \@PossibleBondTypeChanges;
				
		return \@PossibleBondTypeChanges;

  }

}


#DEPRECATED!!!
sub NewEnsureCorrectBondTypes {
	
	my $Self = shift;
	
	my $Atoms = $Self->getAtoms;
	
	my $Bonds = $Self->getBonds;
	
	#Not really atoms with incorrect protonation states more just atoms with non zero protontation
	my @AtomsWithIncorrectProtonationStates = sort { @{$a->getBonds} <=> @{$b->getBonds} } grep { $_->getProtonationState != 0} @$Atoms;
	
	#If are no atoms with protonation states that are non zero, exit out, since its a waste of time 
	return unless @AtomsWithIncorrectProtonationStates;
	
	#Looking for potential bonds to fix, only considering bonds with atoms that can support more than one bond (by element), and contain one atom with
	#a nonzero protonation state since those are likely to change	
	my @BondsWithIncorrectProtonationStates = grep { (! CanAtomOnlySupportOneBond($_->getPrimary) && ! CanAtomOnlySupportOneBond($_->getSecondary)) &&
																									 ($_->getPrimary->getProtonationState != 0 || $_->getSecondary->getProtonationState != 0)  } @$Bonds;
	
	my %BondKeys = map {  $_ => join( " ", sort map { $_->getName } $_->ToAtomArray) } @$Bonds;
	
	my %BondsWithMissingCharge = map {  $BondKeys{$_} => 1 } @BondsWithIncorrectProtonationStates;
		
	my $NumOfAtomsWithIncorrectProtonationStates = @AtomsWithIncorrectProtonationStates;
	
	my $NumOfBonds = @$Bonds;
	
	my $BestScore = 100; my @BestGroups;	
	
	my $NearBestScore = 100; my $NearBestGroup;
	
	my $SolutionExists = 0;
	
	my $BestAtom;
	
	my %ExistingSolutions;
				
	my $CurrentAtom = shift @AtomsWithIncorrectProtonationStates;
	
	my @Groups;
	my @NewBonds = @{$CurrentAtom->getBonds};
	
	my %Visited;
	
  my $SumOfBondsTypesForAtoms = { map { $_->getName => $_->getSumofBondTypes } @$Atoms };
		
  my $BondValues =  { map {  $BondKeys{$_} => $_->getType } @$Bonds } ;
	
	push @Groups, { Atoms => $SumOfBondsTypesForAtoms, BondValues => $BondValues  };
	
	my $Count = 0;
	
	my $Done = 0;
	
	#This algorithm exhustively tries different combinations of 1/2/3 bond types for bonds that can support them to try and find
	#the combination that produces the lowest amount of formal charge for the molecule that is allowed 
	
  while ( !$Done ) {
	
	  my $GroupCount = 0;
	
	  while ( @NewBonds) {
		
		  last unless @NewBonds;
		  
		  my $Bond = shift @NewBonds;
	  	
	    next unless $BondKeys{$Bond}; #joined residues
	
	    next if exists $Visited{ $BondKeys{$Bond} };
				    
		  $Visited{ $BondKeys{$Bond} } = 1;
		    
		  push @NewBonds, (@{$Bond->getPrimary->getBonds}, @{$Bond->getSecondary->getBonds});
		    
		  next unless exists $BondsWithMissingCharge{  $BondKeys{$Bond} };
		   
	    foreach my $Group (@Groups) {
		
		    my @PossibleNext;
		     
			  my $PrimaryMissingCharge   = $Bond->getPrimary->calculateProtonationState($Group->{Atoms}{$Bond->getPrimaryName});
			  my $SecondaryMissingCharge = $Bond->getSecondary->calculateProtonationState($Group->{Atoms}{$Bond->getSecondaryName});
 
        my @Range;

        my $Greater = 0;
        my $Lesser = 0; 

        if($PrimaryMissingCharge >= $SecondaryMissingCharge) {
	
          $Greater = $PrimaryMissingCharge; $Lesser = $SecondaryMissingCharge;
         
				}
				
				else {
					
				  $Greater = $SecondaryMissingCharge; $Lesser = $PrimaryMissingCharge;
					
				}				
				
				do {
										
				  push (@Range, $Greater); $Greater--;	
												
				} while ($Greater >= $Lesser);
			

				foreach my $i (@Range) {

			    my $NewPrimaryMissingCharge   = $Bond->getPrimary->calculateProtonationState($Group->{Atoms}{$Bond->getPrimaryName} - $i);
					my $NewSecondaryMissingCharge = $Bond->getSecondary->calculateProtonationState($Group->{Atoms}{$Bond->getSecondaryName} - $i);

					if($Bond->getType - $i != 0 &&
						 ((abs($PrimaryMissingCharge) > abs($NewPrimaryMissingCharge) && abs($SecondaryMissingCharge) > abs($NewSecondaryMissingCharge)) ||
						 (AllowedProtonationState($Bond->getPrimary, $Group->{Atoms}{$Bond->getPrimaryName} - $i) == 0 && abs($SecondaryMissingCharge) > abs($NewSecondaryMissingCharge)) ||
						 (abs($PrimaryMissingCharge) > abs($NewPrimaryMissingCharge) && AllowedProtonationState($Bond->getSecondary, $Group->{Atoms}{$Bond->getSecondaryName} - $i) == 0 )) ) { 									

					  push @PossibleNext, [$i, $Bond];

				  } 
				
				}

				push @PossibleNext, [0, $Bond];
		
			  $Group->{PossibleNext} = \@PossibleNext;	
						
		  }
		
		  last;
		
		}
		
		$GroupCount = 0;
		
		my @NextGroups;
		
		foreach my $Group (@Groups) {
			
		  my $PossibleNexts = $Group->{PossibleNext};
			
			unless(defined $PossibleNexts) {
				
			  push @NextGroups, $Group; next;
				
			}
			 					    
			foreach my $PossibleNext (@$PossibleNexts) {
				
			  if($PossibleNext->[0] == 0) {
					
				  push @NextGroups, $Group; next;
				
				}
								
				my $NewGroup = { Atoms => { }, BondValues => { }  };
			    
			  while(my ($AtomName, $MissingCharge) = each %{$Group->{Atoms}} ) { $NewGroup->{Atoms}{$AtomName} = $MissingCharge }
					
			  while(my ($Bond, $BondNum) = each %{$Group->{BondValues}} ) { $NewGroup->{BondValues}{$Bond} = $BondNum }
			          
        $NewGroup->{BondValues}{ $BondKeys{ $PossibleNext->[1] } } = -$PossibleNext->[0] + $PossibleNext->[1]->getType;   
  
				$NewGroup->{Atoms}{ $PossibleNext->[1]->getPrimaryName } -= $PossibleNext->[0]; $NewGroup->{Atoms}{ $PossibleNext->[1]->getSecondaryName } -= $PossibleNext->[0]; 
  
	      my $String = "";

	      $String .= $NewGroup->{BondValues}{ $BondKeys{$_ } } foreach (@$Bonds);

	      push @NextGroups, $NewGroup unless exists $ExistingSolutions{$String};

				$ExistingSolutions{$String} = 1;
				 
			}
			
			$GroupCount++;
			
			
		}  
	
	
	  @Groups = @NextGroups;  
	
	  $GroupCount = 0;
	    
	  my $TotalCharge = 0;
		
	  foreach my $Group (@Groups) {
				
		  my $MoleculeCharge = 0;
		
		  foreach my $Atom (@$Atoms) {
			
				$MoleculeCharge += abs($Atom->calculateProtonationState($Group->{Atoms}{$Atom->getName}));

		  }
			
		
		  $TotalCharge += $MoleculeCharge;
				
		  $GroupCount++;
		
	  }
		
  	$Count++;

    my $AvgCharge = @Groups != 0 ? $TotalCharge / @Groups : 0;

    @NextGroups = ();

    if($Count != 0 && $Count % 5  == 0) {
	
	    foreach my $Group (@Groups) {
		 
			  my $MoleculeCharge = 0;

				foreach my $Atom (@$Atoms) {

					$MoleculeCharge += abs($Atom->calculateProtonationState($Group->{Atoms}{$Atom->getName}));

			  }
		
		    push @NextGroups, $Group if $MoleculeCharge < $AvgCharge;
		
	    }   
	
	    @Groups = @NextGroups;   
	
    }
	
	  $Done = 1 if !@NewBonds || $Count > 50;
		
  }
  	
		
	foreach my $Group (@Groups) {
		
	  my $Total = 0;
		
	  my $MoleculeCharge = 0;
	  	
		foreach my $Atom (@$Atoms) {
			
			$Total += AllowedProtonationState($Atom, $Group->{Atoms}{$Atom->getName});
			
			$MoleculeCharge += abs($Atom->calculateProtonationState($Group->{Atoms}{$Atom->getName}));
									
		}

		
		if (abs($Total) < abs($NearBestScore) ) {
			
			$NearBestScore = $Total; $NearBestGroup = $Group;  $BestAtom = $CurrentAtom;
						
		}
		
		next unless $Total == 0;	  		
				
		if( abs($MoleculeCharge) < abs($BestScore) ) {
				
		  $BestScore = $MoleculeCharge; 
		
		  @BestGroups = ($Group);
		
	  }
	
	  elsif( abs($MoleculeCharge) == abs($BestScore) ) {
		
		  push(@BestGroups, $Group);
		
	  }
			
	}
	
	#Keep Charge as close to the outside of the molecule as possible
	
	my $FinalGroup; my $FinalScore = -1;
	
	foreach my $Group (@BestGroups) {
		
		my $CurrentScore = 0;
		
	  foreach my $Atom (@$Atoms) {

			next if $Atom->calculateProtonationState($Group->{Atoms}{$Atom->getName}) == 0;
			
		  if($Atom->getNumOfBonds == 1) { $CurrentScore++ }
		
		  else {
			
			  foreach my $BondedAtom (@{$Atom->getBondedAtoms}) {
				
				  if($BondedAtom->getNumOfBonds == 1) { $CurrentScore += 1 / $Atom->getNumOfBonds }  
				
			  }
			
		  }

		}
		
    if($CurrentScore > $FinalScore ) {
	
	    $FinalScore = $CurrentScore; $FinalGroup = $Group;
	
    }
		
	}
	  
	if(@BestGroups == 0) {
	
	  #Possible error
	  print "No best group determined! for " . $Self->getName . "\n";
	  return;
	
  }	
	
	my $Total = 0; my $MoleculeCharge = 0;
	
	my $BestGroup;
	
	foreach my $Atom (@$Atoms) {

	  $Total += AllowedProtonationState($Atom, $FinalGroup->{Atoms}{$Atom->getName});

		$MoleculeCharge += abs($Atom->calculateProtonationState($FinalGroup->{Atoms}{$Atom->getName}));

	}
	
	foreach my $Bond (@$Bonds) {    
				
		$Bond->UpdateType(  $FinalGroup->{BondValues}{ $BondKeys{$Bond} } );
		
	}
	
}

=head2 AllowedProtonationState

Usage: AllowedProtonationState($Atom, $SumofBondTypes)

Arguments:
  $SumofBondTypes: the sum of the bond types that the atom is taking part in

Synopsis:
  A Wraper function to quickly calculate whether changing a bond number of a atom will yield an allowed protonation state

=cut


sub AllowedProtonationState {
	
	my ($Atom, $SumofBondTypes) = @_;
	
	return $Atom->AmountofChargeToAllowedProtonationState($Atom->calculateProtonationState($SumofBondTypes))
	
}


=head2 CanAtomOnlySupportOneBond

Usage: CanAtomOnlySupportOneBond($Atom);

Arguments:
  $Atom: The Atom which you want to check if its element can only support one bond, such as H, F, CL

Synopsis:
  Checks to see if $Atom can support only one bond, this is used in EnsurecorrectBondTypes, returns 1 if the atom can only support one single bond

=cut

sub CanAtomOnlySupportOneBond {
	
	my $Atom = shift;

  croak "Element " . $Atom->getElement . "does not exists please add it in the Atom.pm in BondNumHash\n" unless exists $Atom::BondNumHash{$Atom->getElement};

  #Get allowable bond numbers
  my $ElementBondNumArray = $Atom::BondNumHash{$Atom->getElement};

  #can support only one single bond
  return 1 if @$ElementBondNumArray == 1 && $ElementBondNumArray->[0] == 1;

  return 0;
	
}

=head2 NewReduceAtomicChargesToLowestPossible

Usage: $MoleculeInstance->NewReduceAtomicChargesToLowestPossible;

Arguments:

Synopsis:
  This algorithm attempts to brute force search all possible BondNumber/Protonation state combinations to see if they reduce
the total net charge of the system. The problem with NewEnsureCorrectBondTypes algorithm is it sometimes comes to suboptimal
solutions that lead to much higher molecular charge. 

=cut

sub NewReduceAtomicChargesToLowestPossible {
	
	my $Self = shift;
	
	my $Atoms = $Self->getAtoms;

	foreach my $Atom (@$Atoms) {
		
		my $TotalCharge = $Atom->getProtonationState;
		
		$TotalCharge += $_->getProtonationState foreach (@{$Atom->getBondedAtoms});
		
		my $AbsTotalCharge = abs($Atom->getProtonationState);
		
		$AbsTotalCharge += abs($_->getProtonationState) foreach (@{$Atom->getBondedAtoms});
		
		#Only consider atoms who between themselves and their bonded partners have a net formal charge
		next if $TotalCharge == 0;
		
		my $CurrentBondNum = $Atom->getSumofBondTypes;
		
		#get all possible bond numbers and protonation states for the current atom
		my $ElementBondNumArray = $Atom::BondNumHash{$Atom->getElement};
		my $ProtonationStateArray = $Atom::AllowedProtonationStates{$Atom->getElement . $Atom->getNumOfBonds};
		
		my $BestGroup; 
		
		#BestScore is the score to beat to change the current bond num / protonation state of an atom, the reason why there are 2 terms and not just 
		#the sum of all the charges in the atom bonding complex is because there are scenarios such as in nitrobenzene where the net amount of signed formal charge
		#changes but the unsigned does not when making one of these bond num / protonation state changes
		#=> It might be interesting to look into do two seperate comparisons of abs(TotalCharge) and $AbsTotalCharge = 8.9.11 JDY
		my $BestScore = abs($TotalCharge) + $AbsTotalCharge; 
				
		foreach my $BondNum (@$ElementBondNumArray) {
		
		  #If the current number of bond is greater than this possible sum of bond types then its a waste of time since it is not acheivable
		  next if $BondNum < $Atom->getNumOfBonds;
			
			foreach my $ProtonationState (@$ProtonationStateArray) {
												
        my $NewSum = $BondNum + $ProtonationState;

        my $Diff = $CurrentBondNum - $NewSum;

        my $CopiedDiff = $Diff;
								
				my %AbleToAcceptCharge;
								
				while($Diff != 0) {
					
					foreach my $Bond (@{$Atom->getBonds}) {
						
						#Cannot have a bond type of 0
						next if $Bond->getType - $Diff == 0;
						
						$Bond->MakeThisAtomPrimary($Atom);
						
						#Can this bonding partner support this change in bond type
						next unless AllowedProtonationState($Bond->getSecondary, $Bond->getSecondarySumofBondTypes - $Diff) == 0;
												
						push @{$AbleToAcceptCharge{$Diff}}, $Bond->getSecondary;
						
					}
					
					$Diff > 0 ? $Diff-- : $Diff++;
					
				}
				
				next unless keys %AbleToAcceptCharge;
				
				#Once there are possible changes in bond num / protonation state available, this part of the algorithm brute force tries to 
				#try every permutation that is allowable 				
				my @Possibles;
				my @Groups;			
				
				while (my ($key, $value) = each %AbleToAcceptCharge) {
									
					foreach my $BondedAtom (@$value) {
						
						 #Each element is the atom plus the change in bond type between itself and the atom of interest 
						 my $Element = [$BondedAtom, $key];
						
						 push @Possibles, $Element; 
						
						 #Groups keep track of how much charge is left over to successfully get the current atom into this bond num / protonation state
						 #It also keeps track of which bonded atoms are being used
						 #Group make up is an array
						 #0: is an array of strings that is in the form of AtomName.ChangeInBondType for each atom in this group 
						 #1: a Hash all of all the atoms that are in this group, this is used for fast checking to make sure an atom has not be used yet
						 #2: the amount of charge left that needs to be satisifed through change in bond types to allow for this bon num / protonation state 
						 #when this value hits 0 its done
						 push @Groups, [ [$BondedAtom->getName . "." . $key] ,{ $BondedAtom => 1 }, $CopiedDiff - $key ];
						
					 }
					
				}
											
				my $ChangedGroup = 1;
				
				while ($ChangedGroup) {
					
					$ChangedGroup = 0;
					
					my %SeenGroups;
					
					my @NewGroups;
					
					foreach my $Group (@Groups) {
						
				    #This group is finished i.e. no other change in bond types need to be made to allow for the current atmos change in bond num / protonation state				
					  if($Group->[2] == 0) { 
						
						  push @NewGroups, $Group; next;
						
						}
		
						foreach my $Possible (@Possibles) {
							
							#my ($PossibleBondedAtom, $)
							
							#Is this atom included in the group already, if so go to the next one
							next if $Group->[1]->{$Possible->[0]};
							
							#Calculate the new amount of missing charge required to be satifised 
							my $NewTotal = $Group->[2] - $Possible->[1];
							
							#This stops from the missing charge for changing signs pretty much just means we overshot our goal
							next if ($NewTotal > 0 && $Group->[2] < 0) || ($NewTotal < 0 && $Group->[2] > 0);
							
							my @NewArray = (@{$Group->[0]}, $Possible->[0]->getName . "." . $Possible->[1]);							
							
							#If this group alredy exists we don't need duplicates
							next if exists $SeenGroups{join(" ", sort @NewArray)};
							
							my @NewGroup = ( \@NewArray, { map { $_ => 1} keys %{$Group->[1]} }, $NewTotal );
							
							$NewGroup[1]->{$Possible->[0]} = 1;
							
							push @NewGroups , \@NewGroup;
							
							$SeenGroups{join(" ", sort @NewArray)} = 1;
							
							$ChangedGroup = 1;
							
						}
						
					}
					
					@Groups = @NewGroups;
					
														
				}
						
				#Scoring of groups, again looking for groups that reduce the amount of formal charge of the atom + its bonding partners, specifically
				#trying to minimize abs(SumofProtonationState) + SumofAbsoluteValuesofProtonationStates
				foreach my $Group (@Groups) {
					
					my $AbsTotal = abs($ProtonationState);
					
					my $Total = $ProtonationState;
					
					my %SeenAtoms;  

          #Split up the first element of each group to yeild the name of each atom in the group and the amount of change in bond type 
          foreach my $String (@{$Group->[0]}) {
		
	          my @spl = split(/\./, $String);
		
		        #go through all of the bonded atoms find the current bonded atom and collect its protonation state information
						foreach my $BondedAtom (@{$Atom->getBondedAtoms}) {
					
					    next if $spl[0] ne $BondedAtom->getName;
					
					    my $BondedAtomProtonationState = $BondedAtom->calculateProtonationState($BondedAtom->getSumofBondTypes - $spl[1]);
					
					    $Total += $BondedAtomProtonationState; $AbsTotal += abs($BondedAtomProtonationState);
															
					    $SeenAtoms{$BondedAtom} = 1;
										
            }

          }

					#Also consider the protonation states of atoms that are not part of the bond num / protonation state change to get the total value of the change
          foreach my $BondedAtom (@{$Atom->getBondedAtoms}) {
	
	          next if $SeenAtoms{$BondedAtom};
				
 						$AbsTotal += abs($BondedAtom->getProtonationState);
		
		        $Total += $BondedAtom->getProtonationState;
		
          }

          if(abs($Total + $AbsTotal) < abs($BestScore)) {
	 
	          $BestScore = abs($Total + $AbsTotal); $BestGroup = $Group;
	
	        }
					
				}
																		
				
			}
			
		}
		
		#If there is a solution better then the current one apply it!, update the necessary bond types stored in the best group
    if(defined $BestGroup) {
		
		  foreach my $String (@{$BestGroup->[0]}) {
			
			  my @spl = split(/\./, $String);
	      
			  foreach my $Bond (@{$Atom->getBonds}) {
		   
		      $Bond->MakeThisAtomPrimary($Atom);
					
					next if $spl[0] ne $Bond->getSecondaryName;
					
					$Bond->UpdateType($Bond->getType - $spl[1] );
		 
		    }
		   
		
	    }
	
	    last;
	
    }
 
	}
	
	
}

=head2 FindResonanceStructures

Usage: $MoleculeInstance->FindResonanceStructures;

Arguments:

Synopsis:
  Finds chemical groups that share formal charge, such as COO-, and Argine 

=cut

sub FindResonanceStructures {
	
	my $Self = shift;
	
	my $Atoms = $Self->getAtoms;
	
	foreach my $Atom (@$Atoms) {
		
		#Only want to consider atoms have have non-zero protonation states. In addition atoms that have already been found to be in resonance do not need 
		#to be reexamined 
    next if $Atom->getProtonationState == 0 || $Atom->getResonance;
		
	  my %SeenAtoms; 
	
	  #Here we are looking for atoms that are 2 bonds away and share the exact same chemical environment 
	  my @Resonance = grep { !$SeenAtoms{$_} ++ }
										grep { $Atom eq $_ || $Atom->DoesLookUpTableMatch($_->getLookUpTable) }
			              map  { @{$_->getBondedAtoms} } @{$Atom->getBondedAtoms};

    #Since the initial atom is included anything under 2 cannot be in resonance other than in a ring, this part considers if its a ring resonance scenario
	  if (@Resonance < 2) {
		
		  my $RingResonance = 0;
		
		  if($Atom->getRingAtom == 1) {
			
			   foreach my $BondedAtom (map { @{$_->getBondedAtoms} } @{$Atom->getBondedAtoms}) {
				
				   next if $Atom->Stringify ne $BondedAtom->Stringify || $Atom eq $BondedAtom;
					
					 @Resonance = ($Atom, $BondedAtom); $RingResonance = 1;
				
			   }
			
		   }
		
		  #Is not in ring resonance, formal charge is simply evenly accounted for over all bonds 
		  if(!$RingResonance) {
		
		    $Atom->setMissingCharge($Atom->getProtonationState / @{$Atom->getBonds}); next; 
		
	    }
		
	  }
	  
	  %SeenAtoms = ();
		
		#This is a clever statement, assuming that there is only one atom that is bonded to all atoms in Resonance then only the atom that has been seen
	  #  $#Resonance times should be the correct atom
		my ($SharedBondedAtom) = grep { $SeenAtoms{$_}++ == $#Resonance } 
		                         map  { @{$_->getBondedAtoms} } @Resonance;
		
		#If cannot find a shared atom something is wrong 		
		next unless defined $SharedBondedAtom;
		
		my $TotalMissingChargeInResonance = $SharedBondedAtom->getProtonationState;
		$SharedBondedAtom->setMissingCharge(0);
		
		my $DistributedCount = @Resonance;
		
		$TotalMissingChargeInResonance += $_->getProtonationState foreach (@Resonance);
	  
	  foreach my $Atom (map { @{$_->getBondedAtoms} } @Resonance) {
		
		  next if $Atom->getName eq $SharedBondedAtom->getName;
		
		  $DistributedCount++;
		
	  }
	
	
    my $DistributedCharge = $TotalMissingChargeInResonance / $DistributedCount;
	  
		foreach my $Atom (@Resonance) { 
			
		  $Atom->setMissingCharge($DistributedCharge); 
		  $Atom->setResonance(1);
		
		}
		
	} 
	
}

=head2 AbleToConnect

Usage: $MoleculeInstance->AbleToConnect;

Arguments:

Synopsis:
  Checks to see if the $MoleculeInstance is able to be connected to another molecule in the Chain object it is in

=cut

sub AbleToConnect {
	
	my $Self = shift;
	
	#A ConnectorAtom is for example N for proteins, if this is not defined then its not going to work 	
	return 0 unless defined $Self->getConnectorAtom;
			
	my $ConnectorAtom = $Self->getConnectorAtom;
	
	#If the ConnectorAtom cannot support another bond then we are done
	return 0 if $ConnectorAtom->getProtonationState > -1;	
	
	#If we don't know what to look for the make an intermolecule connection this is not going to work
	return 0 unless defined $Self->getConnectTo;
	
	my $ConnectTo = $Self->getConnectTo;
	
	return 0 if ref($ConnectTo) eq 'Atom';
	
	return 1;
	
}

=head2 UseBondingInformationFromFile

Usage: $MoleculeInstance->UseBondingInformationFromFile($Bonds);

Arguments:
  $Bonds: The Bond Hashes from File Readin. 

Synopsis:
  Attempts to use the bonding information contained from the File this Molecule was created from

=cut

sub UseBondingInformationFromFile {
	
	my ($Self, $Bonds) = @_;
	
	my $Atoms = $Self->getAtoms;
			
	my $AtomHashByName = { (map { $_->getName => $_ } @$Atoms), Indentifier => 'Name' };
	my $AtomHashByNum  = { (map { $_->getNum => $_ } @$Atoms), Indentifier => 'Num' };
	my @BondObjects;
		
	foreach my $Bond (@$Bonds) {
		
		my $BondedAtomIndentifierList;
		
		#Atoms can be indentified by either number or by name, need to see how they were stored in the Bond information hash from the input
		my $AtomIndentifierHash; 
		
		if(exists $Bond->{'BondAtomNames'}) { 
			
			$BondedAtomIndentifierList = $Bond->{'BondAtomNames'};
						
			$AtomIndentifierHash = $AtomHashByName;
		}
		
		elsif(exists $Bond->{'BondAtomNums'}) {
						
			$BondedAtomIndentifierList = $Bond->{'BondAtomNums'};
			$AtomIndentifierHash = $AtomHashByNum;			
		}
		
		croak "Cannot interpret bond hash, has neither BondAtomNames or BondAtomNums, something is wrong with the bonding info!" if ! defined $AtomIndentifierHash;
				
		my @BondedAtoms;
		my $SkipBond = 0;
		
		#Find both Atom objects that take part in this bond			
		foreach my $AtomIndentifier (@$BondedAtomIndentifierList) {
			
			#Simple case, if either the atom name or number exists in the indentifier hash, we are done with this atom 
			if(exists $AtomIndentifierHash->{$AtomIndentifier}) { 
				
				push @BondedAtoms, $AtomIndentifierHash->{$AtomIndentifier};
				
			}
			
			elsif($AtomIndentifierHash eq $AtomHashByName) {
				
				#Plus signs designate how an atom is to connect in charmm, unfortunately it is rather annoying to parse for connectivity 
				if($AtomIndentifier  =~ /^\+(\S+)/) {
					
					my $NameWithoutPlus = $1;
					
					#If removing the + sign still does not yield an atom name that exists then we have an error
					croak "Atom with + sign does not exists, these are usually connector atoms, something is wrong" unless exists $AtomIndentifierHash->{$NameWithoutPlus};
					
					my $ConnectorAtom = $AtomIndentifierHash->{$NameWithoutPlus};
					
					$Self->setConnectorAtom($ConnectorAtom); 
					
					#Whatever the atom name that the ConnectorAtom is bonded to is the atom we need to find to connect up molecules
					$Self->setConnectTo($ConnectorAtom->getName eq $BondedAtomIndentifierList->[0] ? $BondedAtomIndentifierList->[1] : $BondedAtomIndentifierList->[0]);
					
					#This bond is not going to be recorded since it connecting it to itself will mess things up this is supposed to be an intermolecular bond
					$SkipBond = 1;
					
				}
				
				elsif($AtomIndentifier =~ /^\-(\S+)/) {
										
					my $NameWithoutMinus = $1;
					
					croak "Atom with - sign does not exists, these are usually connector atoms, something is wrong" unless exists $AtomIndentifierHash->{$NameWithoutMinus};
										
					push @BondedAtoms, $AtomIndentifierHash->{$NameWithoutMinus};
					
				}
				
				elsif($Self->getIsPatch == 1) {
					
					#This atom was not actually specified in a RTF ATOM declaration but is referenced in the BOND statement, this atom is just a place holder
					#To hold the new connectivity upon applying this Topology patch
					my $DummyAtom = Atom->New( { $AtomIndentifierHash->{'Indentifier'} => $AtomIndentifier, Dummy => 1 } );
					
					$DummyAtom->Initiate;
																				
					push @$Atoms, $DummyAtom;
					push @BondedAtoms, $DummyAtom;
				}
								
			}
			
			else { 
				
				return 0 if $Parameters->getExitifNotInitiated == 0;
				
				croak "$AtomIndentifier Atom does not exist yet is referenced in the Bonds!";
								
			}
			
			
		}
		
		next if $SkipBond;
				
		croak "This Bond does not have 2 atoms in " . $Self->getName . " !!" if @BondedAtoms != 2;
		
		croak "This Bond does not have a defined Type (singe, double, etc)!" if ! exists $Bond->{'Type'};
		
		#Fix types that are not digits, such as Ar (Aromatic) in MOL2 format		
		$Bond->{'Type'} = 1 if $Bond->{'Type'} !~ /^\d+/;
		
		push @BondObjects, Bond->New(@BondedAtoms, $Bond->{'Type'});
		
	}
	
	$Self->setBonds(\@BondObjects);
	
	return 1;
	
}

=head2 PatchTopologyMoleculeWithPres

Usage: $MoleculeInstance->PatchTopologyMoleculeWithPres($PresMolecule);

Arguments:
  $PresMolecule: The Patch Molecule that one wishes to apply to $MoleculeInstance

Synopsis:
  Apply a topological patch to an molecule object

=cut

sub PatchTopologyMoleculeWithPres {
	
	my ($Self, $PresMolecule) = @_;
	
	my $Atoms = $Self->getAtoms;
	
	#Copy Patch molecule to avoid referencing nightmares
	$PresMolecule = $PresMolecule->Copy;
		
	my $Deletes = $PresMolecule->getDeletes;

  #Remove atoms from molecule object listed in the PRES delete statements	using atom names 
  foreach my $DeleteHash (@$Deletes) {
		
	  next if uc $DeleteHash->{'Type'} ne 'ATOM';
		  
    foreach (0 ..@$Atoms-1) {
	 
	    next if $Atoms->[$_]->getName ne $DeleteHash->{'Name'};
	
	    $_->UnBondAtoms foreach (@{$Atoms->[$_]->getBonds});
	
	    splice(@$Atoms, $_, 1); last
	    
    }
  
  }

  my $PresAtoms = $PresMolecule->getAtoms;

  #Go through all Atoms contained in Patch molecule and see if any of them share the name of the Molecule atoms, if they do then copy over the bonds from the
	#  Patch Atom object to the Molecule Added object.
	foreach my $PresAtom (@$PresAtoms) {
		
	  my ($SelfAtom) = grep { $PresAtom->getName eq $_->getName } @$Atoms;
	
	  #If there is no atom named the same as $PatchAtom it must be a new one and needs to be added to the Molecule object
		if(! defined $SelfAtom) { $PresAtom->setMolecule($Self); push @$Atoms, $PresAtom; next }
				
		#Dummy atoms have no info and are just place holders, if #PatchAtom is a real atom object its Type and Charge are important and need to be saved
	  if(!$PresAtom->getDummy) { 
		  $SelfAtom->setType($PresAtom->getType);
			$SelfAtom->setCharge($PresAtom->getCharge);
		}   
		
		#Bond all atoms that are bonded to the Patch Atom to the Molecule Atom to perserve connectivity
		my @DeadBonds;
		
	  foreach my $Bond (@{$PresAtom->getBonds}) {
			
		  Bond->New($SelfAtom, $Bond->getBondingPartner($PresAtom), $Bond->getType);
						
			push @DeadBonds, $Bond;
			
	  }
	
	  foreach (@DeadBonds) { $_->UnBondAtoms}
		
	}
	
	#Collect all bonds to put in internal Bond array
	$Self->GatherBonds;

  #Update total charge of molecule by summing all the atom charges
	my $NewMoleculeCharge = 0;
	foreach (@$Atoms) { $NewMoleculeCharge += $_->getCharge;}
	
  $Self->setCharge(sprintf("%.2f", $NewMoleculeCharge));	

  $Self->setName($Self->getName . "," . $PresMolecule->getName);

	
}

=head2 GatherBonds

Usage: $MoleculeInstance->GatherBonds;

Arguments:

Synopsis:
  Collects all unique Bond objects from Atom objects and puts in them the internal Bond array

=cut

sub GatherBonds {
	
	my $Self = shift;
	
	my $Atoms = $Self->getAtoms;

	my %SeenBonds;
	
	#Get all unique bonds
	my @Bonds = grep { !$SeenBonds{join(" ", sort $_->ToAtomArray)}++ } map { @{$_->getBonds} } @$Atoms; 

	$Self->setBonds(\@Bonds);

}

=head2 BondAtomsUSingIdealBondLengths

Usage: $MoleculeInstance->BondAtomsUSingIdealBondLengths;

Arguments:

Synopsis:
  Uses Bond length file to bond atoms if no bonding rules are given in Molcule formated file

=cut

sub BondAtomsUsingIdealBondLengths {

  my $Self = shift;

  my $Atoms = $Self->getAtoms;

  #Load in Ideal Bond Length Rules
  #They are in the format by element and maximum distance
  if (scalar(keys %$DistanceBondingRules) == 0) {

    croak "Cannot load Ideal Bond Lengths since variable BondLengthFilePath is not set" unless defined $Parameters->getBondLengthFilePath;

    croak "Cannot load Ideal Bond Lengths since BondLengthFilePath is set but the file does not exist" unless -e $Parameters->getBondLengthFilePath; 
  		
	  open(FILE, $Parameters->getBondLengthFilePath);
	  $DistanceBondingRules = { map { my @spl = split /\s+/, $_; join(" ", sort @spl[0,1]) => $spl[2]  } <FILE> }; 
	  close(FILE);
  	
  }

  my @Bonds;

  #Calculate distance between 2 Atoms and bond them if they are within the distance constraint range
  foreach my $i (0 .. @$Atoms-1) {
	
	  my $CoordinatesI = $Atoms->[$i]->getCartesianCoordinates;
	
	  croak("Atom does not have valid (X,Y,Z) coordinates: " . $Atoms->[$i]->DebugInfo) if ! defined $CoordinatesI->[0] || ! defined $CoordinatesI->[1] || ! defined $CoordinatesI->[2];
	
    foreach my $j ($i+1 .. @$Atoms-1) {
				
	    next if ! exists $DistanceBondingRules->{join(" ", sort ($Atoms->[$i]->getElement, $Atoms->[$j]->getElement))}; 
		    	
      my $CoordinatesJ = $Atoms->[$j]->getCartesianCoordinates;

      croak("Atom does not have valid (X,Y,Z) coordinates: " . $Atoms->[$j]->DebugInfo) if ! defined $CoordinatesJ->[0] || ! defined $CoordinatesJ->[1] || ! defined $CoordinatesJ->[2];
 		  
      my $Distance = sqrt ( ($CoordinatesI->[0] - $CoordinatesJ->[0])**2 + ($CoordinatesI->[1] - $CoordinatesJ->[1])**2 + ($CoordinatesI->[2] - $CoordinatesJ->[2])**2);

      my $DistanceConstraint = $DistanceBondingRules->{join " ", sort ($Atoms->[$i]->getElement, $Atoms->[$j]->getElement)};

      next if $Distance > $DistanceConstraint;

      push @Bonds, Bond->New($Atoms->[$i], $Atoms->[$j], 1);		  

    }

  }

  #Attempt to catch common overbondings 
  foreach my $Atom (@$Atoms) {

    if(($Atom->getElement eq "C" && $Atom->getSumofBondTypes > 4) || ($Atom->getElement eq "O" && $Atom->getSumofBondTypes  > 2) ||
       ($Atom->getElement eq "N" && $Atom->getSumofBondTypes > 4) || ($Atom->getElement eq "H" && $Atom->getSumofBondTypes  > 1)  ) {  

      foreach my $BondedAtom (@{$Atom->getBondedAtoms}) {
	
	      print $BondedAtom->getName . " " . $BondedAtom->getNum . " " . $BondedAtom->getMolecule->getName . "\n";
	
      }

	    croak($Atom->getName . " " . $Atom->getNum .  " in " . $Self->getName . " " . $Self->getNum . ": has too many bonds for its element, structure is warped") if $Parameters->getCheckElementBondingNumber;

	  }
	}


  $Self->setBonds(\@Bonds);

}

=head2 Copy

Usage: $MoleculeInstance->Copy;

Arguments:

Synopsis:
  Creates a deep copy of $MoleculeInstance, Generates new Bond objects for each bond

=cut

sub Copy {
	
	my $Self = shift;
	
	my $SelfAtoms = $Self->getAtoms;
	my $SelfBonds = $Self->getBonds;
	
	my $CopiedMolecule = Molecule->New;
  
  my %CopiedAtomNames;
  my @CopiedBonds;

  #Go through each bond contained in the current Molecule object and copy each atom involved in the bond and create a new bond with the same properties 
  foreach my $SelfBond (@$SelfBonds) {
	  
	  my @Atoms = $SelfBond->ToAtomArray;
	 
	  foreach my $Atom (@Atoms) { 
	
	    if(defined $CopiedAtomNames{$Atom->getName}) { $Atom = $CopiedAtomNames{$Atom->getName} } 
	
	    else {
		
		    $Atom = $Atom->Copy;
		    $Atom->setMolecule($CopiedMolecule);
		    $CopiedAtomNames{$Atom->getName} = $Atom;
		
	    }
	
    }

    push @CopiedBonds, Bond->New(@Atoms, $SelfBond->getType);

  }

  #If there are some atoms that are not in any bonds then they should be copied too, this only happens in PRES topology molecules
  foreach (@$SelfAtoms) {
	
	  $CopiedAtomNames{$_->getName} = $_ unless defined $CopiedAtomNames{$_->getName};
  
  }

  my @CopiedAtoms = values(%CopiedAtomNames);

	my %DoNotCopy = ( _Bonds => 1, _Atoms => 1); 
	
	while ( my ($Field, $Value) = each(%$Self)) {
		
		next if exists $DoNotCopy{$Field};
		
		$CopiedMolecule->{$Field} = $Value;
		
	}
	
	if( defined $CopiedMolecule->getConnectorAtom) {
		
		my $ConnectorAtom = $CopiedMolecule->getConnectorAtom;
	
	  foreach (@CopiedAtoms) { $CopiedMolecule->setConnectorAtom($_) if $ConnectorAtom->getName eq $_->getName }
	
  }
	
	$CopiedMolecule->setAtoms(\@CopiedAtoms);
	$CopiedMolecule->setBonds(\@CopiedBonds);
  
  return $CopiedMolecule;
	
	
}


=head2 ToStringInPdbFormat

Usage: $MoleculeInstance->ToStringInPdbFormat;

Synopsis:
  creates a string that is in Pdb file format to be written to a file

=cut

sub ToStringInPdbFormat {
	
	my ($Self, $Atoms) = @_;
	
	my $String = "";
	
	$Atoms = $Self->getAtoms unless defined $Atoms;
	
	my $CenteredAtomName = sub {
	
	  my $Name = shift;
	
	  if(length($Name) == 1)    { return "  $Name   " }
	  elsif(length($Name) == 2) { return "  $Name  "  }
	  elsif(length($Name) == 3) { return "  $Name "   }
	  elsif(length($Name) == 4) { return " $Name "    }
	  elsif(length($Name) == 5) { return " $Name"     }
	  else { croak("Atom Name to long to print to PDB file\n") }
	  
	  
	
	};
	
	foreach my $Atom (@$Atoms) {
    $String .= "ATOM" . sprintf("%+7s", $Atom->getNum) . &$CenteredAtomName($Atom->getName) . sprintf("%+4s", $Self->getName) . sprintf("%+5s", $Self->getNum) .  sprintf("%+12s", sprintf("%.3f", $Atom->getX)) . sprintf("%+8s",  sprintf("%.3f", $Atom->getY)) . sprintf("%+8s",  sprintf("%.3f", $Atom->getZ)) . sprintf("%+6s", "1.00") . sprintf("%+6s", "0.00") . sprintf("%+10s", $Self->getSegmentId) . "\n";
  }

  return $String;
	
}

sub WriteToPdb {
	
	my ($Self, $Name) = @_;
	
	$Name = $Self->getFileName if ! defined $Name;
	
	open(FILE, ">$Name.pdb");
	
	print FILE $Self->ToStringInPdbFormat;
	print FILE "\nEND\n";
	
	close(FILE);
	
	
}

sub PlaceHydrogens {
	
	my ($Self, $Atom, $Angle, $Length, $NeededHydrogens) = @_;
	
	print $Atom->getName . " " . $NeededHydrogens . "\n";
	
	my ($x, $y, $z);
		
	my $CosRadAngle = cos($Angle*(3.14/180));
	
  my $BondedAtoms = $Atom->getBondedAtoms;

  my @DiffVectors;

  my @HydrogenCoordinates;  
  my @Lengths;

  my @Magnitudes;

  foreach my $BondedAtom (@$BondedAtoms) {
	
	  my @DiffVector = map { $Atom->getCartesianCoordinates->[$_] - $BondedAtom->getCartesianCoordinates->[$_] } (0 .. 2);
		
		$x += -$DiffVector[0];
		$y += -$DiffVector[1];
		$z += -$DiffVector[2];
		
		push @Magnitudes, Magnitude(\@DiffVector);
		
	  push @DiffVectors, \@DiffVector;	
	
  }

  my @Current = ($x,$y,$z);
  
  foreach my $Coord (@Current) { $Coord /= @DiffVectors }

  foreach (0 .. $NeededHydrogens-1) {
	
	  push @HydrogenCoordinates, [$Current[0],$Current[1],$Current[2]];
	  push @Lengths, $Length;
	
  }

  my $Attempts = 0;

  my $LastScore = Score([@DiffVectors, @HydrogenCoordinates], [@Magnitudes, @Lengths], $CosRadAngle);

  my $CurrentScore = 0;

  while( $Attempts < 20000) {
	
	  my $HydrogenToMove = int(rand($NeededHydrogens));
		  	
	  my @Current = @{$HydrogenCoordinates[$HydrogenToMove]};
	
	  my $rand = int(rand(3));
	
	  my $rand2 = int(rand(1000));
	
	  my @Next;
	
		my $Move = 0.002;
		
		if($rand2 > 500 ) { $Move = -$Move }
	
	  foreach my $i (0 .. 2) {
		
		  if($i == $rand) {
			
		    push @Next, $Current[$i] + $Move;
			
		  }
		
		  else {
			
			  push @Next, $Current[$i];
			
		  }
		
	  }
	
	  $HydrogenCoordinates[$HydrogenToMove] = \@Next; 
	
	  $CurrentScore = Score([@DiffVectors, @HydrogenCoordinates], [@Magnitudes, @Lengths], $CosRadAngle);
				
	  if($CurrentScore > $LastScore) {
		
 			$HydrogenCoordinates[$HydrogenToMove] = \@Current; 		
	  }
	
	  else {
		
		  $LastScore = $CurrentScore;
		
	  }
	
	  $Attempts++;
		
 }

  my $HydrogenCount = grep { $_->getElement eq 'H' } @{$Self->getAtoms};

  my $HighNum = 0;

  foreach my $Atom (@{$Self->getAtoms}) {
	
	  if($Atom->getNum > $HighNum) { $HighNum = $Atom->getNum }
	
  }

  foreach my $Current (@HydrogenCoordinates) {

    my $Magnitude = Magnitude($Current);

    my $Scale = $Length / $Magnitude;

    my @Coords =  map { $Atom->getCartesianCoordinates->[$_] - $Current->[$_]*$Scale  } (0 .. 2);

    my $NewAtom = Atom->New({Name => "H$HydrogenCount", X => $Coords[0], Y => $Coords[1], Z => $Coords[2], Num => ++$HighNum });

    $HydrogenCount++;

    $NewAtom->Initiate;

    my $Bond = Bond->New($Atom,$NewAtom,1);

    push @{$Self->getAtoms}, $NewAtom;

    push @{$Self->getBonds}, $Bond;
  
  }

}

sub Score {
	
	my ($Vectors, $Magnitudes, $CosRadAngle) = @_;
	
	my $Score = 0;
		
	foreach my $i (0 .. @$Vectors-1) {
		
		foreach my $j ($i+1 .. @$Vectors-1) {
						
			my $Alpha = $Magnitudes->[$i]*$Magnitudes->[$j]*$CosRadAngle;
			
			my $Dot = DotProduct($Vectors->[$i], $Vectors->[$j]);
			
			#print $Alpha . " " . $Dot . " " . $Magnitudes->[$i] . " " . $Magnitudes->[$j] . "\n";
			
			$Score += ($Alpha - $Dot)**2;
			
		}
		
	}
	
	#$Score /= (@$Vectors -1);
	
	return $Score;
	
}


1;
__END__
