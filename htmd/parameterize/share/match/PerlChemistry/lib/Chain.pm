package Chain;

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
use Molecule;

require Exporter;

our @ISA = qw(BaseObject Exporter);

# Items to export into callers namespace by default. Note: do not export
# names by default without a very good reason. Use EXPORT_OK instead.
# Do not simply export all your public functions/methods/constants.

# This allows declaration	use BaseObject ':all';
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
  'TOS'  => 'Factor',
  'PIP'  => 'Factor',

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

Usage: Chain->New;

Arguments:
  $Class: should be 'Chain'

Synopsis:
  creates a new Chain object

=cut

sub New {
	
	my ($Class, $StructuralHash) = @_;
	
	
	my $Self = {
		
		_Bonds	        	=> undef,
		_FileName       	=> undef,
		_FirstPatch				=> 'NONE',
		_LastPatch		  	=> 'NONE',
		_Molecules      	=> [ ], 
		_MASS							=> 0,
		_ROTBONDS					=> 0,
		_RotatableBonds		=> [ ],
		_Indentification  => undef,
		_Name		        	=> undef,
		_Num		      		=> undef,
		_SegmentId      	=> undef,
		_Type   		    	=> undef,
		
	};
	
	bless $Self, $Class;
	
  $Self->BuildFromStructuralHash($StructuralHash) if defined $StructuralHash;
	
	return $Self;
	
}

=head2 BuildFromStructuralHash

Usage: $ChainInstance->BuildFromStructuralHash;

Arguments:
  $StructuralHash: should be 'Chain'

Synopsis:
  Loads data from StructuralHash built from loading in from file

=cut

sub BuildFromStructuralHash {
	
	my ($Self, $StructuralHash) = @_;
	
	croak("This Hash is not of a Chain") if delete $StructuralHash->{'Structure'} ne ref($Self);
	
	my $Molecules;
		
	
	while(my ($Field, $Value) = each(%$StructuralHash)) {
		
		if($Field eq 'Molecules') { $Molecules = $Value; next } 
					
		croak("$Field is not a Field of the Chain Object:\n") if ! exists $Self->{"_" . $Field};

		$Self->{"_" . $Field} = $Value;
		
	}
	
	#return if ScreenMolecule($StructuralHash) == 1;
		
	my @MoleculeObjects;
	my $TempName = "";
					
	while( my ($MoleculeNumber, $MoleculeHash) = each(%$Molecules)) {
		
		$MoleculeHash->{'Chain'} = $Self;
		
		my $MoleculeObject = Molecule->New($MoleculeHash);
		
		return 0 if $MoleculeObject->getFailed;
		
		$MoleculeObject->setFileName($Self->getFileName) if defined $Self->getFileName;
		
		$MoleculeObject->setChain($Self);
						
		$TempName .= $MoleculeObject->getName . "-";
		
		push @MoleculeObjects, $MoleculeObject;
		
	}
	
	$Self->setName(substr($TempName,0,-1)) if ! defined $Self->{_Name};
	
	$Self->setMolecules(\@MoleculeObjects);	
		
	return 1;
			
}


sub ScreenMolecule {
	
	my $ChainObject = shift;
	
	my %ElementWeight = ( H => "1.008000",  C => "12.01100", N => "14.00700", O => "15.99900", F => "18.99800", P => "30.974000", S => "32.06000", CL => "35.45000", BR => "79.90400", I => "126.90447"  );
	
	my $TotalWeight = 0;

  my $NotGood = 0;
	
  foreach my $Molecule (values %{$ChainObject->{'Molecules'}}) {
   
    my $Atoms = $Molecule->{'Atoms'};

    foreach my $Atom (@$Atoms) {
	
	   if(! exists $ElementWeight{uc($Atom->{'Element'})} ) {
		
		  $NotGood = 1; last
		
	   }
	
	   $TotalWeight += $ElementWeight{uc($Atom->{'Element'})};
		
    }
 
	
	  last;	
  }

  #print $NotGood . " " . $TotalWeight . "\n";

  return 1 if ($NotGood or $TotalWeight > 600);
	
  return 0;	
	
}


=head2 IdentifyThisChain

Usage: $ChainInstance->IdentifyThisChain;

Arguments:

Synopsis:
  Try and figure out what this chain is, protein, nucleic, ligand ..

=cut

sub IdentifyThisChain {
	
	my $Self = shift;
	
	my $Molecules = $Self->getMolecules;
	
	my $HashedIdentification;
	
	$Self->CheckForMissNamedAtoms;
	
	foreach (@$Molecules) {
		
		my $Identification = exists $ResidueDefinitions->{$_->getName} ? $ResidueDefinitions->{$_->getName} : 'Ligand';
				
    $HashedIdentification->{$Identification}++;
		
	#	print $_->getName . " " . $Identification . " ";
		
	}
	
	#print "\n";
		
	my $BestMatch = 0; my $Best = 0;
	
	while( my ($Identification, $Count) = each(%$HashedIdentification)) {
		
		if($Count > $Best) {
			
			$BestMatch = $Identification;
			$Best = $Count; 
			
		}
		
		#print $Identification . " => " . $Count . "\n";
		
	}
	
	if(($BestMatch eq 'Ion' || $BestMatch eq 'Factor') && exists $HashedIdentification->{'Ligand'}) {
		
		return 'Ligand';
		
	} 
			
	if(($BestMatch eq 'Protein' || $BestMatch eq 'Nucleic') && $Best == 1) {
		
		return 'Ligand';
		
	}
	
	return $BestMatch;

	
}

sub CalculateChainInfo {
	
	my $Self = shift;
	
	$Self->CalculateMass;
	
	
}

sub CalculateRotatableBonds {
	
	my $Self = shift;
	
	my @RotatableBonds;
	
	my %SeenBonds;
	
	foreach my $Molecule (@{$Self->getMolecules}) {
		
		foreach my $Bond (@{$Molecule->getBonds}) {
     
			#Only SingleBonds
		  next if $Bond->getType > 1;				

      my @Atoms = $Bond->ToAtomArray;      

			#Only bonds with atoms with more than 1 bond each
      next if $Atoms[0]->getNumOfBonds == 1 ||  $Atoms[1]->getNumOfBonds == 1;      			

 			#Only Bonds not in Rings
			if(scalar(@{$Atoms[0]->getRings}) > 0 && scalar(@{$Atoms[1]->getRings}) > 0) {
								
				my %TypeRingsforAtom1 = map { substr($_,-1) => 1 } @{$Atoms[0]->getRings};
				my %TypeRingsforAtom2 = map { substr($_, 1) => 1 } @{$Atoms[1]->getRings};	
				
				next if $TypeRingsforAtom1{'A'} && $TypeRingsforAtom2{'A'}; 
				
				next if $TypeRingsforAtom1{'N'} && $TypeRingsforAtom2{'N'}; 
				
				
			}			  			 

			#Not Terminal Atoms
			my $TerminalAtom = 0;
			
			foreach my $Atom (@Atoms) {
				
				my $SingleBondedAtoms = 0;
				
				foreach my $BondedAtom (@{$Atom->getBondedAtoms}) {
					
					$SingleBondedAtoms++ if $BondedAtom->getNumOfBonds && $BondedAtom->getElement =~ /H|CL|BR|F/ == 1;
					
				}
				
				if ($SingleBondedAtoms == scalar(@{$Atom->getBondedAtoms}) - 1) {
					
					#Exceptions? 
					
					$TerminalAtom = 1;
					
					#Methyl Exceptions
				  #if($Atom->getState eq "C.4" && $SingleBondedAtoms == 3) {
					
					 # my $BondingPartner = $Bond->getBondingPartner($Atom);
					
					 # $TerminalAtom = 0 if $BondingPartner->getElement !~ /^O|S|N/;
					
				 # }
					
					
					
				}
				
				last if $TerminalAtom;
	
			}
			
			next if $TerminalAtom;      
			
			my %Normal = ('C' => 4, 'N' => 3);
			
			next if scalar(@{$Atoms[0]->getRings}) == 0 && exists $Normal{$Atoms[0]->getElement} && $Atoms[0]->getNumOfBonds ne $Normal{$Atoms[0]->getElement} && 
							scalar(@{$Atoms[1]->getRings}) == 0 && exists $Normal{$Atoms[0]->getElement} && $Atoms[1]->getNumOfBonds ne $Normal{$Atoms[1]->getElement};
		
		  #Need to catch trigonal planar Nitrogens bound to C.3 atoms;
									
			my $NotRotatable = 0;
			
			foreach my $Atom (@Atoms) {
				
				my $BondedPartner = $Bond->getBondingPartner($Atom);
				
				if($Atom->getState eq "N.3") {
														
				  if(IsNitrogenPlanar($Atom) && $BondedPartner->getElement eq "C") {
					
					  $NotRotatable = 1; last
					
				  }

					
				}
				
				if($Atom->getState eq "C.3" && @{$Atom->getRings} > 0) {
					
					my $AromaticRing = 0;
					
					foreach my $Ring (@{$Atom->getRings}) { $AromaticRing = 1 if $Ring =~ /A/ }
					
					next unless $AromaticRing;
										
					#Amidine

					if($BondedPartner->getState eq "C.3") {
						
						my $TerminalNs = 0;
					
					  foreach my $BondedAtom (@{$BondedPartner->getBondedAtoms}) {
						
						  if($BondedAtom->getElement eq "N" ) {
							
							  my $SingleBondedAtoms = 0;
							
							  foreach my $NBondedAtoms (@{$BondedAtom->getBondedAtoms}) {
								
								  $SingleBondedAtoms++  if $BondedAtom->getNumOfBonds == 1;
								
							  }
							
							  $TerminalNs++ if  scalar(@{$BondedAtom->getBondedAtoms}) - 1
							
						  }
						
					  }
					
						if($TerminalNs == 2 ) {

						  $NotRotatable++; last;

						}
					
				  }
				
				  #Sulfonic acid
				
				  elsif($BondedPartner->getState eq "S.4") {
					
					   my $SingleBondedAtoms = 0;
					
					  foreach my $BondedAtom (@{$BondedPartner->getBondedAtoms}) {
					
					    $SingleBondedAtoms++ if $BondedAtom->getElement eq "O.1" && $BondedAtom->getNumOfBonds == 1;					
						
					  }
					
					  if($SingleBondedAtoms == 3) { $NotRotatable = 1; last}
					
				  }
				
				  #Nitro
				
				  elsif($BondedPartner->getElement eq "N") {
					
					  my $SingleBondedAtoms = 0;

						foreach my $BondedAtom (@{$BondedPartner->getBondedAtoms}) {

							$SingleBondedAtoms++ if $BondedAtom->getElement eq "O.1" && $BondedAtom->getNumOfBonds == 1;					

						}

						if($SingleBondedAtoms >= 1) { $NotRotatable = 1; last}
					
				  }
				
				  #Carboxylic
					
					 elsif($BondedPartner->getState eq "C.3") {

					   my $SingleBondedAtoms = 0;

						 foreach my $BondedAtom (@{$BondedPartner->getBondedAtoms}) {

						   $SingleBondedAtoms++ if $BondedAtom->getElement eq "O.1" && $BondedAtom->getNumOfBonds == 1;					

						 }

						 if($SingleBondedAtoms == 2) { $NotRotatable = 1; last}

					  }
					
				}
				
			}
			
			next if $NotRotatable;
			
			push @RotatableBonds, $Bond;
			
		}
		
	}
	
	$Self->setRotatableBonds(\@RotatableBonds);
	
  return @RotatableBonds;
   
}

sub IsNitrogenPlanar {
	
	my $NitrogenAtom = shift;
	
	my $Coordinates =  $NitrogenAtom->getCartesianCoordinates;
	
	my $CutOff = 0.60; #Just a guess!
	
  foreach my $i (0 .. @$Coordinates-1) {
	
	  my $WithinCutOff = 0;
	
	  foreach my $BondedAtom (@{$NitrogenAtom->getBondedAtoms}) {
		
		  my $BondedCoords = $BondedAtom->getCartesianCoordinates;
		
		  my $Distance = ($Coordinates->[$i] - $BondedCoords->[$i]);
		
		  $WithinCutOff++ if abs($Distance) <= $CutOff;
		
	  }
	
	  if($WithinCutOff == 3) {
		
		  return 1;
		
	  }
	
	
  }
	
	
}

sub CalculateMass {
	
	my $Self = shift;
	
	my $Mass = 0;
	
	foreach my $Molecule (@{$Self->getMolecules}) {
		
		foreach my $Atom (@{$Molecule->getAtoms}) {
	
	    $Mass += $Atom::ElementWeight{$Atom->getElement};
			
		}
		
	}
	
	$Self->setMASS($Mass);
	
	return $Mass;
	
}

sub CheckForMissNamedAtoms {
	
	my $Self = shift;
	
	foreach my $Molecule (@{$Self->getMolecules}) {
		
		if(exists $CommonMissnamedResidues->{$Molecule->getName}) {
			
			$Molecule->setName($CommonMissnamedResidues->{$Molecule->getName});
			
			if($Molecule->getName eq "CAL") {
				
				foreach my $Atom (@{$Molecule->getAtoms}) {
					
					$Atom->setName("CAL");
					
				}
				
			}
			
		}
		
		foreach my $Atom (@{$Molecule->getAtoms}) {
			
			if(exists $CommonMissnamedAtoms->{$Atom->getName} && $Atom->getNameChanged == 0) {
				
				$Atom->setName($CommonMissnamedAtoms->{$Atom->getName});
				$Atom->setNameChanged(1);
				
			}
			
		}
		
	}
	
	
}

sub GetFirstAndLastResidueByNum {
	
	my $Self = shift;
	
	
	my $First = 100000; my $Last = -1;
	my $Firsti = 0; my $Lasti = 0;
	
	foreach my $i (0 .. @{$Self->getMolecules}-1) {
				
		if($Self->getMolecule($i)->getNum < $First) { 
		
		  $First = $Self->getMolecule($i)->getNum;
		  $Firsti = $i;
			
		}
		
		if($Self->getMolecule($i)->getNum > $Last) {
			
			$Last = $Self->getMolecule($i)->getNum;
			$Lasti = $i;
		}
		
	}
	
	my $FirstResidue = $Self->getMolecule($Firsti);
	my $LastResidue = $Self->getMolecule($Lasti);
	
	return ($FirstResidue, $LastResidue);
	
	
}

sub GuessPatchesForChain {
	
	my $Self = shift;
	
	my $First = 100000; my $Last = -1;
	my $Firsti = 0; my $Lasti = 0;
	
	foreach my $i (0 .. @{$Self->getMolecules}-1) {
				
		if($Self->getMolecule($i)->getNum < $First) { 
		
		  $First = $Self->getMolecule($i)->getNum;
		  $Firsti = $i;
			
		}
		
		if($Self->getMolecule($i)->getNum > $Last) {
			
			$Last = $Self->getMolecule($i)->getNum;
			$Lasti = $i;
		}
		
	}
	
	my $FirstResidue = $Self->getMolecule($Firsti);
	my $LastResidue = $Self->getMolecule($Lasti);
	
	if($Self->getIndentification eq 'Nucleic') {
		
		#5' end
		
		my %FirstAtomNameHash = map { $_->getName => $_ } @{$FirstResidue->getAtoms};
		
		#5MET	
		if(! exists $FirstAtomNameHash{"O5'"} && ! exists $FirstAtomNameHash{'P'}) {
			
			$Self->setFirstPatch('5MET');
			
		}
		
		#5TER			
		elsif(! exists $FirstAtomNameHash{'P'}) {
			
			$Self->setFirstPatch('5TER');
			
		}
		
		#5PHO/5POM
		elsif(exists $FirstAtomNameHash{'P'}) {
			
		  my $PAtom = $FirstAtomNameHash{'P'};
			
			my $BondedAtoms = $PAtom->getBondedAtoms;
			
			my $Methyl = 0;
			
			foreach my $BondedAtom (@$BondedAtoms) {
				
				if($BondedAtom->getState eq "C.4") { $Methyl = 1; last}
				
			}
			
			if($Methyl) {
				
				$Self->setFirstPatch("5POM");
				
			}
			
			else {
			 
			  $Self->setFirstPatch('5PHO');
			
			}
			
		}
		
		#3' end 
		
		my %LastAtomNameHash = map { $_->getName => $_ } @{$LastResidue->getAtoms};
		
		#3TER
		
		if(! exists $LastAtomNameHash{"O3'"}) {
			
			$Self->setLastPatch('3TER'); return;
			
		}
		
		my $O3primeAtom = $LastAtomNameHash{"O3'"};
		
		my $BondedAtoms = $O3primeAtom->getBondedAtoms;
		
		my $PAtom;
		
		foreach my $BondedAtom (@$BondedAtoms) {
			
			if($BondedAtom->getElement eq "P") { $PAtom = $BondedAtom; last }
			
		}
		
		if(! defined $PAtom) {
			
			$Self->setLastPatch('3TER'); return;
			
		}
		
		$BondedAtoms = $PAtom->getBondedAtoms;
		
	  my $SingleBondedOAtoms;
		
  	foreach my $Atom (@$BondedAtoms) {
			
			$SingleBondedOAtoms++ if $Atom->getElement eq "O" && $Atom->getNumOfBonds == 1;
			
		}
		
	  if($SingleBondedOAtoms == 3) {
			
		  $Self->setLastPatch('3PO3');
			
	  }
	
	  else {
		
		  $Self->setLastPatch('3PHO');
		
		  #Add in 3POM
		
		}
		
	}
	
  elsif($Self->getIndentification eq 'Protein') {
			
	  if($FirstResidue->getName eq "PRO") {
		
		  $Self->setFirstPatch('PROP');
		
	  }
	
	  elsif($FirstResidue->getName eq "GLY") {
		
		  $Self->setFirstPatch('GLYP');
		
	  }
	
	  else {
		
		  $Self->setFirstPatch('NTER');
		
	  }
	
	  $Self->setLastPatch('CTER');
	
  }


}

sub StripOffWaterMolecules {
	
	my $Self = shift;
	
	my $Molecules = $Self->getMolecules;
	
	my %RemovedWaters;
		
	foreach (@$Molecules) {
		
		my $Identification = exists $ResidueDefinitions->{$_->getName} ? $ResidueDefinitions->{$_->getName} : 'Ligand';
		
    if($Identification eq 'Water') {
	
	    $RemovedWaters{$_} = 1;
	
    }
		
	}
	
	my @PrunedMolecules = grep { ! exists $RemovedWaters{$_} } @$Molecules;
	
	$Self->setMolecules(\@PrunedMolecules);
	
	
}

sub StripOffIons { 
	
	my $Self = shift;
	
	my $Molecules = $Self->getMolecules;
	
	my %RemovedIons;
		
	foreach (@$Molecules) {
		
		my $Identification = exists $ResidueDefinitions->{$_->getName} ? $ResidueDefinitions->{$_->getName} : 'Ligand';
		
    if($Identification eq 'Ion') {
	
	    $RemovedIons{$_} = 1;
	
    }
		
	}
	
	my @PrunedMolecules = grep { ! exists $RemovedIons{$_} } @$Molecules;
	
	my @Ions = grep { exists $RemovedIons{$_} } @$Molecules;
	
	$Self->setMolecules(\@PrunedMolecules);
	
	return @Ions;
	
	
}

sub StripOffCrystalizationFactors { 
	
	my $Self = shift;
	
	my $Molecules = $Self->getMolecules;
	
	my %RemovedFactors;
		
	foreach (@$Molecules) {
		
		my $Identification = exists $ResidueDefinitions->{$_->getName} ? $ResidueDefinitions->{$_->getName} : 'Ligand';
		
    if($Identification eq 'Factor') {
	
	    $RemovedFactors{$_} = 1;
	
    }
		
	}
	
	my @PrunedMolecules = grep { ! exists $RemovedFactors{$_} } @$Molecules;
	
	my @Factors = grep { exists $RemovedFactors{$_} } @$Molecules;
	
	$Self->setMolecules(\@PrunedMolecules);
	
	return @Factors;
	
	
}

sub RemoveUnknownMolecules {
	
	my $Self = shift;
	
	my $Molecules = $Self->getMolecules;
	
	my %RemovedUnknowns;
		
	foreach (@$Molecules) {
		
		my $Identification = exists $ResidueDefinitions->{$_->getName} ? $ResidueDefinitions->{$_->getName} : 'Ligand';
		
    if($Identification eq 'Ligand') {
	
	    $RemovedUnknowns{$_} = 1;
	
	    print "Molecule: " . $_->getName . " was removed from Receptor Chains since it is unknown\n";
	
    }
		
	}
	
	my @PrunedMolecules = grep { ! exists $RemovedUnknowns{$_} } @$Molecules;
	
	my @Factors = grep { exists $RemovedUnknowns{$_} } @$Molecules;
	
	$Self->setMolecules(\@PrunedMolecules);
	
	return @Factors;
	
	
}





=head2 ToStringInPdbFormat

Usage: $ChainInstance->ToStringInPdbFormat;

Synopsis:
  creates a string that is in Pdb file format to be written to a file

=cut

sub ToStringInPdbFormat {
	
	my $Self = shift;
	
	my $Molecules = $Self->getMolecules;
	
	#Make SegmentId is consistent!
	
	foreach my $Molecule (@$Molecules) {
		
		$Molecule->setSegmentId($Self->getSegmentId);
		
	}
	
	my $String = "";
	
	$String .= $_->ToStringInPdbFormat foreach( sort { $a->getNum <=> $b->getNum } @$Molecules);
	
	return $String;
	
}

=head2 Initiate

Usage: $ChainInstance->Intiate;

Synopsis:
  Setups $ChainInstance

=cut

sub Initiate {
	
	my $Self = shift;
	
	my $Molecules = $Self->getMolecules;
		
	foreach (0 .. @$Molecules-1) {
				
	  next if $_ == 0;	
		
	  next unless $Molecules->[$_]->AbleToConnect;
						
	  $Self->ConnectMolecules($_, $_-1);
					
	}
	
	#Initiate Molecules 
	
	foreach (@$Molecules) {
		
		$_->Initiate;
		
	}
	
  #Check to see if Molecule parameter ConnecterAtom is defined and connected to another residue 


}

=head2 ConnectMolecules

Usage: $ChainInstance->ConnectMolecules($MoleculeIndex1, $MoleculeIndex2);

Synopsis:
  Setups $ChainInstance

=cut

sub ConnectMolecules {
	
	my ($Self, $MoleculeIndex1, $MoleculeIndex2) = @_;

	my @Molecules = map { $Self->getMolecule($_) } ($MoleculeIndex1, $MoleculeIndex2);
	
	my $ConnectorAtom = $Molecules[0]->getConnectorAtom;
	my $ConnectTo = $Molecules[0]->getConnectTo;		
	
	my ($ConnectToAtom) = grep { $_->getName eq $ConnectTo } @{$Molecules[1]->getAtoms};
		
	confess("Cannot connect molecules $MoleculeIndex1 and $MoleculeIndex2 in this chain because $ConnectTo atom does not exist in $MoleculeIndex2") unless defined $ConnectToAtom;
		
			
	my $ConnectorBond = Bond->New($ConnectorAtom, $ConnectToAtom, 1);
		
	#$_->getLookUpTable->Update foreach (map { @{$_->getAtoms}} @Molecules);
	
	#$_->GatherBonds foreach (@Molecules);
	
	#push @{$Molecules[0]->getAtoms}, $ConnectToAtom;
	#push @{$Molecules[1]->getAtoms}, $ConnectorAtom;
			
	#push @{$_->getBonds}, $ConnectorBond foreach (@Molecules);	
	
}

=head2 AddMolecule

Usage: $ChainInstance->AddMolecule($Molecule);

Synopsis:
  Adds a Molecule Object to the interal Molecules array

=cut

sub AddMolecule {
	
	my ($Self, $Molecule) = @_;
	
	my $Molecules = $Self->getMolecules;
	
	my $Included = 0;
	
  foreach (@$Molecules) { $Included = 1 if $_ eq $Molecule }

  return if $Included;
	
	$Self->setName($Self->getName . "-" . $Molecule->getName);
	
	push @$Molecules, $Molecule;
	
}


sub PatchTopologyChainWithPres {
	
	my ($Self, $MoleculeIndex1, $MoleculeIndex2, $ConnectionMethod) = @_;
	
	my @Molecules = map { $Self->getMolecule($_) } ($MoleculeIndex1, $MoleculeIndex2);
	
	confess($ConnectionMethod->getName . "Is not a Patch Molecule") if ! $ConnectionMethod->getIsPatch;
	
	$ConnectionMethod = $ConnectionMethod->Copy;
		
	my $Deletes = $ConnectionMethod->getDeletes;
	
	foreach my $DeleteHash (@$Deletes) {

	  next if uc $DeleteHash->{'Type'} ne 'ATOM';
	
	  my ($MolNum, $AtomName) =  ( $DeleteHash->{'Name'} =~ /^(\d)(\S+)/); 
		
	  my $Atoms = $Molecules[$MolNum-1]->getAtoms;

    foreach (0 ..@$Atoms-1) {

	    next if $Atoms->[$_]->getName ne $AtomName;

	    $_->UnBondAtoms foreach (@{$Atoms->[$_]->getBonds});

	    splice(@$Atoms, $_, 1); last

    }

  }

  my $PresAtoms = $ConnectionMethod->getAtoms;

  foreach my $PresAtom (@$PresAtoms) {
	
	  my ($MolNum, $AtomName) = ( $PresAtom->getName =~ /^(\d)(\S+)/); 
		
		my ($SelfAtom) = grep { $_->getName eq $AtomName } @{$Molecules[$MolNum-1]->getAtoms};
	  
	  if (! defined $SelfAtom) {
		
		  $PresAtom->setName($AtomName);
			push @{$Molecules[$MolNum-1]->getAtoms}, $PresAtom;
			next; 
		
	  }
	
	  $SelfAtom->setCharge($PresAtom->getCharge); $SelfAtom->setType($PresAtom->getType);
	  
	  foreach my $PresBond (@{$PresAtom->getBonds}) {
 
      my $PresAtomBondPartner = $PresBond->getBondingPartner($PresAtom);

      my ($PartnerMolNum, $PartnerAtomName) = ( $PresAtomBondPartner->getName =~ /^(\d)(\S+)/); 
		  
			if($PartnerMolNum eq $MolNum) { Bond->New($SelfAtom, $PresAtomBondPartner, $PresBond->getType) }
			
			else {
				
				my $OtherMoleculesAtoms = $Molecules[ $MolNum == $MoleculeIndex1 ? $MoleculeIndex2 : $MoleculeIndex1]->getAtoms;
				
				my ($ExistingBondPartnerAtom) = grep { $_->getName eq $PartnerAtomName } @$OtherMoleculesAtoms;
				
				Bond->New($ExistingBondPartnerAtom, $SelfAtom, $PresBond->getType) if( defined $ExistingBondPartnerAtom && !$ExistingBondPartnerAtom->IsBondedTo($SelfAtom));			
				
			}
			
			$PresBond->UnBondAtoms;

    }	
	  	
	}
	
  $_->GatherBonds foreach (@Molecules);
	
}


=head2 Copy

Usage: $ChainInstance->Copy;

Synopsis:
  Makes a Deep Copy of $ChainInstance by rebuilding each Molecule Object

=cut

sub Copy { 
	
	my $Self = shift;
	
	my $Molecules = $Self->getMolecules; 
	
	my $CopiedChain = Chain->New;
	
	while ( my ($Field, $Value) = each(%$Self)) {
		
		next if $Field eq '_Molecules';
		
		$CopiedChain->{$Field} = $Value;
		
	}
	
	my @CopiedMolecules = map { $_->Copy } @$Molecules;
	
	$CopiedChain->setMolecules(\@CopiedMolecules);
		
	foreach (0 .. @CopiedMolecules-1) {
				
	  next if $_ == 0;	
		
	  next unless $CopiedMolecules[$_]->AbleToConnect;
						
	  $Self->ConnectMolecules($_, $_-1);
					
	}
	
	
	return $CopiedChain;
	
}


sub getAtoms { 
	
  my $Self = shift;

  my @Atoms;

  foreach my $Molecule (@{$Self->getMolecules}) {
	
	  push @Atoms, @{$Molecule->getAtoms};
	
  }

  return \@Atoms;
	
}

=head2 getMolecule

Usage: $ChainInstance->getMolecule($Index);

Arguments: 
  $Index: the Array Index value of the specific molecule you want, 0 would be the first

Synopsis:
  Returns a specific molecule out of the Molecule Chain;

=cut

sub getMolecule {
	
	my ($Self, $Index) = @_;
	
	confess("Index is not defined, cannot getMolecule") unless defined $Index;

	my $Molecules = $Self->getMolecules;
	
	return defined $Molecules->[$Index] ? $Molecules->[$Index] : confess("There is no molecule at $Index");
	
}

=head2 WriteToPdbFile

Usage: $ChainInstance->WriteToPdbFile;

Synopsis:
  creates a string that is in Pdb file format to be written to a file

=cut

sub WriteToPdbFile {
	
	my ($Self, $FileName) = @_;

	if(! defined $FileName) { $FileName = $Self->getFileName }
	
	if(! defined $FileName) { $FileName = "test" }
	
	my $FileString = $Self->ToStringInPdbFormat;
	
	open(FILE, ">$FileName.pdb");
	
	print FILE $FileString;
	print FILE "END\n";
	
	close(FILE);
	
	
}


1;
__END__
