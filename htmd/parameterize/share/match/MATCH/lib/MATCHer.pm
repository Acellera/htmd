package MATCHer;
use Math::Round qw(round);

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

use AtomTyper;
use AtomTypeSubstituter;
use AtomCharger;
use MoleculeParameterizer;

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
	  _AtomCharger					 => undef,
	  _AtomTyper						 => undef,
	  _AtomTypeSubstituter   => undef,
	  _ForceFieldName        => undef,
	  _ImproperTypes				 => undef,
	  _MoleculeParameterizer => undef,
  };
	
	bless $Self, $Class;	
	return $Self;
	
}

sub Initiate {
	
	my $Self = shift;
	
	my $AtomTyper = AtomTyper->New;
	$AtomTyper->Initiate;

	my $AtomTypeSubstituter = AtomTypeSubstituter->New;
	$AtomTypeSubstituter->Initiate;

	my $AtomCharger = AtomCharger->New;
	$AtomCharger->Initiate($AtomTypeSubstituter);

	my $MoleculeParameterizer = MoleculeParameterizer->New;
	$MoleculeParameterizer->Initiate($AtomTypeSubstituter);
	
	$Self->SetUpImproperList();
	
	$Self->setAtomTyper($AtomTyper);
	$Self->setAtomTypeSubstituter($AtomTypeSubstituter);
	$Self->setAtomCharger($AtomCharger);
	$Self->setMoleculeParameterizer($MoleculeParameterizer);
	
}

sub SetUpImproperList {
	
	my $Self = shift;
	
	return if ! defined $Parameters->getImproperFilePath;
	
	my @ImproperTypes;
	
	open(FILE, $Parameters->getImproperFilePath);
	
	my @ImprContents = <FILE>;
	
	close(FILE);

	foreach my $Line (@ImprContents) {

	  my @ArrayedLine = split (/\s+/, $Line);

	  next if @ArrayedLine != 3;

	  my $Type = Type->New(@ArrayedLine);

	  $Type->Initiate;

	  push @ImproperTypes, $Type;

	}

	@ImproperTypes = sort { length($b->getString) <=> length($a->getString) } @ImproperTypes;
	
	$Self->setImproperTypes(\@ImproperTypes);
	
}

sub BuildTopologyFileForMolecule {
	
	my ($Self, $Molecule, $PreviousTopologyName) = @_; 
	
	my $ReturnValue = $Self->getAtomTyper->TypeAtomsInMolecule($Molecule);
	
	return $ReturnValue if $ReturnValue != 1;
	
	$ReturnValue = $Self->getAtomCharger->ChargeAtomsInMolecule($Molecule);

	return $ReturnValue if $ReturnValue != 1;

  #NOT necessary!!
	$Self->WriteMoleculeToTopologyFile($Molecule, $PreviousTopologyName);
	
	return 1;
	
}

sub BuildParameterFileForMolecule {
	
  my ($Self, $Molecule, $PreviousParameterName) = @_; 
	
	$Self->getMoleculeParameterizer->NewParameterizeMolecule($Molecule);
	
}

sub BuildBothTopologyAndParameterFileForMolecule {
	
	my ($Self, $Molecule, $PreviousTopologyName) = @_;
	
	my $ReturnValue = $Self->BuildTopologyFileForMolecule($Molecule, $PreviousTopologyName);
	
	return $ReturnValue if $ReturnValue != 1;
	
	$Self->BuildParameterFileForMolecule($Molecule);
	
	return 1;
	
}

=head2 WriteToTopologyFile

Usage: $MATCHer->WriteMoleculeToTopologyFile;

Synopsis:
  Creates a RTF File to be used to load in molecule into CHARMM

=cut

sub WriteMoleculeToTopologyFile {
	
	my ($Self, $Molecule, $TopologyFileName) = @_;
	
	
	$TopologyFileName = $Molecule->getFileName unless defined $TopologyFileName;
	
	$TopologyFileName = "UNK" unless defined $TopologyFileName;
	
	$Molecule->setName("UNK") if ! defined $Molecule->getName;
	
	my $FileName = defined $Molecule->getFileName ? $Molecule->getFileName : "UNK";
			
	my $Atoms = $Molecule->getAtoms;
	my $Bonds = $Molecule->getBonds;
	
	my $AppendingTopologyFilePath = $Parameters->getAppendingTopologyFilePath;
	
	#print $AppendingTopologyFilePath . " " . $TopologyFileName . "\n";
	#my $AddHydrogenBondToplogyInformation = $Parameters->getAddHydrogenBondToplogyInformation;
	
	my %SeenTypes;
	my @UniqueTypes = grep { ! $SeenTypes{$_}++ } map { $_->getType } @$Atoms;
	
	#Save the Atomic Element and Element Weight for a given type using the information from the atoms that contain the type of interest
	my %ElementInformationForTypes = map { $_->getType => [ $_->getElementWeight, $_->getElement] } @$Atoms;
				
	open(FILE, $AppendingTopologyFilePath);
	my @AppendingTopologyFileContents = <FILE>;
	close(FILE);
	
	#Each Type in a CHARMM Topology has a MASS number types that are to be added need a MASS number for CHARMM to use them
	my %MassNumberForType =  map  { my @spl = split(/\s+/, $_); $spl[2] => $spl[1] } 
										       grep { $_ =~ m/^MASS/ } @AppendingTopologyFileContents;
	
	my @SortedMassNumbers = sort { $b <=> $a } values(%MassNumberForType);
  my $LastMassNumber = $SortedMassNumbers[0];

  my %UniqueTypesMassValues = map { if(exists $MassNumberForType{ShortenType($_)}) {$_ => $MassNumberForType{ShortenType($_)}}
 																	  else																			   	 {$_ => -1} } @UniqueTypes;

	#Split Topology file to be used on the last MASS statement so new ones can be added easily
	my ($AllMassStatementsInAppendingTopologyFile, $LastMassLine, $EverythingElseInAppendingTopologyFile) = split(/(MASS\s+$LastMassNumber\s*\w+\s+\d+\.\d+\s+\w.*\n)/, join("", @AppendingTopologyFileContents));

  open(FILE, ">top_$TopologyFileName.rtf");

  print FILE $AllMassStatementsInAppendingTopologyFile;
  print FILE $LastMassLine;

  #Add new MASS statements to new Topology file
  foreach (sort keys %UniqueTypesMassValues) { 

	  next if $UniqueTypesMassValues{$_} != -1;
	
	  $UniqueTypesMassValues{$_} = ++$LastMassNumber; 
	  print FILE "MASS" . sprintf("%+6s", $UniqueTypesMassValues{$_}) . " " . sprintf("%-6s" , ShortenType($_)) . join(" ", @{$ElementInformationForTypes{$_}}) . "\n"; 

  }

  print FILE $EverythingElseInAppendingTopologyFile;

  close(FILE);

	open(FILE, ">" . $FileName.  ".rtf");
  print FILE "* Charmm RTF built by MATCH \n*\n  22     0\n";

	while(my ($UniqueType, $MassNumber) = each(%UniqueTypesMassValues)) {
	  print FILE "MASS" . sprintf("%+6s", $UniqueTypesMassValues{$UniqueType}) . " " . sprintf("%-6s" , ShortenType($UniqueType)) . join(" ", @{$ElementInformationForTypes{$UniqueType}}) . "\n"; 
  }
  
  print FILE "\nAUTO ANGLES DIHE\n\n";  
#  print FILE "RESI  " . sprintf("%-6s", $Molecule->getName) . sprintf("%.06f",round($Molecule->getCharge)) . "\n";
  print FILE "RESI  MOL " . sprintf("%.06f",round($Molecule->getCharge)) . "\n"; # MJH -hardcode name
  print FILE "GROUP\n";
  
  foreach my $Atom (@$Atoms) { 
	  print FILE "ATOM " . sprintf("%-4s", $Atom->getName) .  " " . sprintf("%-4s", ShortenType($Atom->getType)) . " " . sprintf("%10s", sprintf("%.06f", $Atom->getCharge)) . "\n";
	}

  foreach my $Bond (@$Bonds) {
	  print FILE "BOND " . sprintf("%-5s", $Bond->getPrimaryName) . sprintf("%-5s", $Bond->getSecondaryName) . "\n"
  }

  my @Impropers = $Self->FindImproperAtoms($Molecule);

  foreach my $Improper (@Impropers) {
	
	  print FILE "IMPR ";
	
	  $Molecule->setImpropers(\@Impropers);
	
	  foreach my $Atom (@$Improper) {
		
		  print FILE sprintf("%-5s", $Atom->getName);
		
	  } 
	
	  print FILE "\n";
	
  } 

  print FILE $Self->getHydrogenBondTopologyInformation($Molecule) if $Parameters->getAddHydrogenBondToplogyInformation;

  print FILE "PATCH FIRST NONE LAST NONE\n\nEND\n";
  close(FILE);
	
	
}

sub getHydrogenBondTopologyInformation {
	
  my ($Self, $Molecule) = @_;	

	my $Atoms = $Molecule->getAtoms;
	
	my $StringToBeAdded;
	
	foreach my $Atom (@$Atoms) {
		
		if($Atom->getElement eq "N") {
			
			my $BondedAtoms = $Atom->getBondedAtoms;
			
			my @HydrogenAtomBondToAtom = grep { $_->getElement eq "H" } @$BondedAtoms;
			
			if(@$BondedAtoms == 2) { $StringToBeAdded .= "ACCE " . $Atom->getName . "\n"; }
			
			elsif(@$BondedAtoms > 2 && $#HydrogenAtomBondToAtom > -1 ) {
				
				foreach my $HydrogenAtom (@HydrogenAtomBondToAtom) { $StringToBeAdded .= "DONO " . $HydrogenAtom->getName . " " . $Atom->getName . "\n" }
			
			}
		}
		
		elsif($Atom->getElement eq "O") {
			
			my $BondedAtoms = $Atom->getBondedAtoms;
			
			if(@$BondedAtoms == 2)    { $StringToBeAdded .= "ACCE " . $Atom->getName . "\n" }
			
			elsif(@$BondedAtoms == 1) { $StringToBeAdded .= "ACCE " . $Atom->getName . " " . $BondedAtoms->[0]->getName . "\n" }
		  
		}
	}
	
	return $StringToBeAdded;
	
}



sub FindImproperAtoms {
	
	my ($Self,$Molecule) = @_;
	
	my $Atoms = $Molecule->getAtoms;
	
	my @GuessedImpropers;
	
	return unless $Self->getImproperTypes;
	
	my $ImproperTypes = $Self->getImproperTypes;
	
	foreach my $Atom (@$Atoms) {
		
		my $Matched = 0;
		my $ImproperNum;
			
	  foreach my $ImproperType (@$ImproperTypes) {
		
      next if $ImproperType->getName ne $Atom->getType;
		
		  next unless $ImproperType->getLookUpTable->AreAllNodesSharedBetween([$Atom->getLookUpTable]);
    
	    $Matched = 1; $ImproperNum = $ImproperType->getBondNumber; last;
		
	  }
	
	  next if $Matched == 0;
	 
	  next if $ImproperNum == 0;
	
	  #print $Atom->getName . "\n";
	
	  my $BondedAtoms = $Atom->getBondedAtoms;
	
	  my %InAtomRing;
	  
	  foreach my $Atom (@$BondedAtoms) {
		
		  $InAtomRing{$Atom->getName} = 0;
		
	  }
	
	  if($Atom->getRingAtom) {
		
		  my $Rings = $Molecule->getRings;
				
		  foreach my $Ring (@$Rings) {
			
			  my $IsAtomInTheRing = 0;
			
			  foreach my $RingAtom (@$Ring) {
				
				  if($RingAtom eq $Atom) {
					
					  $IsAtomInTheRing = 1;
					  last; 
					
				  }
				
			  }
			
			  if($IsAtomInTheRing) {
				
				  foreach my $RingAtom (@$Ring) {
					
					  $InAtomRing{$RingAtom->getName} = 1;
					
				  }
				
			  }
			
		  }
			
	  }
	
	  my @SortedAtoms = sort { $InAtomRing{$b->getName} <=> $InAtomRing{$a->getName} || $b->getSumofBondTypes <=> $a->getSumofBondTypes } @$BondedAtoms;
		
		unshift (@SortedAtoms, $Atom);
		
	  my $Improper = \@SortedAtoms;
	
	  push @GuessedImpropers, $Improper;
	
	  if($ImproperNum == 2) {
		
		  my $NewImproper = [$Improper->[0], $Improper->[2], $Improper->[1], $Improper->[3]];
		
		  push @GuessedImpropers, $NewImproper;
		
	  }
	
	  elsif($ImproperNum == 3) {
		
		  my $NewImproper = [$Improper->[0], $Improper->[1], $Improper->[3], $Improper->[2]];

			push @GuessedImpropers, $NewImproper;
		
	  }
	
	  #print join(" ", ($Atom->getName, map { $_->getName} @SortedAtoms )) . "\n";
	
	}
	
	return @GuessedImpropers;
	
}
