#!/usr/bin/env perl

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

#Load MATCH Libraries 
use lib $ENV{'MATCH'} . "/lib"; 

use MATCHParameters;
use MATCHBaseObject;
use MATCHFunctions ':func';	
use Storable;

use BaseObject ':vars';

our $VERSION = '0.01';

#Exported Variables

 
#Non Exported Variables


# Preloaded methods go here.

srand(time ^ ($$ + ($$ << 15)));

my $TopologyFile;
my $PatchFile;
my $CompiledFile;
my $IsParameter = 0;

foreach my $i (0 .. @ARGV-1) {
	
	if(lc($ARGV[$i]) eq "-topologyfile") {
		
		$TopologyFile = $ARGV[$i+1]; $IsParameter = 1;
		
	}
	
	elsif(lc($ARGV[$i]) eq "-patchfile" ) {
		
		$PatchFile = $ARGV[$i+1];  $IsParameter = 1;
		
	}
	
	elsif(lc($ARGV[$i]) eq "-compiledfile") {
		
		$CompiledFile = $ARGV[$i+1]; $IsParameter = 1;
		
	}
	
	elsif(lc($ARGV[$i]) =~ /^-h/) {
		
		PrintUsage();
		exit 1;
		
	}
	
	else {
		
		if($IsParameter) { $IsParameter = 0; next }
		
		croak($ARGV[$i] . "is not recognized\n");
		
	}
	
}


my $DefaultParameters = MATCHParameters->New;

$DefaultParameters->Initiate;

$DefaultParameters->setExitifNotInitiated(0);

$Parameters = $DefaultParameters;


my $Chains;

if(defined $TopologyFile) {
	
	$Chains = SetupMoleculesFromTopology($TopologyFile, $PatchFile);
	
}

elsif(defined $CompiledFile) {
	
	$Chains = retrieve($CompiledFile);
	
}

else {
	
  print "No file specified!!\n";
  PrintUsage();	
  exit 1;
	
}

print "made it!\n";

my $Name = defined $TopologyFile ? $TopologyFile : $CompiledFile;

my @spl = split /\//, $Name;

$Name = pop @spl;

@spl = split /\./, $Name;

if(@spl > 1) {
	
  pop @spl;

}

$Name = join(".", @spl);


my @Molecules = grep { $_->getName ne "DMTT" && $_->getName ne "TMAO" } map { @{$_->getMolecules } } @$Chains;

my @Bonds = grep { $_->getMolecule->getName ne "DMTT" && $_->getMolecule->getName ne "TMAO" } map { @{$_->getBonds} } map { @{$_->getMolecules} } @$Chains;

my @Atoms =  grep { $_->getMolecule->getName ne "DMTT" && $_->getMolecule->getName ne "TMAO" } map { @{$_->getAtoms} } map { @{$_->getMolecules} } @$Chains;

my $AllIncrementsForTwoTypesBonded;

SolveAllNormalIncrements(\@Bonds, 1);

foreach my $Molecule (@Molecules) {

  my $MoleculeBonds = $Molecule->getBonds;
	my @NotSolvedMoleculeBonds = grep { ! $_->IsIncrementSolved } @$MoleculeBonds;
	
	foreach my $Bond (@NotSolvedMoleculeBonds) {
		
		next if $Bond->IsIncrementSolved;
		    
		foreach my $AtomInBond ($Bond->ToAtomArray) {
    
		  my $Type = IsAtomBondedtoTwoAtomsWithSameType($AtomInBond);
		
		  if($Type ne "0" && ($Type ne $AtomInBond->getType || sprintf("%.5f", $AtomInBond->getCharge) == 0)) { 
				 
	      $Bond->MakeThisAtomPrimary($AtomInBond);
	
				 my ($OtherBond ) = grep { $_->getSecondary ne $Bond->getSecondary } 
			 									    map  { $_->MakeThisAtomPrimary($AtomInBond); $_ } 
														grep { ! $_->IsIncrementSolved } @{$AtomInBond->getBonds};
		    
					
	
	       next unless defined $OtherBond;
	
	       next if $Bond->getSecondaryMissingCharge != $OtherBond->getSecondaryMissingCharge;	      
		
		     CalculateSymmetricalIncrement($Bond, $OtherBond);
		
				 AddBondToHash($Bond); AddBondToHash($OtherBond); 

			   SolveAllNormalIncrements($MoleculeBonds, 1);
			   
	    }
    
    }
		
	}

}

foreach my $Molecule (@Molecules) {
	
	my $MoleculeBonds = $Molecule->getBonds;
	my @NotSolvedMoleculeBonds = grep { ! $_->IsIncrementSolved } @$MoleculeBonds;
 
	my $BestScore = 10; my $BestBond = 0; my $BestIncrement = 0;		
	  
  foreach my $Bond (@NotSolvedMoleculeBonds) {
	
	  my $HashKey = GenerateHashKeyForBondIncrement($Bond);
	    
    next unless exists $AllIncrementsForTwoTypesBonded->{ $HashKey};

    foreach my $BondIncrement (map { my @spl = split(/\s+/, $_); [@spl] } keys %{$AllIncrementsForTwoTypesBonded->{$HashKey}}) {
			
		  $Bond->setIncrement([$BondIncrement->[0], $BondIncrement->[1]]);
			$Bond->ApplyIncrement(-$BondIncrement->[0], -$BondIncrement->[1] );
			   
			SolveAllNormalIncrements(\@NotSolvedMoleculeBonds);
			
			my $Score = 0;
		  
		  foreach my $NotSolvedBond (@NotSolvedMoleculeBonds) {
		  
		    next if ! $NotSolvedBond->IsIncrementSolved;
			  
			  my $BondScore = 10;
			  
			  unless( exists $AllIncrementsForTwoTypesBonded->{ GenerateHashKeyForBondIncrement($NotSolvedBond) }) {
			   
			    $BondScore = 0; next; 
			
			  }
						
			  my @BondIncrementsForThisBond =  map { my @spl = split(/\s+/, $_); [@spl] } keys %{$AllIncrementsForTwoTypesBonded->{$HashKey}};
			  
			  foreach my $BondIncrementForThisBond (@BondIncrementsForThisBond) {
			  
			    my $CurrentIncrement = $NotSolvedBond->getIncrement;
			
					my $TestScore = abs($CurrentIncrement->[0] - $BondIncrementForThisBond->[0]) + abs($CurrentIncrement->[1] - $BondIncrementForThisBond->[1]);
				
					$BondScore = $TestScore if $TestScore < $BondScore;
			
				}
				
				$Score += $BondScore;
			  
			  $_->UnSetIncrement foreach(@NotSolvedMoleculeBonds);
			  
			  if($Score < $BestScore) { 
				  
				  $BestScore = $Score; $BestBond = $Bond; $BestIncrement = $BondIncrement;
				
				}
			  
		
		  }
		
		}
			
	}
	
	
		
	if($BestScore != 10 && !$BestBond->IsIncrementSolved ) {
								
	  $BestBond->setIncrement([$BestIncrement->[0], $BestIncrement->[1]]);
		$BestBond->ApplyIncrement(-$BestIncrement->[0], -$BestIncrement->[1] );
		
		#AddBondToHash($BestBond);
		
		SolveAllNormalIncrements(\@NotSolvedMoleculeBonds, 1);
		
	}

		
}

foreach my $Molecule (@Molecules) {
		
	my $MoleculeBonds = $Molecule->getBonds;
	my @NotSolvedMoleculeBonds = grep { ! $_->IsIncrementSolved } @$MoleculeBonds;

	foreach my $Bond (@NotSolvedMoleculeBonds) {

		next if $Bond->IsIncrementSolved;
		
    next if $Bond->getPrimaryCharge != -$Bond->getSecondaryCharge;

    $Bond->setIncrement([$Bond->getPrimaryCharge, $Bond->getSecondaryCharge]);
		$Bond->ApplyIncrement(-$Bond->getPrimaryCharge, -$Bond->getSecondaryCharge );
		
		SolveAllNormalIncrements($MoleculeBonds, 1);
 	
	}
	
}


print "Increments that might require refining increment rules\n";

my @IncrementToPrintToFile;
foreach my $BondIncrement (keys %$AllIncrementsForTwoTypesBonded) {
  my @Increments =  keys %{$AllIncrementsForTwoTypesBonded->{$BondIncrement}};
  
   
	my ($Largest) = sort { @{$AllIncrementsForTwoTypesBonded->{$BondIncrement}{$b}} <=> @{$AllIncrementsForTwoTypesBonded->{$BondIncrement}{$a}} } @Increments;
	my @Types = split(/\s+/, $BondIncrement);
	my @LargestIncrement = split(/\s+/, $Largest);
	
	print "Default Increment: " . $BondIncrement . " " . $Largest . "\n" if @Increments > 1;
	
	foreach my $Incr (@Increments) {
		
		my @spl = split(/\s+/, $Incr);
		
		next if $spl[0] == $LargestIncrement[0] && $spl[1] == $LargestIncrement[1];
		
		print "  Refined Increment " . $Incr . "\n";
		
		print "  Bonds that require this increment rule:\n";
		
		foreach my $Bond ( @{$AllIncrementsForTwoTypesBonded->{$BondIncrement}{$Incr}} ) {
			
			print "    " . sprintf("%-10s", $Bond->getPrimaryName) . sprintf("%-10s", $Bond->getSecondaryName) .  sprintf("%-10s", $Bond->getMolecule->getName) . "\n"; 
			
		} 
		
		print "\n";
		
	}
	

	if($Types[0] eq $Types[1] && $LargestIncrement[0] != 0) { 

	  push(@IncrementToPrintToFile, [ @Types, "0.00000", "0.00000"]);

  }

	else {
	 	
	  push(@IncrementToPrintToFile, [ @Types, @LargestIncrement]);
	
	}
	
	
	  #print sprintf("%-8s", $Types[0]) . sprintf("%-8s", $Types[1]) . sprintf("%+9s", $LargestIncrement[0]) . sprintf("%+9s", $LargestIncrement[1]) . "\n";
}	
	
open(INCR_FILE, ">$Name.incr");
	
foreach ( sort { $a->[0] cmp $b->[0] || $a->[1] cmp $b->[1] } @IncrementToPrintToFile) {
  print INCR_FILE sprintf("%-8s", $_->[0]) . sprintf("%-8s", $_->[1]) . sprintf("%+9s", $_->[2]) . sprintf("%+9s", $_->[3]) . "\n";
}
	
close(INCR_FILE);


my $Charge = 0;

foreach my $Atom (@Atoms) {
	
	$Charge += abs($Atom->getCharge);
	
}

print "Bonds that did not have their increments extracted:\n";

foreach my $Bond (@Bonds) {
	
	next if $Bond->IsIncrementSolved;
	
	#print $Bond->getPrimaryName . " " . $Bond->getSecondaryName . " " . $Bond->getPrimaryCharge . " " . $Bond->getSecondaryCharge . " " . $Bond->getPrimaryMissingCharge . " " . $Bond->getSecondaryMissingCharge . " " .$Bond->getMolecule->getName . "\n";

}

while (my ($Types, $Key) = each %$AllIncrementsForTwoTypesBonded ) {
	
	#print $Types . "\n";
	
	while( my ($Increments, $Bonds) = each %{$AllIncrementsForTwoTypesBonded->{$Types}} ) {
		
	  #print $Types . " " . $Increments . " " . scalar(@$Bonds) . "\n";
		
		
	}
		
}

#print $Charge . "\n";


sub AddBondToHash {

	my ($Bond) = @_;
	
	my $BondHashKey = GenerateHashKeyForBondIncrement($Bond);
	
	my $IncrementKey = join(" ", map { if($_ == 0) { "0.00000" } else { $_ } } map { sprintf("%.5f", $_)} @{$Bond->getIncrement});
		
	if(! exists $AllIncrementsForTwoTypesBonded->{$BondHashKey}{$IncrementKey}) { $AllIncrementsForTwoTypesBonded->{$BondHashKey}{$IncrementKey} = [ $Bond] }
	else																																		   	{ push(@{$AllIncrementsForTwoTypesBonded->{$BondHashKey}{$IncrementKey}}, $Bond); }

}

sub SolveAllNormalIncrements {

	my ($Bonds, $AddToHash) = @_;

	my $SolvedAllPossibleBonds  = 0;
	
	while(!$SolvedAllPossibleBonds) {
		
		$SolvedAllPossibleBonds = 1;
		
		foreach (@$Bonds) {
			
			if($_->getPrimary->getNumOfNotSolvedBonds == 1 && !$_->IsIncrementSolved) { 
				
				CalculateIncrement($_); 
				AddBondToHash($_) if $AddToHash;
				$SolvedAllPossibleBonds = 0; 
			
			}
			elsif($_->getSecondary->getNumOfNotSolvedBonds == 1 && !$_->IsIncrementSolved) {
			
			  $_->SwapPrimaryAtom; 
				CalculateIncrement($_); 
				AddBondToHash($_) if $AddToHash;	
				$SolvedAllPossibleBonds = 0;
			
			}
			
		}
 	}	
}


sub CalculateIncrement {
	my $Bond = shift @_;
	my $Relationship;
	
	if($Bond->getPrimaryCharge < 0.0000000000001 and $Bond->getPrimaryCharge > -0.000000000000001) { $Bond->setPrimaryCharge(0) }
	
  if($Bond->getPrimaryMissingCharge == 0 and $Bond->getSecondaryMissingCharge == 0) { $Relationship = sub { return -$_[0] } }	

  elsif($Bond->getPrimaryMissingCharge == 0) { $Relationship = sub { return -$_[0] + $Bond->getSecondaryMissingCharge }}

  elsif($Bond->getSecondaryMissingCharge == 0) { $Relationship = sub { return -$_[0] + $Bond->getPrimaryMissingCharge }}
	

  #if($Bond->getMolecule->getName eq "SDS") {
  #  print $Bond->getPrimaryName . " " . $Bond->getPrimaryType . " " . $Bond->getSecondaryName . " " . $Bond->getSecondaryType . " " . $Bond->getPrimaryCharge . " " . $Bond->getPrimaryMissingCharge . " " . $Bond->getSecondaryMissingCharge . " " . $Bond->getMolecule->getName . " " . $Bond->getPrimaryResonance . " " . $Bond->getSecondaryResonance .  "\n";
  #}


  if(defined $Relationship) {
	 # if($Bond->getMolecule->getName eq "BAMI") {
	  #  print $Bond->getPrimaryName . " " . $Bond->getPrimaryType . " " . $Bond->getSecondaryName . " " . $Bond->getSecondaryType . " " . $Bond->getPrimaryCharge . " " . &$Relationship($Bond->getPrimaryCharge) . " " . $Bond->getMolecule->getName .  "\n";
   # }
    $Bond->setIncrement([$Bond->getPrimaryCharge, &$Relationship($Bond->getPrimaryCharge)]); 
	  $Bond->ApplyIncrement(-$Bond->getPrimaryCharge, -&$Relationship($Bond->getPrimaryCharge));     
  }
}

sub CalculateSymmetricalIncrement { 
  my ($Bond, $OtherBond) = @_; 
	my $Relationship;
	
	if($Bond->getPrimaryCharge < 0.0000000000001 and $Bond->getPrimaryCharge > -0.000000000000001) { $Bond->setPrimaryCharge(0) }

  if($Bond->getPrimaryMissingCharge == 0 and $Bond->getSecondaryMissingCharge == 0) { $Relationship = sub { return -$_[0] / 2  } }	

  elsif($Bond->getPrimaryMissingCharge == 0) { $Relationship = sub { return (-$_[0] * .5 + $Bond->getSecondaryMissingCharge)  }}

  elsif($Bond->getSecondaryMissingCharge == 0) { $Relationship = sub { return (-$_[0] * .5 + $Bond->getPrimaryMissingCharge)  }}

  if(defined $Relationship) {
    #print $Bond->getPrimaryName . " " . $Bond->getPrimaryType . " " . $Bond->getSecondaryName . " " . $Bond->getSecondaryType . " " . ($Bond->getPrimaryCharge / 2) . " " . (&$Relationship($Bond->getPrimaryCharge)  ) . " " . $Bond->getMolecule->getName . "\n";
    $Bond->setIncrement([$Bond->getPrimaryCharge * .5, &$Relationship($Bond->getPrimaryCharge)]);
	  $Bond->ApplyIncrement(0, -&$Relationship($Bond->getPrimaryCharge)); 
	
	  #print $OtherBond->getPrimaryName . " " . $OtherBond->getPrimaryType . " " . $OtherBond->getSecondaryName . " " . $OtherBond->getSecondaryType . " " . ($OtherBond->getPrimaryCharge/2) . " " . &$Relationship($OtherBond->getPrimaryCharge) . " " . $OtherBond->getMolecule->getName . "\n";
	  $OtherBond->setIncrement([$OtherBond->getPrimaryCharge * .5, &$Relationship($OtherBond->getPrimaryCharge)]);
		$OtherBond->ApplyIncrement(-$OtherBond->getPrimaryCharge, -&$Relationship($OtherBond->getPrimaryCharge)); 
	}	
	
	else {
		print "made it!\n";
	}
}

sub GenerateHashKeyForBondIncrement  { 
	
  my @SortedTypes = sort map { $_->getType} $_[0]->ToAtomArray;
  
  if($SortedTypes[0] ne $_[0]->getPrimaryType ) { $_[0]->SwapPrimaryAtom }

  return join(" ", @SortedTypes);
}

sub IsAtomBondedtoTwoAtomsWithSameType {
	
	my $Atom = shift;
	
	return 0 if $Atom->getNumOfNotSolvedBonds != 2;
	
	my @NotSolvedBonds = map { if($Atom eq $_->getPrimary ) { $_->getSecondary } 
														 else													{ $_->getPrimary   } }
											 grep { ! $_->IsIncrementSolved } @{$Atom->getBonds};

  #print $Atom->getName . " " . $Atom->getMolecule->getName . "\n";

 if($NotSolvedBonds[0]->getType eq $NotSolvedBonds[1]->getType) { return $NotSolvedBonds[0]->getType}
  else												                                  { return 0 }
	
}

sub PrintUsage {
	
	print "Usage:\n";
	print "BuildBondIncrement -TopologyFile TopologyFilePath [-PatchFile PatchFilePath]\n";
	print "BuildBondIncrement -CompiledFile CompiledFile\n";
	
}



