#!/usr/bin/env perl

=head1 NAME

MATCH.pl

=head1 SYNOPSIS

MATCH.pl generates topology and parameter files for CHARMM force fields

=head2 EXPORT

NONE

=head1 AUTHOR

Joseph Yesselman, E<lt>jyesselm@umich.eduE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2010 by Joseph Yesselman

This library is free software; you can redistribute it and/or modifycd
it under the same terms as Perl itself, either Perl version 5.10.0 or,
at your option, any later version of Perl 5 you may have available.

=head1 FUNCTIONS

=cut

use strict;
use warnings;
use Carp;

#Get MATCH libraries 
use lib $ENV{'MATCH'} . "/lib"; 

use MATCHBaseObject;
use BaseObject ':vars';

use MATCHParameters;
use MoleculeFileHandler;

use MATCHer;
use Storable;

#Setup Default Parameters

my $DefaultParameters = MATCHParameters->New;

$DefaultParameters->Initiate;

$Parameters = $DefaultParameters;

$DefaultParameters->Initiate("../resources/" . $ARGV[0] . ".par");

my $Chains =  retrieve("StoredTopologys/" . $ARGV[0] . ".dat");

my $OutputFile = $ARGV[1];

open(FILE, $OutputFile);

my @FileContents = <FILE>;

close(FILE);

my @AtomLines = grep { $_ =~ /^\s+\S+\s+\S+\s+\S+\s+\n/ } @FileContents;

my %MessedAtoms; 

foreach my $Line (@AtomLines) {
	
	my @spl = split /\s+/, $Line;
	
	$MessedAtoms{$spl[2]} = [ ] if ! exists $MessedAtoms{$spl[2]};
	
	push @{$MessedAtoms{$spl[3]}}, ($spl[1], $spl[2]);
		
}

my $AtomTypeSubstituter = AtomTypeSubstituter->New;

$AtomTypeSubstituter->Initiate;

my $AtomTyper = AtomTyper->New;

$AtomTyper->Initiate($AtomTypeSubstituter);

my $AtomCharger = AtomCharger->New;

$AtomCharger->Initiate($AtomTypeSubstituter);

my $MoleculeParameterizer = MoleculeParameterizer->New;

$MoleculeParameterizer->Initiate($AtomTypeSubstituter);

$Parameters->{_ExitifNotInitiated} =  0;
$Parameters->{_ExitifNotTyped}     =  0;
$Parameters->{_ExitifNotCharged}   =  0;
$Parameters->{_SubstituteIncrements} = 0;

my %SavedPrintOuts;
my %TypeDiff;

my $TotalDiff = 0;

foreach my $Chain (@$Chains) {
	
	#print $Chain->getName . "\n";
		
	my @Atoms = map { @{$_->getAtoms } } @{$Chain->getMolecules};
	
	next if @Atoms < 2; #Ions
	
	my %MessedAtomHash;

  if(exists $MessedAtoms{$Chain->getName}) {
	
	  my $MessedAtomArray = $MessedAtoms{$Chain->getName};
	
	  foreach my $AtomName (@$MessedAtomArray) {
		
		  $MessedAtomHash{$AtomName} = 1;
				
	  }
	 
	
  }
	  
  my $CopiedChain = $Chain->Copy;

  foreach my $Atom (map { @{$_->getAtoms } } @{$CopiedChain->getMolecules}) {
	
	  $Atom->setCharge(0);
	
  } 

  next if $AtomCharger->ChargeAtomsInChain($CopiedChain) != 1;

  my $TotalForMolecule = 0;


  

  foreach my $MoleculeNum (0 .. @{$CopiedChain->getMolecules}-1) {
	
	  #print $Chain->getMolecule($MoleculeNum)->getName . "\n";
	
	  my %AtomChargeHash = map { $_->getName => $_->getCharge } @{$Chain->getMolecule($MoleculeNum)->getAtoms};
	  my %CopiedAtomChargeHash = map { $_->getName  => $_->getCharge } @{$CopiedChain->getMolecule($MoleculeNum)->getAtoms};

    foreach my $Atom (@{$Chain->getMolecule($MoleculeNum)->getAtoms} ) {
	
	    if($MessedAtomHash{$Atom->getName}) {
		
		    $SavedPrintOuts{$Atom->getType} = [ ] unless exists $SavedPrintOuts{$Atom->getType};
				
				my $Diff = abs($AtomChargeHash{$Atom->getName} - $CopiedAtomChargeHash{$Atom->getName});
				
				$TotalDiff += $Diff;
				
				$TypeDiff{$Atom->getType} += $Diff;
					
		    push @{$SavedPrintOuts{$Atom->getType}}, [$Atom, $Diff] if $Diff > 0.01;
		    #print $Atom->getName . " " . $Chain->getName . "  " . ($AtomChargeHash{$Atom->getName} - $CopiedAtomChargeHash{$Atom->getName})  . "\n" if abs($AtomChargeHash{$Atom->getName} - $CopiedAtomChargeHash{$Atom->getName}) > 0.01; 
		
	    }
	
    }

  }

  #print "Total for Molecule: " . $TotalForMolecule . "\n";
	
}

foreach my $Type (sort keys %SavedPrintOuts) {
	
	print $Type . " " . $TypeDiff{$Type} . "\n";
	
	foreach my $AtomInfo (@{$SavedPrintOuts{$Type}}) {
		
		print " "  . $AtomInfo->[0]->getName . " " . $AtomInfo->[0]->getMolecule->getName . " " . $AtomInfo->[1] . "\n";
		
	}
	
}

