#!/usr/bin/env perl

use 5.010000;
use strict;
use warnings;
use Carp;
use Cwd qw(abs_path);

#Load MATCH Libraries 
use lib $ENV{'MATCH'} . "/lib"; 

use MATCHBaseObject;
use BaseObject ':vars';

use MATCHParameters;
use MATCHFunctions ':func';
use AtomCharger;
use AtomTyper;
use AtomTypeSubstituter;
use Storable;

our $VERSION = '0.01';

#Exported Variables


#Non Exported Variables


# Preloaded methods go here.

#Initate Random Number generator
srand(time ^ ($$ + ($$ << 15)));

#Deal With Initial Parameter Setup, This will be handled in main .pl script
my $DefaultParameters = MATCHParameters->New;

$DefaultParameters->Initiate("../resources/" . $ARGV[0] . ".par");

$Parameters = $DefaultParameters;

#Set Variables

$Parameters->{_ExitifNotInitiated}   = 0;
$Parameters->{_ExitifNotTyped}       = 0;
$Parameters->{_ExitifNotCharged}     = 0;
$Parameters->{_SubstituteIncrements} = 0;

my $AtomTypeSubstituter = AtomTypeSubstituter->New;

$AtomTypeSubstituter->Initiate;

my $AtomTyper = AtomTyper->New;

$AtomTyper->Initiate($AtomTypeSubstituter);

my $AtomCharger = AtomCharger->New;

$AtomCharger->Initiate($AtomTypeSubstituter);

my $Chains =  retrieve("StoredTopologys/" . $ARGV[0] . ".dat");

my $PrimaryType = $ARGV[1];
my $SecondaryType = $ARGV[2];

my @Atoms = map { @{$_->getAtoms} } map { @{$_->getMolecules} } @$Chains;

my @TypeAtoms = grep { $_->getType eq $PrimaryType } @Atoms;

my @Final;

##### get topologyfile 

my %ResiLineNumbers;

open(FILE, $ENV{"MATCH"} . "/resources/top_all36_cgenff_new/top_all36_cgenff.rtf");

my @TopFileContents = <FILE>;

close(FILE);

my $Count = -1;

foreach my $Line (@TopFileContents) {
	
	$Count++;
		
  next if uc($Line) !~ /RESI\s+(\S+)\s+(-?\d+\.\d+)/;

  $ResiLineNumbers{$1} = $Count;
	
}

#exit 1;

###


foreach my $Atom (@TypeAtoms) {
	
	my $Bonded = $Atom->getBondedAtoms;
	
	foreach my $BondedAtom (@$Bonded) {
		
		if($BondedAtom->getType eq $SecondaryType) {
			
			push @Final, [$Atom, $BondedAtom]; last
			
		}
		
	}
	
}

my $TotalDiff = 0;

my %Seen;
	
my $MoleculeCount = 0;

foreach my $Atoms (@Final) {
	
	my %Seen;
	
  my @Bonds = @{$Atoms->[0]->getMolecule->getBonds};

  my $CopiedMolecule = $Atoms->[0]->getMolecule->Copy;
 
  foreach my $Atom (@{$Atoms->[0]->getMolecule->getAtoms}) {
	
	  $Atom->setCharge(0);
	
  }

  foreach my $Bond (@Bonds) {

    my $Print = 0; 

    if($Bond->getPrimaryType eq $PrimaryType && $Bond->getSecondaryType eq $SecondaryType) {
	
	    $Print = 1;
	
    }

    elsif($Bond->getPrimaryType eq $SecondaryType && $Bond->getSecondaryType eq $PrimaryType) {
	
	    $Print = 1; $Bond->SwapPrimaryAtom;
	
    }

    $AtomCharger->ChargeBond($Bond, $Print);	

  }

  my %Diffs;

  my $PrintOut = "";

  my %CopiedAtomChargeHash = map { $_->getName => $_->getCharge } @{$CopiedMolecule->getAtoms};

  foreach my $Atom (@{$Atoms->[0]->getMolecule->getAtoms}) {
	
	 if(abs($Atom->getCharge - $CopiedAtomChargeHash{$Atom->getName}) > 0.01) {
	
	    $PrintOut .= $Atom->getName . " " . $Atom->getType . " " . ($Atom->getCharge - $CopiedAtomChargeHash{$Atom->getName}) . "\n";
	
	    $Diffs{$Atom->getType} = 1;
	
	 }
	
	  $TotalDiff += abs($Atom->getCharge - $CopiedAtomChargeHash{$Atom->getName});
	
  }

  #print $Atoms->[0]->getMolecule->getName. "\n";

  if(!$Seen{$Atoms->[0]->getMolecule->getName} && exists $Diffs{$PrimaryType} && exists $Diffs{$SecondaryType} ) {
		 
	  $MoleculeCount++;
	
	  $Seen{$Atoms->[0]->getMolecule->getName} = 1;
	
	  my $Start = $ResiLineNumbers{$Atoms->[0]->getMolecule->getName};
	
	  next unless defined $Start;
	  
	  my $Done = 0;
	
	  print $TopFileContents[$Start];
	
	  while(!$Done) {
		
	    if((uc($TopFileContents[$Start]) =~ /^ATOM/ || uc($TopFileContents[$Start]) =~ /^GROUP/ || $TopFileContents[$Start] =~ /^ /) && $TopFileContents[$Start] =~ /\!/) {
		
		    print $TopFileContents[$Start];
		
	    }
	
	    elsif(uc($TopFileContents[$Start]) =~ /$PrimaryType/ || uc($TopFileContents[$Start]) =~ /$SecondaryType/) {
		
		    print $TopFileContents[$Start];
		
  	  }
	
	    $Start++;
	
	   if(uc($TopFileContents[$Start]) =~ /RESI/ || $Start > @TopFileContents) { $Done = 1 }
	    
	
	  }
	
	  print $PrintOut . "\n";
	  
	
  }
	
}

print "Total Diff " . $TotalDiff . " " . $MoleculeCount . "\n";
