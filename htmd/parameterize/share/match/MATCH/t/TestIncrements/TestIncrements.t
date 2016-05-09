#!/usr/bin/env perl

use strict;
use warnings;

use Test::More tests => 1;
use lib $ENV{'MATCH'} . "/lib"; 
use Carp;
use Storable;

use MATCHBaseObject;
use BaseObject ':vars';

use MATCHParameters;
use MATCHFunctions ':all';
use MoleculeFileHandler;

use AtomTyper;
use AtomCharger;
use AtomTypeSubstituter;
use MoleculeParameterizer;

#Deal With Initial Parameter Setup, This will be handled in main .pl script

my $DefaultParameters = MATCHParameters->New;

$DefaultParameters->Initiate("../../resources/" . $ARGV[0] . ".par");
#$DefaultParameters->Initiate;

$Parameters = $DefaultParameters;

$Parameters->{_ExitifNotInitiated} =  0;
$Parameters->{_ExitifNotTyped}     =  0;
$Parameters->{_ExitifNotCharged}   =  0;

my $Chains =  retrieve("../StoredTopologys/" . $ARGV[1] . ".dat");

my $AtomTypeSubstituter = AtomTypeSubstituter->New;

$AtomTypeSubstituter->Initiate;

my $AtomTyper = AtomTyper->New;

$AtomTyper->Initiate($AtomTypeSubstituter);

my $AtomCharger = AtomCharger->New;

$AtomCharger->Initiate($AtomTypeSubstituter);

my $MoleculeParameterizer = MoleculeParameterizer->New;

$MoleculeParameterizer->Initiate($AtomTypeSubstituter);

$Verbosity = 2;

my $TotalChargeDiff = 0;
my $Count = 0;
my $MoleculeCount = 0;
my $TotalCount = 0;
my $TotalError = 0;

my ($Name1) = ($ARGV[0] =~ /\d+\_(\S+)$/);
my ($Name2) = ($ARGV[1] =~ /\d+\_(\S+)$/);

my ($xsum, $ysum, $xysum, $xsqsum, $ysqsum);

open(FILE, ">FF_" . $Name1 . "_" . $Name2. ".dat");
 
my %NotAllowed = map { $_ => 1 } qw(CRBZ MEOI CYIN CPDE EIND); 
  
foreach my $Chain (@$Chains) {

  #next if $Chain->getName eq "HEME";
	
	next if $NotAllowed{$Chain->getName};
	
	print $Chain->getName . "\n";
		
	my @Atoms = map { @{$_->getAtoms } } @{$Chain->getMolecules};
	
	next if @Atoms < 2; #Ions
	
	#next if $AtomTyper->TypeAtomsInChain($Chain) != 1;
  
  my $CopiedChain = $Chain->Copy;

  foreach my $Atom (map { @{$_->getAtoms } } @{$CopiedChain->getMolecules}) {
	
	  $Atom->setCharge(0);
	
  } 

  next if $AtomCharger->ChargeAtomsInChain($CopiedChain) != 1;

  $MoleculeCount++;

  my $TotalForMolecule = 0;

  my $String;

  foreach my $MoleculeNum (0 .. @{$CopiedChain->getMolecules}-1) {
	
	  #print $Chain->getMolecule($MoleculeNum)->getName . "\n";
	
	  my %AtomChargeHash = map { $_->getName => $_->getCharge } @{$Chain->getMolecule($MoleculeNum)->getAtoms};
	  my %CopiedAtomChargeHash = map { $_->getName  => $_->getCharge } @{$CopiedChain->getMolecule($MoleculeNum)->getAtoms};

    foreach my $Atom (@{$Chain->getMolecule($MoleculeNum)->getAtoms} ) {
	
	    #print FILE  $AtomChargeHash{$Atom->getName} .  " " . $CopiedAtomChargeHash{$Atom->getName}  . "\n";
	    	
	    $TotalError += abs($AtomChargeHash{$Atom->getName} - $CopiedAtomChargeHash{$Atom->getName});
	
	    if(abs($AtomChargeHash{$Atom->getName} - $CopiedAtomChargeHash{$Atom->getName} ) > 0.01 ) {
	      $String .= $Atom->getName . " " . $Atom->getType . " " . $AtomChargeHash{$Atom->getName} .  " " . $CopiedAtomChargeHash{$Atom->getName} . "\n";
      }
	
	    if($AtomChargeHash{$Atom->getName} != 0 && $CopiedAtomChargeHash{$Atom->getName} != 0) {
	
	      $TotalChargeDiff += abs($AtomChargeHash{$Atom->getName} - $CopiedAtomChargeHash{$Atom->getName} ) / abs($AtomChargeHash{$Atom->getName});
	 
	      print FILE  $Atom->getName . " " . $Atom->getMolecule->getName . " " .  $Atom->getType . " " . $AtomChargeHash{$Atom->getName} .  " " . $CopiedAtomChargeHash{$Atom->getName}  . "\n";
		   		
      }

      elsif($AtomChargeHash{$Atom->getName} == 0 && $CopiedAtomChargeHash{$Atom->getName} != 0) {	      
	   	
	      $Count--;
	
      }

      elsif($AtomChargeHash{$Atom->getName} == 0 && $CopiedAtomChargeHash{$Atom->getName} == 0) {
	
	      print FILE  $Atom->getName . " " . $Atom->getMolecule->getName . " " .  $Atom->getType . " " . $AtomChargeHash{$Atom->getName} .  " " . $CopiedAtomChargeHash{$Atom->getName}  . "\n";
		
      }

      $xsum += $AtomChargeHash{$Atom->getName}; $ysum += $CopiedAtomChargeHash{$Atom->getName}; $xysum += $AtomChargeHash{$Atom->getName}*$CopiedAtomChargeHash{$Atom->getName};
		  $xsqsum += ($AtomChargeHash{$Atom->getName}) **2;  $ysqsum += ($CopiedAtomChargeHash{$Atom->getName}) **2;
	
	    $TotalForMolecule += abs($AtomChargeHash{$Atom->getName} - $CopiedAtomChargeHash{$Atom->getName} ) ;
	  	
	    $Count++;
	
	    $TotalCount++;
	    
	
    }

  }

  $String .= "Total for Molecule: " . $TotalForMolecule . "\n";

  print $String if $TotalForMolecule > 0.09;
	
}

close(FILE);


print $TotalCount . " " . $Count . "\n";

my $r = 1;

if($TotalCount != 0 ) {

  eval { $r = ($xysum - ($xsum*$ysum)/$TotalCount) / (sqrt($xsqsum - $xsum**2 / $TotalCount) * sqrt($ysqsum - $ysum**2 / $TotalCount)) };
  if($@) {
	
	  $r = 1;
	
  }

}

#$TotalChargeDiff *= 100;

open(FILE, ">>molecule.txt");

print FILE "$Name1 $Name2 MOL:" . $MoleculeCount . " NUM:" . $TotalCount .  " USE:" . $TotalError / $TotalCount . " PUSE:"  . ($TotalChargeDiff / $Count)*100 . " R2:" . $r**2 . "\n" if $Count != 0;

close(FILE);




