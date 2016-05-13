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
use AtomTyper;
use Storable;

our $VERSION = '0.01';

#Exported Variables


#Non Exported Variables


# Preloaded methods go here.

srand(time ^ ($$ + ($$ << 15)));

#Deal With Initial Parameter Setup, This will be handled in main .pl script
my $DefaultParameters = MATCHParameters->New;

#$DefaultParameters->Initiate("Top22Parameters.par");
$DefaultParameters->Initiate("../resources/" . $ARGV[0] . ".par");

$Parameters = $DefaultParameters;

my $AtomTyper = AtomTyper->New;

$AtomTyper->Initiate;

#my $Chains = SetupMoleculesFromTopology("../resources/top_all22_prot/top_all22_prot.inp",
# 																		    "../resources/top_all22_prot/top_all22s_prot.patches");

my $Chains = retrieve("StoredTopologys/" . $ARGV[0] . ".dat");

#my $Chains = retrieve("../t/StoredTopologys/top_all22_prot.dat");

my @Atoms = map { @{$_->getAtoms} } map { @{$_->getMolecules} } @$Chains;


my %DoNotType =  ( 'CT1X' => 1, 'CT2X' => 1, 'CT3X' => 1, 'HA' => 1 ); 


foreach my $Atom (@Atoms) {
	
	#next if $Atom->getMissingCharge == 0;
	
	#print $Atom->getName . " " . $Atom->getType . " " . $Atom->getMolecule->getName . " " . $Atom->getMissingCharge . "\n";
	
}


my @ElementAtoms = grep { $_->getState eq "C.3"} @Atoms;

#my @ElementAtoms = grep { $_->getElement eq "F" } @Atoms;


my @UniqueTypes = keys %{ {map { $_->getType => 1 } @ElementAtoms} };

print join(" ", @UniqueTypes) . "\n";

my $Types = $AtomTyper->getTypes;

foreach my $Type (@$Types) {
	
	#print $Type->getName . " " . $Type->getString . "\n";
	
}

#print "START!\n\n";

my $CG2DCFlag = 0;

foreach my $Atom (@ElementAtoms) {
		
	my $Flag =0;
			
	#print  $Atom->getType . " " . $Atom->getName . " " .  $Atom->getMolecule->getName . " " . $Atom->Stringify . "\n";
	
	#next if $Atom->getMolecule->getName ne "PRAC";
	
	foreach my $Type (@$Types) {
		
		#print  $Type->getName . "\n";
		
		#next if $Type->getName ne "CG2DC3";	
					
		if($Type->getLookUpTable->AreAllNodesSharedBetween([$Atom->getLookUpTable])) {
			
			if($Type->getName eq "CG2DC?" && ($Atom->getType eq "CG2DC1" || $Atom->getType eq "CG2DC2")) {
				
				$Atom->{_TestType} = "CG2DC?";
				
				$CG2DCFlag = 1;
				$Flag = 1;
				last;
				
			}
			
			$Atom->{_TestType} = $Type->getName;
			
			#print $Atom->getType . " " . $Type->getName . " " . $Atom->getName . " " .  $Atom->getMolecule->getName . " " . $Atom->Stringify . " " . $Type->getString . "\n" if $Atom->getType ne $Type->getName;
			
			$Flag = 1;
			
			last;
			
		}
		
	}
	
	if(!$Flag) {
	
	  print "NOT MATCHED: " . $Atom->getType . " " . $Atom->getName . " " .  $Atom->getMolecule->getName . " " . $Atom->Stringify . "\n";
	
  }
	
}


exit 1;
#exit 1 if $CG2DCFlag == 0;

my @Chains;
my %SeenCG2DC;

my %AcceptableEnd = map { $_ => 1 } qw(CG2D1O CG2D2O CG2DC3);

foreach my $Atom (@ElementAtoms) {
	
	next if $Atom->getTestType ne "CG2DC?" || exists $SeenCG2DC{$Atom};
		
	$SeenCG2DC{$Atom} = 1;
		
	my @Chain = ($Atom);
	
	my $Finished = 0;
	
	while(!$Finished) {
	
	  $Finished = 1;
	
	  my $StartBonds = $Chain[0]->getBondedAtoms;
	
	  foreach my $Bonded (@$StartBonds) {
		
		  $Bonded->{_TestType} = $Bonded->getType if ! defined $Bonded->getTestType;
		
		  next if $Bonded->getTestType ne "CG2DC?" || exists $SeenCG2DC{$Bonded};
		
		  unshift @Chain, $Bonded; $Finished = 0;
		
		  $SeenCG2DC{$Bonded} = 1;
		 
		  last;
		
	  }
	
	  my $EndBonds = $Chain[-1]->getBondedAtoms;
	
	  foreach my $Bonded (@$EndBonds) {
		
		  $Bonded->{_TestType} = $Bonded->getType if ! defined $Bonded->getTestType;
		
		  next if $Bonded->getTestType ne "CG2DC?" || exists $SeenCG2DC{$Bonded};
		
		  push @Chain, $Bonded; $Finished = 0;
		
		  $SeenCG2DC{$Bonded} = 1;
		 
		  last;
		
	  }  	
		
	}
	
	my $Start = ""; my $End = "";
	
	foreach my $Bonded (@{$Chain[0]->getBondedAtoms}) {
		
		if(exists $AcceptableEnd{$Bonded->getType}) {
		
		  $Start = $Bonded->getType;
			
		}
		
	}

	foreach my $Bonded (@{$Chain[-1]->getBondedAtoms}) {
		
		if(exists $AcceptableEnd{$Bonded->getType}) {
			
			$End = $Bonded->getType;
			
		}
		
	}
	
	push @Chains, \@Chain;
	
	if (@Chain == 1) {
		
		$Chain[0]->setTestType("CG2DC1");
				
	}
	
	elsif($Start eq "CG2DC3" && $End eq "CG2DC3") {
		
		if($Chain[0]->getName ge $Chain[-1]->getName) {
						
			@Chain = reverse(@Chain);
			
		}
		
		$Chain[0]->setTestType("CG2DC2");
		
	}
	
	elsif($End eq "CG2DC3") {
		
	  @Chain = reverse(@Chain);
	
		$Chain[0]->setTestType("CG2DC2");
		
	}
	
	elsif($Start eq "CG2D1O") {
		
		$Chain[0]->setTestType("CG2DC1");
		
	}
	
	elsif($End eq "CG2D1O") {
		
		@Chain = reverse(@Chain);
		
		$Chain[0]->setTestType("CG2DC1");
		
	}
	
	$Chain[0]->setTestType("CG2DC2") if $Chain[0]->getTestType eq "CG2DC?";
	
	AlternateCG2DC(@Chain);
	
	my $String = "";
	
	my $Print = 0;
	
	$String .= $Atom->getMolecule->getName . " ";
	
	foreach my $Bonded (@{$Chain[0]->getBondedAtoms}) {
		
		if(exists $AcceptableEnd{$Bonded->getType}) {
			
			$String .= $Bonded->getName . " " . $Bonded->getType . " ";
		
		  $String .= $Bonded->getBondwithAtom($Chain[0])->getType . " ";
		
		  $Start = $Bonded->getType;
			
		}
		
	}
	
	foreach my $i (0 .. @Chain-1) {
		
		next if $Chain[$i]->getType eq $Chain[$i]->getTestType;
		
		$Print = 1;
		
		$String .= $Chain[$i]->getName . " " . $Chain[$i]->getType . " " .  $Chain[$i]->getTestType . " ";
		
		if($i < @Chain-1) {
			
			$String .= $Chain[$i]->getBondwithAtom($Chain[$i+1])->getType . " ";
			
		}
		
	}
	
	foreach my $Bonded (@{$Chain[-1]->getBondedAtoms}) {
		
		if(exists $AcceptableEnd{$Bonded->getType}) {
			
			$String .= $Bonded->getName . " " . $Bonded->getType . " ";
			
			$String .= $Bonded->getBondwithAtom($Chain[-1])->getType . " ";
		  
			$End = $Bonded->getType;
			
		}
		
	}
	
  print $String . "\n" if $Print;
	
	
}

sub AlternateCG2DC {
	
	my (@Chain) = @_;
	
	foreach my $i (1 .. @Chain-1) {
		
    my $Same; my $Alt;
		
		if($Chain[$i-1]->getTestType eq "CG2DC1") {
			
			$Same = "CG2DC1"; $Alt = "CG2DC2";
			
		} 
		
		else {
			
			$Same = "CG2DC2"; $Alt = "CG2DC1";
			
		}

    my $CurrentBondType = $Chain[$i-1]->getBondwithAtom($Chain[$i])->getType;

		if($CurrentBondType == 2) {
			
			$Chain[$i]->setTestType($Same);
			
		}
		
		else {
			
			$Chain[$i]->setTestType($Alt);
			
		}
		
	}
	
}


#If only one then CG2DC1
#If connected to CG2D1O => CG2DC1
#If not only one and connected to CG2DC3 => CG2DC2, except MECH
#if only 2 and both ends are CG2DC3, then CG2DC1 & CG2DC2
 
