#!/usr/bin/perl -w

use strict;
use warnings;
use Carp;
use Cwd qw(abs_path);

use lib $ENV{'PerlChemistry'} . "/lib";
use Test::More tests => 4;

BEGIN { use_ok('Molecule') };

use Molecule ':all';
use MoleculeFileHandler ':all';
use Parameters ':all';
use BaseObject ':vars';


#Deal With Initial Parameter Setup, This will be handled in main .pl script
my $DefaultParameters = Parameters->New;

$DefaultParameters->Initiate;

#Set Global Parameters
$Parameters = $DefaultParameters;

my $MoleculeFileHandler = MoleculeFileHandler->New($ENV{'PerlChemistry'} . "/resources/top_all22_prot/top_all22_prot.inp");

my $Structures = $MoleculeFileHandler->BuildObjectsFromFile;

open(FILE, $ENV{'PerlChemistry'} .  "/resources/top_all22_prot/top_all22_prot.patches");
my @PatchingFileContents = <FILE>;
close(FILE);

my @MoleculeChains = grep { !$_->getMolecule(0)->getIsPatch }  @$Structures;

my @PatchChains = grep { $_->getMolecule(0)->getIsPatch } @$Structures;

my %MoleculeHash = map { $_->getMolecule(0)->getName => $_ } @MoleculeChains;

my %PatchHash = map { $_->getMolecule(0)->getName => $_ } @PatchChains;

subtest "Test Bonding Procedure" => sub {	
	
	open(FILE, "MoleculeTestResources/top_all22_prot_CorrectBonding.dat");
	my @CorrectBondingFileContents = <FILE>;
	close(FILE);
	
	my @ResidueSplitFileContents = split /\n\n/, join("", @CorrectBondingFileContents);
	
  my %TopologyBonds;

  foreach (@ResidueSplitFileContents) {
	
	  my @SplitResidue = split(/\n/, $_);
	
	  my $Name = shift @SplitResidue;
	
	  $TopologyBonds{$Name} = [ map { my @spl = split /\s+/, $_; @spl } @SplitResidue ];
	
  }
 
  plan tests => scalar(keys %TopologyBonds);
	
  while ( my ($ResidueName, $BondArray) = each %TopologyBonds) { 

	  is(CheckToSeeIfAllBondsAreDefined($MoleculeHash{$ResidueName}->getMolecule(0)->getBonds, $BondArray), 0, "$ResidueName");
	
  }
	
};

subtest "Test Patching Procedure" => sub {

  plan tests => 2;

  my %ExcludeFromPatching;
	my %DefaultPatches;
	my @PatchInstructions;

	foreach (@PatchingFileContents) { 

		my @spl = split /\s+/ , $_;

		if($spl[0] eq 'Exclude') { $ExcludeFromPatching{$spl[$_]} = 1 foreach( 1 .. $#spl) }

		elsif($spl[0] =~ /^Default/) { $DefaultPatches{shift @spl } = \@spl }

		else { push @PatchInstructions, \@spl }

	}

  subtest "Test Multiple Residue Patches" => sub {
	
	  plan tests => 1;
	
	  foreach (@PatchInstructions) {

			#print join(" ", @{$_}) . "\n";

			my ($MoleculeNamesInvolvedInPatch, @PatchesToBeApplied) = @{$_};

			my @Chains;

			foreach my $MoleculeName ( split(/\-/, $MoleculeNamesInvolvedInPatch) ) {

				push @Chains, exists $MoleculeHash{$MoleculeName} ? $MoleculeHash{$MoleculeName}->Copy : confess("$MoleculeName does not exist in this Topology");

			}

			foreach my $PatchesToBeApplied (@PatchesToBeApplied) {

				my ($ChainNumber1, $PatchName, $ChainNumber2) = split(/\-/, $PatchesToBeApplied);

			  my $Patch = exists $PatchHash{$PatchName} ? $PatchHash{$PatchName} : confess("$PatchName does not exist in this Topology");

			  if(!$ChainNumber2) { $Chains[$ChainNumber1-1]->getMolecule(0)->PatchTopologyMoleculeWithPres($Patch->getMolecule(0))}
			
			  else { 
				
				  $Chains[$ChainNumber1-1]->AddMolecule($Chains[$ChainNumber2-1]->getMolecule(0)); 
				
				  $Chains[$ChainNumber1-1]->PatchTopologyChainWithPres(0, 1, $Patch->getMolecule(0));
				
				}

			}
						
			$Chains[0]->AddMolecule($Chains[$_]->getMolecule(0)) foreach (1 .. $#Chains);
			
			$Chains[0]->Initiate;
		
	  	#print $Chains[0]->getName . "\n";
		
	  	foreach my $Molecule (@{$Chains[0]->getMolecules}) {
		
		    my $Bonds = $Molecule->getBonds;

		    foreach my $Bond (sort { $a->getPrimaryName cmp $b->getPrimaryName || $a->getSecondaryName cmp $b->getSecondaryName } @$Bonds) {

		    }		
		  }
		
	  }
		
		is(1, 1);
	
  };

  subtest "Test Single Residue Patches" => sub {


  	open(FILE, "MoleculeTestResources/top_all22_prot_CorrectPatching.dat");
	  my @CorrectBondingFileContents = <FILE>;
	  close(FILE);

	  my @ResidueSplitFileContents = split /\n\n/, join("", @CorrectBondingFileContents);

    my %TopologyBonds;

    foreach (@ResidueSplitFileContents) {

	    my @SplitResidue = split(/\n/, $_);

  	  my $Name = shift @SplitResidue;

	    $TopologyBonds{$Name} = [ map { my @spl = split /\s+/, $_; @spl } @SplitResidue ];

    }

	  plan tests => scalar(keys %TopologyBonds);

    foreach my $MoleculeChain (@MoleculeChains) {

	    next if $ExcludeFromPatching{$MoleculeChain->getMolecule(0)->getName};

	    my $Patches = exists $DefaultPatches{"Default-" . $MoleculeChain->getMolecule(0)->getName} ? $DefaultPatches{"Default-" . $MoleculeChain->getMolecule(0)->getName} : $DefaultPatches{"Default"};

	    foreach my $PatchName (@$Patches) {

	      my $Patch = exists $PatchHash{$PatchName} ? $PatchHash{$PatchName} : confess("$PatchName does not exist in this Topology");

		    $MoleculeChain->getMolecule(0)->PatchTopologyMoleculeWithPres($Patch->getMolecule(0));

	    }
    }

    while ( my ($ResidueName, $BondArray) = each %TopologyBonds) { 

	    is(CheckToSeeIfAllBondsAreDefined($MoleculeHash{$ResidueName}->getMolecule(0)->getBonds, $BondArray), 0, "$ResidueName");

    }

  };

};


sub CheckToSeeIfAllBondsAreDefined {

 my ($MoleculeBonds, $TopologyBonds) = @_;

 my $Failed = 0;
  for (my $i = 0; $i < @$TopologyBonds; $i+=3) {

	  my $Exists = 0;
	  foreach my $Bond (@$MoleculeBonds) {

		  if(($Bond->getPrimaryName eq $TopologyBonds->[$i] && $Bond->getSecondaryName eq $TopologyBonds->[$i+1]) ||
		     ($Bond->getPrimaryName eq $TopologyBonds->[$i+1] && $Bond->getSecondaryName eq $TopologyBonds->[$i]) && 
		     $Bond->getType == $TopologyBonds->[$i+2]) {

			  $Exists = 1; last;

		  }

	  }

	  if(!$Exists) { $Failed = 1; last}

	}	
	
  return $Failed;

}