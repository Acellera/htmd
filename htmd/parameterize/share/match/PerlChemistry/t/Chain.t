#!/usr/bin/perl -w

use strict;
use warnings;

use lib $ENV{'PerlChemistry'} . "/lib";

use Test::More tests => 1;

BEGIN { use_ok('Chain') };

use Chain;
use MoleculeFileHandler;

my ($TopologyPath, $PatchingFile) = ("../resources/top_all22_prot/top_all22_prot.inp", "../resources/top_all22_prot/top_all22_prot.patches");

my $MoleculeFileHandler = MoleculeFileHandler->New($TopologyPath);

my $Structures = $MoleculeFileHandler->BuildObjectsFromFile;

open(FILE, $PatchingFile);
my @PatchingFileContents = <FILE>;
close(FILE);

my @MoleculeChains = grep { !$_->getMolecule(0)->getIsPatch }  @$Structures;

my @PatchChains = grep { $_->getMolecule(0)->getIsPatch } @$Structures;

my %MoleculeHash = map { $_->getMolecule(0)->getName => $_ } @MoleculeChains;

my %PatchHash = map { $_->getMolecule(0)->getName => $_ } @PatchChains;

my %ExcludeFromPatching;
my %DefaultPatches;
my @PatchInstructions;

foreach (@PatchingFileContents) { 

	my @spl = split /\s+/ , $_;

	if($spl[0] eq 'Exclude') { $ExcludeFromPatching{$spl[$_]} = 1 foreach( 1 .. $#spl) }

	elsif($spl[0] =~ /^Default/) { $DefaultPatches{shift @spl } = \@spl }

	else { push @PatchInstructions, \@spl }

}

my @FinishedChains;
	
foreach (@PatchInstructions) {
		
  my ($MoleculeNamesInvolvedInPatch, @PatchesToBeApplied) = @{$_};

	my @Chains;
	
	foreach my $MoleculeName ( split(/\-/, $MoleculeNamesInvolvedInPatch) ) {

	  push @Chains, exists $MoleculeHash{$MoleculeName} ? $MoleculeHash{$MoleculeName}->Copy : confess("$MoleculeName does not exist in this Topology");

	}
	
	next if @Chains < 2;
			
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

  print $Chains[0]->getName . "\n";
			
	push @FinishedChains, $Chains[0];	
			
}  
