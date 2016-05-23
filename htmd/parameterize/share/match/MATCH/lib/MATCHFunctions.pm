package MATCHFunctions;

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
use BaseObject ':vars';
use MoleculeFileHandler;

require Exporter;

our @ISA = qw(BaseObject Exporter);

# Items to export into callers namespace by default. Note: do not export
# names by default without a very good reason. Use EXPORT_OK instead.
# Do not simply export all your public functions/methods/constants.

# This allows declaration	use BaseObject ':all';
# If you do not need this, moving things directly into @EXPORT or @EXPORT_OK
# will save memory.
our %EXPORT_TAGS = ( 'all'  => [ qw(ApplyStringNamingShortCuts SetupMoleculesFromTopology LoadAtomicParametersFromFile ShortenType)],
											 
										 'vars' => [ qw () ],
										
										 'func' => [ qw (ApplyStringNamingShortCuts SetupMoleculesFromTopology LoadAtomicParametersFromFile ShortenType)]);
							

our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

our @EXPORT = qw(
	
);

our $VERSION = '0.01';

#Exported Variables

our $ShortenTypeNameLookUp;

#Non Exported Variables


# Preloaded methods go here.

sub ShortenType {
	
  my $Type = shift;

  return $Type if uc($Parameters->getForceField) ne "CGENFF";
	
	if(! defined $ShortenTypeNameLookUp) {
		
		my $ShortenTypeNameFilePath = $Parameters->getShortenTypeFilePath;
		
		open(FILE, $ShortenTypeNameFilePath);
		
		$ShortenTypeNameLookUp = { map { my @spl = split /\s+/, $_; $spl[0] => $spl[1] } <FILE> };
		
		close(FILE);
		
	}
	
	
	return exists $ShortenTypeNameLookUp->{$Type} ? $ShortenTypeNameLookUp->{$Type} : $Type ;
	
}


sub LoadAtomicParametersFromFile {
	
	my $ParameterFilePath = shift;
	
  my $MoleculeFileHandler = MoleculeFileHandler->New($ParameterFilePath);
	
	return $MoleculeFileHandler->getSections;
	
}


sub SetupMoleculesFromTopology {
	
	my ($TopologyPath, $PatchingFile) = @_;
	
	my $MoleculeFileHandler = MoleculeFileHandler->New($TopologyPath);

	my $Structures = $MoleculeFileHandler->BuildObjectsFromFile;

  my @PatchingFileContents;
	
  if(defined $PatchingFile) {

	  open(FILE, $PatchingFile);
	  @PatchingFileContents = <FILE>;
	  close(FILE);
	
  }

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
				
		push @FinishedChains, $Chains[0];	
			
	}  

	foreach my $MoleculeChain (@MoleculeChains) {
		
		#next if $MoleculeChain->getMolecule(0)->getName ne "ALLOSE";
		
	#	print $MoleculeChain->getMolecule(0)->getName . "\n";
				
    next if $ExcludeFromPatching{$MoleculeChain->getMolecule(0)->getName};

    my $Patches = exists $DefaultPatches{"Default-" . $MoleculeChain->getMolecule(0)->getName} ? $DefaultPatches{"Default-" . $MoleculeChain->getMolecule(0)->getName} : $DefaultPatches{"Default"};

    foreach my $PatchName (@$Patches) {

      my $Patch = exists $PatchHash{$PatchName} ? $PatchHash{$PatchName} : confess("$PatchName does not exist in this Topology");

	    $MoleculeChain->getMolecule(0)->PatchTopologyMoleculeWithPres($Patch->getMolecule(0));

    }

    $MoleculeChain->Initiate;

    if($MoleculeChain->getMolecule(0)->getName eq "HSP,NTER,CTER") {
	
	    #print "made it\n"; 
	
	   # my %AtomsByName = map { $_->getName => $_ } @{$MoleculeChain->getMolecule(0)->getAtoms};
	
	    #print "START!!!\n";
	
	    #print $AtomsByName{'NE2'}->DoesLookUpTableMatch($AtomsByName{'ND1'}->getLookUpTable) . "\n";
	
	    foreach my $Atom (@{$MoleculeChain->getMolecule(0)->getAtoms}) {
		
		  #  print $Atom->getName . " " . $Atom->getNumOfNotSolvedBonds . " : ";
		
		    foreach my $Bonded ( @{$Atom->getBondedAtoms}) {
			
			 #   print $Bonded->getName . " ";
			
		    }
		
		   #print "\n";
		
	    }
	   # exit 1;
	
    }

    push @FinishedChains, $MoleculeChain;

  }

	
 
  return \@FinishedChains;

	
	
}


sub ApplyStringNamingShortCuts {
	
	my ($ShortCuts, $Lines) = @_;
	
	my %ShortCutHash;
	
	#Short cuts are defined currently in .trans and .types files in the format ShortCut = Value
	#Short cuts are used to decluter Type strings
	foreach my $ShortCut (@$ShortCuts) {
		
	  chomp($ShortCut);
	
		my @spl = split /\s+\=\s+/, $ShortCut;
	 
		$ShortCutHash{$spl[0]} = $spl[1]; 
	
	}
	
	#Go through all the no short cut lines and find where a short cut is used and replace it with its values
	foreach my $Line (@$Lines) {
	
	  foreach my $ShortCut (keys %ShortCutHash) {
	  
	  	next if $Line !~ $ShortCut;
	  	
	  	my $Value = $ShortCutHash{$ShortCut};
	  		
	  	$Line =~ s/$ShortCut/$Value/g;
	  
	    #print $Line; 
	  
	  }
	
	}
		
	return $Lines;
	
}

1;
__END__
