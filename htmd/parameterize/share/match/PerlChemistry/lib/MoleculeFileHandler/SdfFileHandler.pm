package MoleculeFileHandler::SdfFileHandler;

=head1 SYNOPSIS

MoleculeFileHandler::SdfFileHandler - Handles loading in SDF formated files  

=head1 DESCRIPTION

use MoleculeFileHandler::SdfFileHandler;
@ISA = qw(MoleculeFileHandler);

This filehandler also currently opens MOL format files too

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

#use 5.010000;
use strict;
use warnings;
use Carp;
use MoleculeFileHandler;

require Exporter;

our @ISA = qw(MoleculeFileHandler Exporter);

our $VERSION = '0.01';

#Exported Variables


#Non Exported Variables


# Preloaded methods go here.

sub New {
	
	#Use MoleculeFileHandler->New
	
}

=head2 Initiate;

Usage: $SdfFileHandlerInstance->Initiate;

Arguments:

Synopsis:
  Setups $SdfFileHandlerInstance must be called to setup the neccessary LineParser objects to process the SDF file

=cut

sub Initiate {
	
	my $Self  = shift;
  
	my $EndPattern = qr/\$\$\$\$/;

  my $Atom = qr/^\s*(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(\w+)\s+\d\s+\d\s+\d\s+\d/;
  my @AtomVariables = qw(X Y Z Element);

  my $Bond = qr/^\s*(\d+)\s+(\d+)\s+(\d+)\s+\d+\s+\d+\s+\d+\s+\d+\s*\n/;
	my @BondVariables = qw(BondAtomNums Type);
	
	my $BondBackUp = qr/^\s*(\d+)\s+(\d+)\s+(\d+)\s+\d\s*\n/;
	my @BondVariablesBackUp = qw(BondAtomNums Type);
		
	my $AtomLineParser = LineParser->New('Atom',
																    	 sub { $_[0] =~ $Atom },
																	     sub { return { map {$AtomVariables[$_] => $_[0]->[$_] } (0 .. $#AtomVariables) } } );
	
	my $BondLineParser = LineParser->New('Bond',
																	     sub {  return () if $_[0] !~ $Bond;
																		          return ([$1, $2], $3) },
																	     sub { return { map {$BondVariables[$_] => $_[0]->[$_] } (0 .. $#BondVariables) } } );
	
	my $BondBackUpLineParser = LineParser->New('Bond',
																			 sub {  return () if $_[0] !~ $BondBackUp;
																						  return ([$1, $2], $3) },
																			 sub { return { map {$BondVariablesBackUp[$_] => $_[0]->[$_] } (0 .. $#BondVariablesBackUp) } } );																
		
	my $EndOfChainParser = LineParser->New('EndChain',
																			 sub { $_[0] =~ $EndPattern },
																			 sub {  } );	

  $Self->{_Patterns} = [ $AtomLineParser,
												 $BondLineParser,
												 $BondBackUpLineParser,
												 $EndOfChainParser];											
}

=head2 ReadFile;

Usage: $SdfFileHandlerInstance->ReadFile($File);

Arguments:

Synopsis:
  This is where all the data is collected from File

=cut

sub ReadFile {
		
	my ($Self, $FilePath) = @_;
	
	open(FILE, $FilePath);
	
	my @Chains;
	my $Finished = 0;
	
	my $Count = 0;
	
	while(!$Finished) {
	  push @Chains, $Self->ReadChainFromFile(*FILE); 
	 		
	  $Finished = 1 if($Chains[-1]->[0] eq 'EOF');
  }
	
	close(FILE);
		
	my @FinishedChains;
	
	my $ResidueCount = 0;

	foreach my $Chain (@Chains) {

    my $Atoms = $Chain->[1]->{'Atom'};

    #Makes sure that if the $$$$ divider is not at the end it doesnt create a dummy chain with no atoms
    next if ! defined $Atoms;

    $Count = 1;

    #Atom names in SDF/MOL files are just the element such as every carbon atom is called C, this gives each atom a unique name
    #In the format Element + Num
    foreach my $Atom (@$Atoms) {
	
		  $Atom->{'Name'} = $Atom->{'Element'} . $Count++;
		
    }

    my $ResidueIndentifier = "UNK" . " " . $ResidueCount;
		
		#Preprocess molecule information since there is only one molecule
		$Chain->[1]->{'Molecules'}{$ResidueIndentifier}{'Structure'} = 'Molecule';
	  $Chain->[1]->{'Molecules'}{$ResidueIndentifier}{'Bonds'} = delete $Chain->[1]->{'Bond'};
	  $Chain->[1]->{'Molecules'}{$ResidueIndentifier}{'Atoms'} = delete $Chain->[1]->{'Atom'};
	  $Chain->[1]->{'Molecules'}{$ResidueIndentifier}{'Name'} =  "UNK";
	  $Chain->[1]->{'Molecules'}{$ResidueIndentifier}{'Num'} = $ResidueCount;
	
	  push @FinishedChains, $Chain->[1];
	
	  $ResidueCount++;
	
	}
	
	$Self->setChains(\@FinishedChains);
	
  return 1;
	
}

=head2 ReadChainFromFile;

Usage: $SdfFileHandlerInstance->ReadChainFromFile($File);

Arguments:

Synopsis:
  Collects the data on a specific Chain, this allows the break up of the data

=cut

sub ReadChainFromFile {
	
		my ($Self, $FH) = @_;
		
		my $Patterns = $Self->getPatterns;

		my $Sections;
				
		my $CurrentRepresentation = "";
		
		my $DoneReading = 0;
		
		while(my $Line = <$FH>) {
			
		  if(eof($FH)) { $CurrentRepresentation = 'EOF'; last; }
												
			foreach my $LineParser (@$Patterns) {
				
				next if ! $LineParser->Parse($Line);
				
				my $DataHash = $LineParser->ReturnResults;
								
				push @{$Sections->{$LineParser->getRepresentation}}, $DataHash if $DataHash;
											
			  $CurrentRepresentation = $LineParser->getRepresentation;
			  			
				last;
			  					
			}		
											
			last if $CurrentRepresentation eq 'EndChain' or $CurrentRepresentation eq 'EndOfComplex' ;
		
		}		
			
		return [$CurrentRepresentation, $Sections];
}






1;
__END__
