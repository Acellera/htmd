package MoleculeFileHandler::Mol2FileHandler;
use Data::Dumper;
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
use MoleculeFileHandler;

require Exporter;

our @ISA = qw(MoleculeFileHandler Exporter);

# Items to export into callers namespace by default. Note: do not export
# names by default without a very good reason. Use EXPORT_OK instead.
# Do not simply export all your public functions/methods/constants.

# This allows declaration	use BaseObject ':all';
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


# Preloaded methods go here.

sub New {
	
}

=head2 BuildFromFile;

Usage: $MoleculeFileHandlerInstance->BuildFromFile;

Arguments:

Synopsis:
  Creates the neccessary Objects to represent the Data in the file, i.e. if its only one molecule, this function will return a new molecule object 
with all the data stored in it, see MoleculeFileHandler.t if interested

=cut

sub LoadDataFromFile {
	
}

=head2 Initiate;

Usage: $Mol2FileHandlerInstance->Initiate;

Arguments:

Synopsis:
  Setups $Mol2FileHandlerInstance must be called

=cut

sub Initiate {
	
	my $Self  = shift;
	
	my $EndPattern = qr/\$\$\$\$/;

  my $AtomWithCharge = qr/^\s*(\d+)\s+(\S+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+\S+\s+(\d+)\s+(\S+)\s+(-?\d+\.\d+)/;
  my @AtomWithChargeVariables = qw(Num Name X Y Z ResidueNum ResidueName Charge);

  my $Atom = qr/^\s*(\d+)\s+(\S+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+\S+\s+(\d+)\s+(\S+)/;
  my @AtomVariables = qw(Num Name X Y Z ResidueNum ResidueName);

  my $Bond = qr/^\s*(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s*\n/;
	my @BondVariables = qw(Num BondAtomNums Type);
	
	my $AtomWithChargeLineParser = LineParser->New('Atom',
																    	 sub { $_[0] =~ $AtomWithCharge },
																	     sub { return { map {$AtomWithChargeVariables[$_] => $_[0]->[$_] } (0 .. $#AtomWithChargeVariables) } } );
	
	my $AtomLineParser = LineParser->New('Atom',
																    	 sub { $_[0] =~ $Atom },
																	     sub { return { map {$AtomVariables[$_] => $_[0]->[$_] } (0 .. $#AtomVariables) } } );
	
	my $BondLineParser = LineParser->New('Bond',
																	     sub {  return () if $_[0] !~ $Bond;
																		          return ($1, [$2, $3], $4) },
																	     sub { return { map {$BondVariables[$_] => $_[0]->[$_] } (0 .. $#BondVariables) } } );
																	
		
	my $EndOfChainParser = LineParser->New('EndChain',
																			 sub { $_[0] =~ $EndPattern },
																			 sub {  } );	

  $Self->{_Patterns} = [ $AtomWithChargeLineParser,
 												# $AtomLineParser,
												 $BondLineParser,
												 $EndOfChainParser];
												
}

=head2 ReadFile;

Usage: $Mol2FileHandlerInstance->ReadFile($File);

Arguments:

Synopsis:
  This is where all the data is collected from File

=cut

sub ReadFile {
		
	my ($Self, $FilePath) = @_;
	
	open(FILE, $FilePath);
	
	my @Chains;
	my $Finished = 0;
		
	while(!$Finished) {
	  push @Chains, $Self->ReadChainFromFile(*FILE); 
	  $Finished = 1 if($Chains[-1]->[0] eq 'EOF');
  }
	
	close(FILE);
	
	my @FinishedChains;

	foreach my $Chain (@Chains) {

    my $Atoms = $Chain->[1]->{'Atom'};

    my $ResidueIndentifier = $Atoms->[0]->{'ResidueName'} . " " . $Atoms->[0]->{'ResidueNum'};
		
		$Chain->[1]->{'Molecules'}{$ResidueIndentifier}{'Atoms'} = $Atoms;
		$Chain->[1]->{'Molecules'}{$ResidueIndentifier}{'Structure'} = 'Molecule';
	  $Chain->[1]->{'Molecules'}{$ResidueIndentifier}{'Bonds'} = delete $Chain->[1]->{'Bond'};
	  $Chain->[1]->{'Molecules'}{$ResidueIndentifier}{'Name'} =  $Atoms->[0]->{'ResidueName'};
	  $Chain->[1]->{'Molecules'}{$ResidueIndentifier}{'Num'} = 1;
	
	  push @FinishedChains, $Chain->[1];
    	
	}
	
	$Self->setChains(\@FinishedChains);
  return 1;
	
}

=head2 ReadChainFromFile;

Usage: $Mol2FileHandlerInstance->ReadChainFromFile($File);

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

print Dumper ( $Sections );
			
		return [$CurrentRepresentation, $Sections];
}






1;
__END__
