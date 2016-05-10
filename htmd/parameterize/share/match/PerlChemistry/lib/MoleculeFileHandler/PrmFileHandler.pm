package MoleculeFileHandler::PrmFileHandler;

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
use MoleculeFileHandler;
use Data::Dumper;

require Exporter;

our @ISA = qw(MoleculeFileHandler Exporter);

# Items to export into callers namespace by default. Note: do not export
# names by default without a very good reason. Use EXPORT_OK instead.
# Do not simply export all your public functions/methods/constants.

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

=head2 Initiate;

Usage: $RtfFileHandlerInstance->Initiate;

Arguments:

Synopsis:
  Setups $RtfFileHandlerInstance must be called

=cut

=head2 Initiate;

Usage: $SdfFileHandlerInstance->Initiate;

Arguments:

Synopsis:
  Setups $SdfFileHandlerInstance must be called

=cut

sub Initiate {
	
	my $Self  = shift;

  my $Bond = qr/^\s*(\w+)\s+(\w+)\s+(\d+\.\d+)\s+(\d+\.\d+)/;
	my @BondVariables = qw(Types Kb b0);

	my $Angle = qr/^\s*(\w+)\s+(\w+)\s+(\w+)\s+(\d+\.\d+)\s+(\d+\.\d+)/;
	my @AngleVariables = qw(Types Ktheta Theta0);

	my $Dihedral = qr/^\s*(\w+)\s+(\w+)\s+(\w+)\s+(\w+)\s+(\d+\.\d+)\s+(1|2|3|4|5|6)\s+(\d+\.\d+)/;
	my @DihedralVariables = qw(Types Kchi n delta);

  my $Improper = qr/^\s*(\w+)\s+(\w+)\s+(\w+)\s+(\w+)\s+(\d+\.\d+)\s+(0)\s+(\d+\.\d+)/;
	my @ImproperVariables = qw(Types Kchi n delta);
	
	my $Nonbond = qr/^\s*(\w+)\s+(?:\d+.\d+|\d+)\s+(-?\d+\.\d+)\s+(\d+\.\d+)/;
	my @NonbondVariables = qw(Type Epsilon Rmin );


	my $BondLineParser = LineParser->New('Bond',
																	     sub {  return () if $_[0] !~ $Bond;
																		          return ([$1, $2], $3, $4) },
																	     sub { return { map {$BondVariables[$_] => $_[0]->[$_] } (0 .. $#BondVariables) } } );

  my $AngleLineParser = LineParser->New('Angle',
																			  sub {  return () if $_[0] !~ $Angle;
																							 return ([$1, $2, $3], $4, $5) },
																				sub { return { map {$AngleVariables[$_] => $_[0]->[$_] } (0 .. $#AngleVariables) } } );

	my $DihedralLineParser = LineParser->New('Dihedral',
																					 sub {  return () if $_[0] !~ $Dihedral;
																					        return ([$1, $2, $3, $4], $5, $6, $7) },
																					 sub { return { map {$DihedralVariables[$_] => $_[0]->[$_] } (0 .. $#DihedralVariables) } } );

	my $ImproperLineParser = LineParser->New('Improper',
																					 sub {  return () if $_[0] !~ $Improper;
																									return ([$1, $2, $3, $4], $5, $6, $7) },
																					 sub { return { map {$ImproperVariables[$_] => $_[0]->[$_] } (0 .. $#ImproperVariables) } } );

  my $NonbondLineParser = LineParser->New('Nonbond',
																				  sub {  return () if $_[0] !~ $Nonbond;
																								 return ($1, $2, $3) },
																					sub { return { map {$NonbondVariables[$_] => $_[0]->[$_] } (0 .. $#NonbondVariables) } } );
	

  $Self->{_Patterns} = [ $BondLineParser,
												 $AngleLineParser,
												 $DihedralLineParser,
												 $ImproperLineParser,
												 $NonbondLineParser];
												
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

	my $Patterns = $Self->getPatterns;

	my $Sections;
			
	my $CurrentRepresentation = "";
	
	my $DoneReading = 0;
	
	while(my $Line = <FILE>) {
		
	  if(eof(FILE)) { $CurrentRepresentation = 'EOF'; last; }
											
		foreach my $LineParser (@$Patterns) {
			
			next if ! $LineParser->Parse($Line);
			
			my $DataHash = $LineParser->ReturnResults;
							
			push @{$Sections->{$LineParser->getRepresentation}}, $DataHash if $DataHash;
										
		  $CurrentRepresentation = $LineParser->getRepresentation;
		  			
			last;
		  					
		}		
											
	}
	
	close(FILE);
	
	$Self->{_Sections} = $Sections;
	
	
}







1;
__END__