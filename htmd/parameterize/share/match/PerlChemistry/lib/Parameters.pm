package Parameters;

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
use BaseObject ':all';

require Exporter;

our @ISA = qw(BaseObject Exporter);

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
	
	my $Class = shift;
	
  my $Self = {
	  _PATH 							       => undef,
	  _BondLengthFilePath        => undef,
		_ExitifNotInitiated	       => 1,
		_CheckAtomProtState        => 1,
		_CheckElementBondingNumber => 1,
		_RenameAtoms               => 1,
	
	};
	
	bless $Self, $Class;
}
	

sub Initiate {
	
  my ($Self, $ParameterFilePath, @CommandLineParameters) = @_;

  my $DefaultPath = $ENV{'PerlChemistry'};

  $ParameterFilePath = $ENV{'PerlChemistry'} . "/resources/DefaultParameters.par" unless defined $ParameterFilePath;

  $Self->ReadParameterFile($DefaultPath, $ParameterFilePath);

  $Self->ProcessComandLineArguments($DefaultPath, @CommandLineParameters);

}

sub ProcessComandLineArguments {
	
	my ($Self, $DefaultPath, @CommandLineParameters) = @_;
	
	foreach my $CommandLineParameter (@CommandLineParameters) {
		
	  $Self->ProcessParameter($DefaultPath, substr($CommandLineParameter->[0],1), $CommandLineParameter->[1]);
		
	}
	
}

sub ReadParameterFile {
	
	my ($Self, $DefaultPath, $ParameterFilePath) = @_;
	
	confess("Path to ParameterFile does not exist") unless -e $ParameterFilePath;
	
	open(FILE, $ParameterFilePath);

  my %ParameterFileContents = map { my ($ParameterName, $Value) = ($_ =~ /\s*(\S+)\s*=\s*(\S+)/ );  $ParameterName => $Value } 
 														  grep { $_ =~ /^[^#]\s*(\S+)\s*=\s*(\S+)/ } <FILE>;
 
 close(FILE);

  while (my ($Parameter, $Value) = each %ParameterFileContents) {
	
    $Self->ProcessParameter($DefaultPath, $Parameter, $Value);  
	
	}
	
}

sub ProcessParameter {
	
	my ($Self, $DefaultPath, $Parameter, $Value) = @_;
	
	if($Value =~ /\// ) {
	
	  unless(-e $Value) {
		
		  confess("Parameter: " . $Parameter. " File does not exist: $Value")  unless(-e "$DefaultPath/$Value");
	    
	    $Value = "$DefaultPath/$Value";
		
	  }
	
  }

  PrintMessage("Warning: $Parameter is an unknown parameter!", 1) unless exists $Self->{"_" . $Parameter};

  $Self->{"_" . $Parameter} = $Value;
	
}




1;
__END__
