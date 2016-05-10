package Complex;

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
use BaseObject;
use Chain;

require Exporter;

our @ISA = qw(BaseObject Exporter);

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

=head2 New

Usage: Complex->New;

Arguments:
  $Class: should be 'Complex'

Synopsis:
  creates a new Complex object which is the Highest level of datastructure organization, it contains all informatin about a system, 
all chains 

=cut

sub New {
	
	my ($Class, $StructuralHash) = @_;
	
	
	my $Self = {
		
		_Chains      => [ ], 
		_FileName		 => undef,
		_Name				 => undef,
		
	};
	
	bless $Self, $Class;
	
  $Self->BuildFromStructuralHash($StructuralHash) if defined $StructuralHash;
	
	return $Self;
	
}

=head2 BuildFromStructuralHash

Usage: $ComplexInstance->BuildFromStructuralHash;

Arguments:
  $StructuralHash: should be 'Complex'

Synopsis:
  Loads data from StructuralHash built from loading in from file

=cut

sub BuildFromStructuralHash {
	
	my ($Self, $StructuralHash) = @_;
	
	confess("This Hash is not of a Complex") if $StructuralHash->{'Structure'} ne ref($Self);
	
	my $Chains;
	
	$Chains = $StructuralHash->{'Chains'} if exists $StructuralHash->{'Chains'};
	
	my @ChainObjects;
			
	while( my ($ChainNumber, $ChainHash) = each(%$Chains)) {
				
		my $ChainObject = Chain->New($ChainHash);
		
		$ChainObject->setNum($ChainNumber);
		
		push @ChainObjects, $ChainObject;
				
	}	
	
	$Self->setChains(\@ChainObjects);
	
	
}

=head2 IdentifyChains

Usage: $ComplexInstance->IdentifyChains;

Arguments:

Synopsis:
  Tries to identify each Chain as RNA, DNA, Proteins, Water, or Ligand, Work in progress!

=cut

sub IdentifyChains {
	
	my $Self = shift;
	
	my $Chains = $Self->getChains;
	
	$_->IdentifyThisChain foreach (@$Chains);
	
}

=head2 getChainsOfType

Usage: $ComplexInstance->getChainsOfType($Type);

Arguments:
  $Type: The type of chain you wish to recover, Protein, Nucleic, Ligand, Water etc

Synopsis:
  returns the chains of specific type that are of interest

=cut

sub getChainsOfType {
	
	my ($Self, $Type) = @_;
	
	my $Chains = $Self->getChains;
	
	my @ChainsOfThisType =  grep { $_->getType && $_->getType eq $Type } @$Chains;
				
	return \@ChainsOfThisType;
	
}




1;
__END__
