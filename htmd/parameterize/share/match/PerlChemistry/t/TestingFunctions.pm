package TestingFunctions;

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

use 5.010000;
use strict;
use warnings;
use Carp;
use BaseObject;
use Molecule;

require Exporter;

our @ISA = qw(BaseObject Exporter);

# Items to export into callers namespace by default. Note: do not export
# names by default without a very good reason. Use EXPORT_OK instead.
# Do not simply export all your public functions/methods/constants.

# This allows declaration	use BaseObject ':all';
# If you do not need this, moving things directly into @EXPORT or @EXPORT_OK
# will save memory.
our %EXPORT_TAGS = ( 'all'  => [ qw(is_atom_copied)],
											 
										 'vars' => [ qw () ],
										
										 'func' => [ qw (is_atom_copied)]);
							

our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

our @EXPORT = qw(
	
);

our $VERSION = '0.01';

#Exported Variables


#Non Exported Variables


sub is_atom_copied {
	
	my ($Atom1, $Atom2) = @_;
	
	my %DoNotCopy = ( _Bonds => 1, _BondedAtoms => 1, _SumofBondTypes => 1, _Molecule => 1, _Chain => 1, _NumOfBonds => 1, _NumOfNotSolvedBonds => 1, _Rings => 1, _UniqueString => 1); 
	
	my $NotSame = 1;
	
	while ( my ($Field, $Value) = each(%$Atom1)) {
		
		next if exists $DoNotCopy{$Field};
				
    if(! exists $Atom2->{$Field} ) {
	
	    $NotSame = 0; last;
	
    }

    next if ! defined $Value and ! defined $Atom2->{$Field};

    next if $Value eq $Atom2->{$Field};

    $NotSame = 0; last;
		
	}
	
	return $NotSame;
	
}

sub is_bond_copied {
	
 	my ($Bond1, $Bond2) = @_;
	

}

