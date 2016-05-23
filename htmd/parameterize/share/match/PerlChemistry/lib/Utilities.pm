package Utilities;

=head1 NAME

=head1 SYNOPSIS

=head1 DESCRIPTION

=head2 EXPORT

Nothing

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

# If you do not need this, moving things directly into @EXPORT or @EXPORT_OK
# will save memory.
our $VERSION = '0.01';

#Exported Variables

our %EXPORT_TAGS = ( 'all'  => [ qw (RotateStructure MoveStructure Magnitude CrossProduct DotProduct ArcCosine)],
											 
										 'vars' => [ qw () ],
										
										 'func' => [ qw (RotateStructure MoveStructure Magnitude CrossProduct DotProduct ArcCosine)]);
							

our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

our @EXPORT = qw(
	
);

sub Magnitude {
	my $Vector = shift;
	my $Total = 0;
	foreach my $Element (@$Vector) {
		$Total += $Element**2;
	}
	return sqrt($Total);
}

sub CrossProduct {
	
	my ($V1, $V2) = @_;
	
	return [$V1->[1]*$V2->[2] - $V1->[2]*$V2->[1],
	 				$V1->[2]*$V2->[0] - $V1->[0]*$V2->[2], 
					$V1->[0]*$V2->[1] - $V1->[1]*$V2->[0] ];
	
}

sub DotProduct {
  my ($Vector1, $Vector2) = @_;
  return $Vector1->[0]*$Vector2->[0] + $Vector1->[1]*$Vector2->[1] + $Vector1->[2]*$Vector2->[2];
}

sub ArcCosine { 
	my $Num = shift;
	return atan2( sqrt(1 - $Num * $Num), $Num ); 
}

sub CaculateDistanceBetweenAtoms {
	
	my ($Atom1, $Atom2) = @_;
	
	my ($Coords1, $Coords2) = ($Atom1->getCartesianCoordinates, $Atom2->getCartesianCoordinates);
	
	#my $Distance = sqrt ( $)
	
	
} 

sub MoveStructure {
	
  my ($Structure, $x, $y, $z) = @_;

  foreach my $atom (@{$Structure->getAtoms}) {
	
    $atom->setX($atom->getX + $x);
    $atom->setY($atom->getY + $y);
    $atom->setZ($atom->getZ + $z);

  } 

}

sub RotateStructure {
	
  my ($Structure, $axis, $a) = @_;
  my @R = GenerateRotationMatrix($axis, $a);
  foreach my $atom (@{$Structure->getAtoms}) { 
    my $x = $R[0][0]*$atom->getX + $R[0][1]*$atom->getY + $R[0][2]*$atom->getZ;
    my $y = $R[1][0]*$atom->getX + $R[1][1]*$atom->getY + $R[1][2]*$atom->getZ;
    my $z = $R[2][0]*$atom->getX + $R[2][1]*$atom->getY + $R[2][2]*$atom->getZ;
    $atom->setX($x);
    $atom->setY($y);
    $atom->setZ($z);
  }
}

sub GenerateRotationMatrix {
  my ($axis, $a) = @_;
  my @u = @{$axis};
  my $length = 0;
  foreach my $v (@u) { $length += $v**2  }
  $length = sqrt($length);
  foreach my $v (@u) { $v = $length == 0 ? 0 : $v /$length  } 

  my @matrix = ([$u[0]**2 + (1 - $u[0]**2)*cos($a), $u[0]*$u[1]*(1-cos($a))- $u[2]*sin($a), $u[0]*$u[2]*(1 - cos($a)) + $u[1]*sin($a)], 
                [$u[0]*$u[1]*(1-cos($a)) + $u[2]*sin($a), $u[1]**2 + (1-$u[1]**2)*cos($a), $u[1]*$u[2]*(1-cos($a)) - $u[0]*sin($a)], 
                [$u[0]*$u[2]*(1-cos($a)) - $u[1]*sin($a), $u[1]*$u[2]*(1-cos($a)) + $u[0]*sin($a), $u[2]**2 + (1 - $u[2]**2)*cos($a)]);

  return @matrix;

}

return 1;