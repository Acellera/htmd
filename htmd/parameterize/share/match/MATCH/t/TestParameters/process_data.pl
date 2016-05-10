#!/usr/bin/perl

use strict;
use warnings;


my @files = qw(nonbond_epsilon nonbond_rmin bond_b0 bond_kb angle_theta0 angle_ktheta dihedral_kchi dihedral_n);

foreach my $filename (@files) {

  open(FILE, "$filename\_partial.dat");

  my @FileContents = <FILE>;

  close(FILE);

  print $filename . " " . join(" ", calculate(@FileContents)) . "\n";

}

sub calculate { 

  my @FileContents = @_;

  my ($xsum, $ysum, $xysum, $xsqsum, $ysqsum);

  my $diff = 0;

  my $percentdiff = 0;

  my $count = 0;

  foreach my $Line (@FileContents) {

    my @spl = split /\s+/, $Line;

    next if @spl < 2;

    next if $spl[0] eq "";
    next if $spl[1] eq "";

    $xsum += $spl[0]; $ysum += $spl[1]; $xysum += $spl[0]*$spl[1]; 

    $xsqsum += $spl[0] ** 2; $ysqsum += $spl[1] ** 2;

    $diff += abs( $spl[0] - $spl[1] );

    $percentdiff += abs($spl[0] - $spl[1]) / abs($spl[0]);

    $count++;

  }
  
#  my $r = ($xysum - ($xsum*$ysum)/$count) / (sqrt(($xsqsum - $xsum**2 / $count) * ($ysqsum - $ysum**2 / $count)));

  my $r =  ($count*$xysum - ($xsum*$ysum)) / (sqrt(($count*$xsqsum - $xsum**2) * ($count*$ysqsum - $ysum**2)));
 $r = $r**2;


 return ($r,$diff/$count, $percentdiff/$count*100, $count); 

}
