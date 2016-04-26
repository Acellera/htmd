#!/usr/bin/perl -w


# MJH 5/2/2015
# This program reads the combined rtf+prm str file from GAMMP
# 
 
use strict;
use File::Slurp;
use Data::Dumper;
use Math::Round qw(nearest);
#use File::Which;;

my $mapping={};
my $mapping_idx=97;

#foreach my $dep ("babel" ) { #, "tleap", "psfgen" ) {
#	if( !which($dep) ) { die "This program depends on $dep"; }
#}

sub r{
	my $y = shift;
	my $x= nearest( 0.001, $y );
	return $x;
}

sub t{
	my $x= shift;
	if( $x eq "X" ) { return "X "; }
	$x = lc($x); 
	if( !defined( ${$mapping}{$x} ) ) {
		${$mapping}{$x} = "t" . chr($mapping_idx);
		$mapping_idx++;
	}
	return ${$mapping}{$x};
}

my @f= read_file( "mol.rtf" ) or die "Cannot read mol.rtf";

my $mass= {};
my $residues = {};
my $curres;
my $end_of_rtf=0;
my $i =0;
for( $i = 0 ; $i<scalar(@f) && !$end_of_rtf; $i++ ) {
	my @tok = split(/\s+/, $f[$i] );
	if( scalar(@tok) ) {
	if( $tok[0] eq "MASS" ) {
		${$mass}{ t($tok[2]) } = $tok[3];
	}
	if( $tok[0] eq "RESI" ) {
		$curres={};
		${$residues}{$tok[1]} = $curres;
		${$curres}{charge} = $tok[2]; 
	}
	if( $tok[0] eq "ATOM" ) {
		my $f = { "type" => t($tok[2]), "charge" => $tok[3] };
		${$curres}{"atom"}{$tok[1]} = $f;
	}
	if( $tok[0] eq "BOND" ) {
		${$curres}{"bond"}{"$tok[1]-$tok[2]"} = $tok[5];
	}
	if( $tok[0] eq "IMPH" ) {
		${$curres}{"bond"}{"$tok[1]-$tok[2]-$tok[3]-$tok[4]"} = 1;
	}
	if( $tok[0] eq "END" ) {
		$end_of_rtf = 1;
	}
	}
}

@f= read_file( "mol.prm" ) or die "Cannot read mol.prm";
$i=0;

my $sec=0;
my $e14fac = 1.0;

my $bonds     = {};
my $angles    = {};
my $dihedrals = {};
my $imprs     = {};
my $nb        = {};

my $done=0;
for( ; $i < scalar(@f) && !$done; $i++ ) {
	if( scalar($f[$i]) && ! ($f[$i]=~ /\s*!/)  && ! ($f[$i] =~ /\*/) && !($f[$i] =~ /^read/) ) {
		my @tok=split(/\s+/, $f[$i] );
		if( scalar(@tok) > 0 ) {
print "$f[$i]";
		if(    $tok[0] eq "BONDS" )     { $sec="B"; }
		elsif( $tok[0] eq "END" )    { $done=1; }   
		elsif( $tok[0] eq "ANGLES" )    { $sec="A"; }   
		elsif( $tok[0] eq "DIHEDRALS" ) { $sec="D"; }   
		elsif( $tok[0] =~ /IMPROPER/ )  { $sec="I";  }   
		elsif( $tok[0] eq "NONBONDED" ) { $sec="N"; if ($f[$i] =~ /-$/) { $i++;} } 
		else {
			if( $sec eq "B" ) {
				${$bonds}{ t($tok[0])."-".t($tok[1]) } = { K => r($tok[2]), r=> r($tok[3]) };
			}
			elsif( $sec eq "A" ) {
				${$angles}{ t($tok[0])."-".t($tok[1])."-".t($tok[2]) } = { K => r($tok[3]), theta=> r($tok[4]) };
			}
			elsif( $sec eq "D" ) {
				${$dihedrals}{ t($tok[0])."-".t($tok[1])."-".t($tok[2])."-".t($tok[3]) } = { K => r($tok[4]), n=> r($tok[5]), delta=>r($tok[6]) };
			}
			elsif( $sec eq "I" ) {
				${$imprs}{ t($tok[0])."-".t($tok[1])."-".t($tok[2])."-".t($tok[3]) } = { K => r($tok[4]), n=> r($tok[5]), delta=>r($tok[6]) };
			}
			elsif( $sec eq "N" ) {
				${$nb}{$tok[0]} = {
					epsilon => r($tok[2]),
					rmin => r($tok[3])
				};
				if( defined($tok[5]) ) { ${$nb}{epsilon14} = r($tok[5]); }
				if( defined($tok[6]) ) { ${$nb}{rmin14}    = r($tok[6]); }
			}
			else { die "Unknown line //$f[$i]//"; }
		}
	}     
	}
}

# Output the frcmod
open ( my $fh, ">", "mol.frcmod")  or die "Cannot open file";


print $fh "Generated from GAMMP output\n";
print $fh "MASS\n";
foreach my $k (sort(keys %{$mass})) {
	print $fh lc($k)." \t". ${$mass}{$k} . " \t0.000\n"; # Don't have pols.
}


print $fh "\nBOND\n";
foreach my $k (sort(keys %{$bonds})) {
	print $fh lc($k)." \t". ${$bonds}{$k}{K}. " \t" . ${$bonds}{$k}{r} ."\n";
}

print $fh "\nANGLE\n";
foreach my $k (sort(keys %{$angles})) {
	print $fh lc($k)." \t". ${$angles}{$k}{K}. " \t" . ${$angles}{$k}{theta} ."\n";
}

print $fh "\nDIHE\n";
foreach my $k (sort(keys %{$dihedrals})) {
	my $x = lc($k);
	$x =~ s/x/X/g;
	print $fh $x." \t1.000\t". ${$dihedrals}{$k}{K}. " \t" . ${$dihedrals}{$k}{delta} ." \t" . ${$dihedrals}{$k}{n}."\n";
}

print $fh "\nIMPROPER\n";
foreach my $k (sort(keys %{$imprs})) {
	my $x = lc($k);
	$x =~ s/x/X/g;
	print $fh $x." \t". ${$imprs}{$k}{K}. " \t" . ${$imprs}{$k}{delta} ." \t" . ${$imprs}{$k}{n}."\n";
}

print $fh "\nNONBON\n";
foreach my $k(sort(keys %{$nb})) {
	print $fh t($k)." \t" . ${$nb}{$k}{rmin}. " \t" . ${$nb}{$k}{epsilon}."\n";
}

close ($fh);



## Next read in the mol2 and PDB file, and correct the charges
# The PDB contains the optimised geometry
# The mol2 contains the bonding information

my @pdb = read_file( "QM-min.pdb" ) or die "Cannot open QM-min.pdb";
my @atom=();
for( $i=0; $i<scalar(@pdb); $i++ ) {
	if( $pdb[$i] =~ /^ATOM/  || $pdb[$i]=~ /^HETATM/) {
		my @a=split(/\s+/, $pdb[$i] );
		push(@atom, \@a);
	}
}

my $bonds_from_mol=`babel -i pdb QM-min.pdb -o mol2 - 2> /dev/null | awk 'BEGIN{a=0}{if(\$1=="@<TRIPOS>BOND"){a=1} if(a==1){print \$0}  }'`;
if( ! ($bonds_from_mol  =~ /^@<TRIPOS>/ ) ) { die "CAnnot make mol2 from mol-qm.pdb"; }

open ( $fh, ">", "QM-min.mol2")  or die "Cannot open mol.mol2";
print( $fh "@<TRIPOS>MOLECULE\n");
print( $fh "MOL\n");
print( $fh " ".scalar(@atom)." 0 0 0 0\n");
print( $fh "SMALL\n");
print( $fh "USER_CHARGES\n\n");
print( $fh "@<TRIPOS>ATOM\n");
for( $i=0; $i<scalar(@atom); $i++ ) {
	my $r = ${$residues}{$atom[$i][3]};
	if( !defined($r) ) { die "Unknown residue $atom[$i][3]" }
	my $at = ${$r}{atom}{$atom[$i][2]};
	printf( $fh "%5s %2s       %9.4f %9.4f %9.4f %2s      1  MOL1     %9.4f\n", $i, $atom[$i][2], $atom[$i][6], $atom[$i][7], $atom[$i][8], ${$at}{type}, ${$at}{charge} );
}
print $fh $bonds_from_mol;
close($fh);


## Output tleap input

#open ($fh, ">", "build_tleap.in" ) or die "Cannot open build_tleap.in";
#print( $fh "loadAmberParams mol.frcmod\n" );
#print( $fh "MOL = loadMol2 mol.mol2\n" );
#print( $fh "saveAmberParm MOL structure.prmtop structure.crd\n" );
#print( $fh "quit\n" );
#close( $fh );
#
#`tleap -f build_tleap.in > build_tleap.log 2>&1`;
#
#if( ! -e "structure.prmtop" ) { die "Did not create structure.prmtop - check build_tleap.log"; }


## Make psfgen input

#open ( $fh, ">", "build_psfgen.in" ) or die "Cannot open build_psfgen.in";
#print( $fh "topology ff.str\n" ); 
#print( $fh "segment A {\n\tpdb mol-qm.pdb\n\tfirst none\n\tlast none\n}\n" );
#print( $fh "coordpdb mol-qm.pdb A\n");
#print( $fh "regenerate angles dihedrals\n");
#print( $fh "writepsf structure.psf\n" );
#print( $fh "writepdb structure.pdb\n" );
#close($fh);

#`psfgen build_psfgen.in < /dev/null > build_psfgen.log 2>&1`;

#if( ! -e "structure.psf" ) { die "Did not create structure.psf - check build_psfgen.log"; }
exit 0;





