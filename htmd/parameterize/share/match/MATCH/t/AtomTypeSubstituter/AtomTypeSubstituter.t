#!/usr/bin/env perl

use strict;
use warnings;

use Test::More tests => 1;
use lib $ENV{'MATCH'} . "/lib"; 
use Carp;

use MATCHBaseObject;
use MoleculeFileHandler;
use BaseObject ':vars';

use MATCHParameters;
use MATCHFunctions ':func';

use AtomCharger;
use AtomTyper;

BEGIN { use_ok('AtomTypeSubstituter') };

#Deal With Initial Parameter Setup, This will be handled in main .pl script
my $DefaultParameters = MATCHParameters->New;

$DefaultParameters->Initiate("../../resources/" . $ARGV[0] . ".par");

$Parameters = $DefaultParameters;

#Get Increments

my $IncrementFilePath = $Parameters->getIncrementFilePath;

open(FILE, $IncrementFilePath);

my @Increments = map { my @spl = split /\s+/, $_; [@spl] }<FILE>;

close(FILE);


my $AtomTypeSubstituter = AtomTypeSubstituter->New;

$AtomTypeSubstituter->Initiate;

my $AvgDiff = 0;

my $Count = 0;

#Test Increments

print "BOND INCREMENTS\n";

open(FILE, ">IncrementSubstitution.dat");

foreach my $Increment (@Increments) {
			
  my $Best = FindIncrementSubstitution($AtomTypeSubstituter, $Increment, @Increments);	
	
	next unless $Best != 0;
	
	#print $Increment->[0] . " " . $Increment->[1] . " "  . $Best->[0] . " " . $Best->[1] . "\n";
	
	#print join(" ", @$Increment ) . " => " . join(" ", @$Best) . "\n";
	
	$Count++;
	
	$AvgDiff += abs($Increment->[2] - $Best->[2]) + abs($Increment->[3] - $Best->[3]);
	
	print FILE $Increment->[2] . " " . $Best->[2] . "\n";
	print FILE $Increment->[3] . " " . $Best->[3] . "\n";
	
		
	#exit 1;
	
}

#print $AvgDiff / ($Count*2) . " " . $Count .  "\n";

close(FILE);

#exit 1;

print "PARAMETERS\n";

#Bonds
my $AvgDiffb0 = 0;
my $AvgDiffKb = 0;
#
my $AvgDiffKtheta = 0;
my $AvgDiffTheta0 = 0;

$Count = 0;

my $ParameterFile = $Parameters->getParameterFilePath;

my $Parameters = LoadAtomicParametersFromFile($ParameterFile);

my $Bonds = $Parameters->{'Bond'};

open(FILE, ">Bondb0Substition.dat");
open(FILE1, ">BondkbSubstition.dat");

foreach my $Bond (@$Bonds) {
		
	my $Best = FindParameterSubstitution($AtomTypeSubstituter, $Bond, $Bonds);
	
	next unless $Best != 0;
		
	$Count++;
		
	$AvgDiffb0 += abs($Bond->{'b0'} - $Best->{'b0'});
	
	print FILE $Bond->{'b0'} . " " . $Best->{'b0'} . "\n";
	
	$AvgDiffKb += abs($Bond->{'Kb'} - $Best->{'Kb'});
	
	print FILE1 $Bond->{'Kb'} . " " . $Best->{'Kb'} . "\n";
 			
}

close(FILE);
close(FILE1);

print "Bonds: Avg Kb: " . $AvgDiffKb / $Count . " Avg b0: " . $AvgDiffb0 / $Count . " " . $Count .  "\n";

$Count = 0;

my $Angles = $Parameters->{'Angle'};

open(FILE, ">AngleKthetaSubstition.dat");
open(FILE1, ">AngleTheta0Substition.dat");

foreach my $Angle (@$Angles) {
		
	my $Best = FindParameterSubstitution($AtomTypeSubstituter, $Angle, $Angles);
	
	next unless $Best != 0;
	
	$Count++;
		
	$AvgDiffKtheta += abs($Angle->{'Ktheta'} - $Best->{'Ktheta'});
	
	#print join(" " ,@{$Angle->{'Types'}}) . " " . join(" " ,@{$Best->{'Types'}}) . " " . abs($Angle->{'Ktheta'} - $Best->{'Ktheta'}) . " " . abs($Angle->{'Theta0'} - $Best->{'Theta0'}) . "\n";
	
	print FILE $Angle->{'Ktheta'} . " " . $Best->{'Ktheta'} . "\n";
	
	$AvgDiffTheta0 += abs($Angle->{'Theta0'} - $Best->{'Theta0'});
	
	print FILE1 $Angle->{'Theta0'} . " " . $Best->{'Theta0'} . "\n";	
	
}

print "Angles: Avg Ktheta: " . 	$AvgDiffKtheta / $Count . " Avg Theta0: " . $AvgDiffTheta0 / $Count . " " . $Count . "\n";

close(FILE);
close(FILE1);

exit 1;

my $Dihedrals = $Parameters->{'Dihedral'};

$Count = 0;

my $AverageNumDiff = 0;
my $PercentCorrectNum = 0;
my $PercentCorrectN = 0;
my $DiffInDelta = 0;

my %SeenDihedrals;
my $DiheCount = 0;

open(FILE, ">DihedralKchi.dat");

foreach my $Dihedral (@$Dihedrals) {
	
	my $Best = FindParameterSubstitution($AtomTypeSubstituter, $Dihedral, $Dihedrals);
	
	next unless $Best != 0;
	
	my $OrigDihedralString = join(" ", @{$Dihedral->{'Types'}});	
	
	my @SameParameters;
	
	foreach my $Dihedral (@$Dihedrals)  {
		
		my $DihedralString = join(" ", @{$Dihedral->{'Types'}});	
	
	  if($DihedralString eq $OrigDihedralString) { push @SameParameters, $Dihedral }	
	
	}
	
	my $BestDihedralString = join(" ", @{$Best->{'Types'}});	
	
	my @BestParameters;
	
	foreach my $Dihedral (@$Dihedrals)  {
		
		my $DihedralString = join(" ", @{$Dihedral->{'Types'}});	
	
	  if($DihedralString eq $BestDihedralString) { push @BestParameters, $Dihedral }	
	
	}
	
	#print @SameParameters . "\n";
	
  $AverageNumDiff += abs(@SameParameters - @BestParameters);
	
	if(@SameParameters == @BestParameters) {
		
		$DiheCount += @SameParameters;
		
		$PercentCorrectNum++;
		
		my %Seen;
		
		my @Pairs;
		
		foreach my $OriginalDihe (@SameParameters) {
			
			foreach my $BestDihe (@BestParameters) {
				
				if($OriginalDihe->{'n'} eq $BestDihe->{'n'} && ! exists $Seen{$BestDihe}) {
					
					print FILE $OriginalDihe->{'Kchi'} . " " .  $BestDihe->{'Kchi'} . "\n";
					
					$DiffInDelta++ if $OriginalDihe->{'delta'} eq $BestDihe->{'delta'};
					
					push @Pairs, [$OriginalDihe, $BestDihe];
					
					$Seen{$BestDihe} = 1;
					
				}
				
			}
			
		}
		
		$PercentCorrectN += @Pairs;
		
	}
	
	$Count++;
	
}

print $AverageNumDiff / $Count . " " . $PercentCorrectNum / $Count . "\n";
print $PercentCorrectN / $DiheCount . " " . $DiffInDelta / $DiheCount . "\n";

close(FILE);

sub FindParameterSubstitution {
	
	my ($AtomTypeSubstituter, $OriginalParameter, $Parameters) = @_;
	
	my @PossibleParameterSubstitutions;
	
	my $OrginalTypes = $OriginalParameter->{'Types'};
	
	my $OriginalParameterString = join(" ", @{$OriginalParameter->{'Types'}});	

	foreach my $Parameter (@$Parameters) {
		
		next if $OriginalParameter eq $Parameter;	
		
		my $ParameterString = join(" ", @{$Parameter->{'Types'}});
		
		next if $OriginalParameterString eq $ParameterString;
		
		my $Types = $Parameter->{'Types'};
		
		my $ReverseTypes = [reverse @$Types];
				
	  if (ElementsAreTheSame($OrginalTypes, $Types) ) {
		
		  push @PossibleParameterSubstitutions, $Parameter;
		
	  }
		
		elsif (ElementsAreTheSame($OrginalTypes, $ReverseTypes)) {
			
			$Parameter->{'Types'} = $ReverseTypes;
			
		  push @PossibleParameterSubstitutions, $Parameter;
			
		}	
		
		
	}
	
	my ($BestMatch, $BestMatchingParameter) = (0,0);
		
	foreach my $Parameter ( @PossibleParameterSubstitutions ) {
		
	  my $Types = $Parameter->{'Types'};
	
	  my $Score = 0;
	
	  my $WeightArray = $WeightIndex{@$OrginalTypes};
	
	  foreach my $i (0 .. @$OrginalTypes-1) {
		
		  $Score += $AtomTypeSubstituter->NewScoreTypeRelation($OrginalTypes->[$i], $Types->[$i]);
		
	  }
	
	  if($BestMatch < $Score)   { $BestMatch = $Score; $BestMatchingParameter = $Parameter}
		
	}
	
	return $BestMatchingParameter;
	
}

sub ElementsAreTheSame {
	
	my ($Types1, $Types2) = @_;
	
	my $Flag = 1;
	
	foreach my $i (0 .. @$Types1-1) {
		
	  if(substr($Types1->[$i],0,1) ne substr($Types2->[$i],0,1) ) {
		
		  $Flag = 0; last;
		
	  }
		
	}
	
	return $Flag;
	
}

sub FindIncrementSubstitution {
	my ($AtomTypeSubstituter, $OrginalIncrement, @Increments) = @_;
	
	my @PossibleIncrementSubstituions;
	
	foreach my $Increment (@Increments) {
		
		next if $OrginalIncrement eq $Increment;
				
		push @PossibleIncrementSubstituions, $Increment if substr($Increment->[0], 0, 1) eq substr($OrginalIncrement->[0], 0, 1) && substr($Increment->[1], 0, 1) eq substr($OrginalIncrement->[1], 0, 1);
		
		push @PossibleIncrementSubstituions, [$Increment->[1], $Increment->[0], $Increment->[3], $Increment->[2]] if substr($Increment->[0], 0, 1) eq substr($OrginalIncrement->[1], 0, 1) && substr($Increment->[0], 0, 1) eq substr($OrginalIncrement->[1], 0, 1)
	
	}
		
	my ($BestMatch, $BestMatchingIncrement) = (0,0);
	
	foreach my $Increment (@PossibleIncrementSubstituions) {
		
		next if $OrginalIncrement->[0] eq $OrginalIncrement->[1] && $Increment->[0] ne $Increment->[1];
		
    next unless ChargeDifferenceAsIncrement($AtomTypeSubstituter, $OrginalIncrement, $Increment);

    my $FowardSubstitutionScore = $AtomTypeSubstituter->NewScoreTypeRelation($OrginalIncrement->[0], $Increment->[0]) + $AtomTypeSubstituter->NewScoreTypeRelation($OrginalIncrement->[1], $Increment->[1]);

    #$FowardSubstitutionScore = 0 if $AtomTypeSubstituter->NewScoreTypeRelation($OrginalIncrement->[0], $Increment->[0]) == 0 || $AtomTypeSubstituter->NewScoreTypeRelation($OrginalIncrement->[1], $Increment->[1]) == 0;

    my $ReversedSubstitutionScore  = $AtomTypeSubstituter->NewScoreTypeRelation($OrginalIncrement->[1], $Increment->[0]) + $AtomTypeSubstituter->NewScoreTypeRelation($OrginalIncrement->[0], $Increment->[1]);
		
		#$ReversedSubstitutionScore = 0 if $AtomTypeSubstituter->NewScoreTypeRelation($OrginalIncrement->[1], $Increment->[0]) == 0 || $AtomTypeSubstituter->NewScoreTypeRelation($OrginalIncrement->[0], $Increment->[1]) == 0;
    
		
    if($BestMatch < $FowardSubstitutionScore)   { $BestMatch = $FowardSubstitutionScore; $BestMatchingIncrement = $Increment}
    if($BestMatch < $ReversedSubstitutionScore) { $BestMatch = $ReversedSubstitutionScore; $BestMatchingIncrement = [$Increment->[1], $Increment->[0], $Increment->[3], $Increment->[2]]}
		
	}
	
  return $BestMatchingIncrement;
	
}


sub ChargeDifferenceAsIncrement {
	
	my ($Self, $Increment1, $Increment2) = @_;
		
	my $ChargeRelationshipLookUp = $Self->getChargeRelationshipLookUp;
		
	my @ChargeRelationshipsForIncrement1 = map { $ChargeRelationshipLookUp->{ $_ } } ($Increment1->[0], $Increment1->[1]);
	
	my @ChargeRelationshipsForIncrement2 = map { $ChargeRelationshipLookUp->{ $_ } } ($Increment2->[0], $Increment2->[1]);
	
		
	my $Foward1 = $Increment1->[2] + $Increment1->[3];
	my $Backward1 = $Increment1->[3] + $Increment1->[2];
	
  my $Foward2 = $Increment2->[2] + $Increment2->[3];
	my $Backward2 = $Increment2->[3] + $Increment2->[2];
	
	return 1 if $Foward1 == $Foward2 && $Backward1 == $Backward2;
	
	return 0;
	
	#print join(" ", (@ChargeRelationshipsForIncrement1, @ChargeRelationshipsForIncrement2)) . "\n";
			
	my $Flag = 0;
	
	foreach my $i (0 .. @ChargeRelationshipsForIncrement1-1) {
		
		$Flag = 0;
		
		foreach my $ChargeRelationship (@{$ChargeRelationshipsForIncrement1[$i]}) {
			
		  foreach my $ChargeRelationship2 (@{$ChargeRelationshipsForIncrement2[$i]}) {
			
			  if($ChargeRelationship == $ChargeRelationship2) {

				  $Flag = 1; 

					return 1 if $ChargeRelationship != 0;

					last;

			  }
			
		  }
		
		  last if $Flag == 1;
			
		}
						
		return 0 if $Flag == 0;
		
	}
		
	return 1;
	
	
}
