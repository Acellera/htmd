#!/usr/bin/env perl

use strict;
use warnings;

use Test::More tests => 1;
use lib $ENV{'MATCH'} . "/lib";  
use Carp;
use Storable;

use MATCHBaseObject;
use BaseObject ':vars';

use MATCHParameters;
use MATCHFunctions ':all';
use MoleculeFileHandler;

use AtomTyper;
use AtomTypeSubstituter;
use MoleculeParameterizer;


#Deal With Initial Parameter Setup, This will be handled in main .pl script

my $DefaultParameters = MATCHParameters->New;

$DefaultParameters->Initiate("../resources/" . $ARGV[0] . ".par");
#$DefaultParameters->Initiate;

$Parameters = $DefaultParameters;

my $Chains =  retrieve("StoredTopologys/" . $ARGV[1] . ".dat");

my $AtomTypeSubstituter = AtomTypeSubstituter->New;

$AtomTypeSubstituter->Initiate;

my $AtomTyper = AtomTyper->New;

$AtomTyper->Initiate($AtomTypeSubstituter);

my $MoleculeParameterizer = MoleculeParameterizer->New;

$MoleculeParameterizer->Initiate($AtomTypeSubstituter);

$Verbosity = 2;

my $TotalChargeDiff = 0;
my $Count = 0;


my $OrigParameters = MATCHParameters->New;

$OrigParameters->Initiate("../resources/" . $ARGV[1] . ".par");

my $OrigParameterFile = $OrigParameters->getParameterFilePath;

my $OrigParameterLists = LoadAtomicParametersFromFile($OrigParameterFile);

my $OrigNonbonded = $OrigParameterLists->{'Nonbond'};
my $OrigBond = $OrigParameterLists->{'Bond'};
my $OrigAngle = $OrigParameterLists->{'Angle'};
my $OrigDihedral = $OrigParameterLists->{'Dihedral'};




#open(FILE, ">test.dat");

open(EPSILON, ">nonbond_epsilon.dat");
open(RMIN, ">nonbond_rmin.dat");

open(B0, ">bond_b0.dat");
open(KB, ">bond_kb.dat");

open(KTHETA, ">angle_ktheta.dat");
open(THETA0, ">angle_theta0.dat");

open(N, ">dihedral_n.dat");
open(KCHI, ">dihedral_kchi.dat");

my $DiffInDelta = 0;
my $DiheCount = 0;
my $PercentCorrectN = 0;

foreach my $Chain (@$Chains) {
			
	my @Atoms = map { @{$_->getAtoms } } @{$Chain->getMolecules};
	
	next if @Atoms < 2; #Ions
	
	next if $Chain->getName eq "CHL1" || $Chain->getName eq "CHNS" || $Chain->getName eq "CHM1";

  my $CopiedChain = $Chain->Copy;
	
	next if $AtomTyper->TypeAtomsInChain($Chain) == -1; 

  foreach my $i (0 .. @{$Chain->getMolecules}-1) {
	
	  my @RequiredParameters = $MoleculeParameterizer->GenerateAllPossibleParametersForMolecule($Chain->getMolecule($i));
	
	  my %AtomNameHash = map { $_->getName => $_->getType } @{$CopiedChain->getMolecule($i)->getAtoms};
	
	  $Chain->getMolecule($i)->setImpropers([]);
	
	 	$MoleculeParameterizer->NewParameterizeMolecule($Chain->getMolecule($i));
	
	  print $Chain->getName . "\n";
	
		my $ParameterFilePath = $ARGV[1] . ".prm";

		my $ParameterLists = LoadAtomicParametersFromFile($ParameterFilePath);

		my $Nonbonded = $ParameterLists->{'Nonbond'};
		my $Bond = $ParameterLists->{'Bond'};
		my $Angle = $ParameterLists->{'Angle'};
		my $Dihedral = $ParameterLists->{'Dihedral'};
		
		
		my @CoupledParameters = ([$Bond, $OrigBond], [$Angle, $OrigAngle], [$Dihedral, $OrigDihedral]);
	
	 
	  foreach my $j (0 .. @RequiredParameters-1) {
		
		  if($j == 0) {
			
		    foreach my $Atom (@{$RequiredParameters[$j]}) {
						
			    my $NewParameter; my $OrigParameter;
			
			    foreach my $Parameter (@$Nonbonded) {
				
				    if($Parameter->{'Type'} eq $Atom->getType) {
					
					    $NewParameter = $Parameter; last;
					
			     	}
				
			    }
			
			    foreach my $Parameter (@$OrigNonbonded) {
				
						if($Parameter->{'Type'} eq $AtomNameHash{$Atom->getName}) {
					
					    $OrigParameter = $Parameter; last;
					
			     	}
				
			    }
			
			    #print $NewParameter->{'Type'} . " " . $OrigParameter->{'Type'} . "\n";
					
					print EPSILON $OrigParameter->{'Epsilon'} . " " . $NewParameter->{'Epsilon'} . "\n";
					print RMIN $OrigParameter->{'Rmin'} . " " . $NewParameter->{'Rmin'} . "\n";
					
					
		    }
		
		  }
		
		  elsif($j < 3) {
						
			  my @CurrentCoupledParameters = @{shift @CoupledParameters};
			
			  my $UniqueParameters = ExtractUniqueParameters($RequiredParameters[$j]);
			
			  foreach my $AtomArray (@$UniqueParameters) {
								
				  my $NewParameter; my $OrigParameter;
							
			    my $TypeString = join(" ", map { $_->getType } @$AtomArray);
			
			    my $ReverseTypeString =  join(" ", reverse map { $_->getType } @$AtomArray);
			
			    my $OrigTypeString = join(" ", map { $AtomNameHash{$_->getName} } @$AtomArray);
				
				  my $OrigReverseTypeString = join(" ", reverse map { $AtomNameHash{$_->getName} } @$AtomArray);	
								
				  foreach my $Parameter (@{$CurrentCoupledParameters[0]}) {
				
				    my $ParameterTypeString = join(" ", @{$Parameter->{'Types'}});
																		
					  if($TypeString eq $ParameterTypeString || $ReverseTypeString eq $ParameterTypeString) {
						
						  $NewParameter = $Parameter; last;
						
					  }
					
			   	}
			
			    foreach my $Parameter (@{$CurrentCoupledParameters[1]}) {

					  my $ParameterTypeString = join(" ", @{$Parameter->{'Types'}});
										
						if($OrigTypeString eq $ParameterTypeString || $OrigReverseTypeString eq $ParameterTypeString) {

							$OrigParameter = $Parameter; last;

						}

				  }
					
				  if($j == 1) {
														
				    print B0 $OrigParameter->{'b0'} . " " . $NewParameter->{'b0'} . "\n";
				    print KB $OrigParameter->{'Kb'} . " " . $NewParameter->{'Kb'} . "\n";
				    
				
			    }
			
			    elsif($j == 2) {
				
				    print KTHETA $OrigParameter->{'Ktheta'} . " " . $NewParameter->{'Ktheta'} . "\n";
						print THETA0 $OrigParameter->{'Theta0'} . " " . $NewParameter->{'Theta0'} . "\n";
				
			
			    }
				
				
			  }			  
			
			
		  }
		
		  else {
			
			  my @CurrentCoupledParameters = @{shift @CoupledParameters};

			  my $UniqueParameters = ExtractUniqueParameters($RequiredParameters[$j]);
			
			  foreach my $AtomArray (@$UniqueParameters) {
								
				  my $NewParameter; my $OrigParameter;
							
			    my $TypeString = join(" ", map { $_->getType } @$AtomArray);
			
			    my $ReverseTypeString =  join(" ", reverse map { $_->getType } @$AtomArray);
			
			    my $OrigTypeString = join( " ", map { $AtomNameHash{$_->getName} } @$AtomArray);
				
				  my $OrigReverseTypeString = join( " ", reverse map { $AtomNameHash{$_->getName} } @$AtomArray);
			
					my @NewParameters; my @OrigParameters;
			 
			    foreach my $Parameter (@{$CurrentCoupledParameters[0]}) {
				
				    my $ParameterTypeString = join(" ", @{$Parameter->{'Types'}});
												
					  if($TypeString eq $ParameterTypeString || $ReverseTypeString eq $ParameterTypeString) {
						
						  push @NewParameters,$Parameter;
						
					  }
					
			   	}
			
			    foreach my $Parameter (@{$CurrentCoupledParameters[1]}) {

					  my $ParameterTypeString = join(" ", @{$Parameter->{'Types'}});
										
						if($OrigTypeString eq $ParameterTypeString || $ParameterTypeString eq $OrigReverseTypeString) {

							push @OrigParameters, $Parameter;

						}

			   	}
			
			    print N scalar(@OrigParameters) . " " . scalar(@NewParameters) . "\n"; 
			
			    if(@OrigParameters == @NewParameters) {
				
				    $DiheCount += @OrigParameters;

						my %Seen;

						my @Pairs;
				
				    foreach my $OriginalDihe (@OrigParameters) {

					 	  foreach my $BestDihe (@NewParameters) {

							  if($OriginalDihe->{'n'} eq $BestDihe->{'n'} && ! exists $Seen{$BestDihe}) {

								  print KCHI $OriginalDihe->{'Kchi'} . " " .  $BestDihe->{'Kchi'} . "\n";

									$DiffInDelta++ if $OriginalDihe->{'delta'} eq $BestDihe->{'delta'};

									push @Pairs, [$OriginalDihe, $BestDihe];

									$Seen{$BestDihe} = 1;

								}

							}

						}

						$PercentCorrectN += @Pairs;
				
				
			    }
			
			  }
			
		  }
		  
				
	  }
	

  }
	
  system("rm " . $ARGV[1] . ".prm");

}

close(EPSILON);
close(RMIN);
close(B0);
close(KB);
close(N);

print $PercentCorrectN / $DiheCount . " " . $DiffInDelta / $DiheCount . "\n";

sub ExtractUniqueParameters {
	
	my $Parameters = shift;
	
	my %Seen; 
	
	my @UniqueParameters;
	
	foreach my $AtomArray (@$Parameters) {
		
		my @Parameters =   map { $_->getType } @$AtomArray;
		
		push @UniqueParameters, $AtomArray if ! $Seen{join(" ", @Parameters)} && ! $Seen{join(" ", reverse @Parameters)};
		
		$Seen{join(" ", @Parameters)} = 1;
		
	}
		
	return \@UniqueParameters;
	
}

#close(FILE);