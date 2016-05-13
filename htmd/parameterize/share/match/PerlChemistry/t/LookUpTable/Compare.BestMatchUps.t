#!/usr/bin/perl -w

use strict;
use warnings;
use Cwd qw(abs_path);
use Benchmark;

use lib '../../lib';
use Test::More tests => 1;

BEGIN { use_ok('LookUpTable') };

use LookUpTable;
use MoleculeFileHandler ':all';
use Parameters ':all';
use BaseObject ':vars';

srand(time ^ ($$ + ($$ << 15)));

sub TESTDetermineBestMatchUps { 
	
  my ($NextLevel, $ReferenceNodes, $MatchingNodes) = @_;

  #print join(" ", map { ref($_) } (@$ReferenceNodes, @$MatchingNodes) ) . "\n";

  my $SizeOfMatchingNodes = @$MatchingNodes;
	my $SizeOfReferenceNodes = @$ReferenceNodes;

	my $Min = $SizeOfMatchingNodes <= $SizeOfReferenceNodes ? $SizeOfMatchingNodes : $SizeOfReferenceNodes;
	
	$_->UpdateChildScoreGrid foreach (@$ReferenceNodes, @$MatchingNodes);
	
	#foreach my $Node (@$ReferenceNodes, @$MatchingNodes) {
		
	#	$Node->UpdateChildScoreGrid;
		
	#}
	
	my @MatchUps;
	
	foreach my $ReferenceNode (@$ReferenceNodes) {
	  push (@MatchUps, map { [$ReferenceNode, $_, 0] } @$MatchingNodes);
	}	
	
	my $Recurse = sub { 
		
		my $List = shift;
		my $Allowed = shift;
		
		my @NextAllowed = grep {$_ % $SizeOfMatchingNodes != $List->[-1] % $SizeOfMatchingNodes && int($_/$SizeOfMatchingNodes) != int($List->[-1]/$SizeOfMatchingNodes) } @$Allowed;
	 
	  return [$List, \@NextAllowed] if ! @NextAllowed;
	
	  my @NewLists = map { [ @$List, $_] } @NextAllowed;
	
	  return map { [$_, \@NextAllowed] } @NewLists;
		
	};
	
	my %ResultsHash;

	foreach my $Start (0 .. $#MatchUps) {

		my @Results = &$Recurse([$Start], [0 .. $#MatchUps]);

		my $Count = 1;

		while($Count < $Min) {
			@Results = map { &$Recurse(@{$_}) } @Results;
			$Count++;
		}

		foreach (@Results) {
			
			my $StringedResult = join(" ", sort @{$_->[0]});
			
		  $ResultsHash{$StringedResult} = $_->[0] if ! exists $ResultsHash{$StringedResult};
	  }
	
	}
	
	my @Permutations = map { [[map { \$_ } @MatchUps[@{$_}]], 0] } values %ResultsHash;
	
	my $PermutationNumber = 0;
	my $BestScore = 0;
	my $Score = 0;
	
	$PermutationNumber = $#Permutations;
	
  
  while( $PermutationNumber > 0) {	
				
		$_->[2] = TESTScoreChildrenUsingStateGrid($_->[0], $_->[1], $NextLevel) foreach (@MatchUps);
		
		foreach (@MatchUps) {
			
			print $_->[0]->getAtomName . " " . $_->[1]->getAtomName . " " . $_->[2] . "\n";
			
		}
		  		
		foreach my $Permutation (@Permutations) {
			$Score = 0;
			$Score += $$_->[2] foreach (@{$Permutation->[0]});
			$BestScore = $Score if $BestScore < $Score;
			$Permutation->[1] = $Score;
		}
		
		my $PermutationPointer = 0;
		
		while ( $PermutationPointer <= $PermutationNumber) {
			
		  if($Permutations[$PermutationPointer]->[1] < $BestScore) { $PermutationNumber--; splice(@Permutations, $PermutationPointer, 1) }
		  else																						         { $PermutationPointer++ }
			
		}
		
		$NextLevel++;
				
		$PermutationNumber = 0 if $NextLevel >= 10 && $PermutationNumber != 0;
				 			
	}
		
	return map { [$$_->[0]->getAtomUniqueString, $$_->[1]] } @{$Permutations[0]->[0]}
	
}


sub TESTScoreChildrenUsingStateGrid($$$) {
	
	my ($Node1, $Node2, $Level) = @_;
	
	my $StateGridForNode1 = $Node1->getChildScoreGrid;
  my $StateGridForNode2 = $Node2->getChildScoreGrid;
	my $Score = 0;
	
	print $Node1->getAtomName . " " . $Node2->getAtomName . "\n";
	
	while(my ($State, $GridElement) = each %$StateGridForNode1) {
	
	  my $StateLevel = substr($State, -2);
	
	  next if $StateLevel > $Level || ! exists $StateGridForNode2->{$State};	
	
	  print $State . " " . $GridElement->getCount . " " . $StateGridForNode2->{$State}->getCount . "\n";
	
		foreach my $NodeReference (@{$GridElement->getNodeReferences }) {
			print $$NodeReference->getAtomName . " "; 
		}
		
		print "\n";
		
		foreach my $NodeReference (@{$StateGridForNode2->{$State}->getNodeReferences }) {
			print $$NodeReference->getAtomName . " "; 
		}
		
		print "\n";
		
    $Score += $GridElement->getCount <= $StateGridForNode2->{$State}->getCount ? $GridElement->getCount : $StateGridForNode2->{$State}->getCount;	
	
	}
	
	return $Score;
	
}


#Deal With Initial Parameter Setup, This will be handled in main .pl script
my $DefaultParameters = Parameters->New;

my ($Path) = (abs_path($0) =~ /(\S+)\/\S+$/);

my $DefaultParametersPath = "../../resources/DefaultParameters.par";

$DefaultParameters->Initiate($Path, $DefaultParametersPath);

#Set Global Parameters
$Parameters = $DefaultParameters;

my $MoleculeFileHandler = MoleculeFileHandler->New("../../resources/top_all22_prot/top_all22_prot.inp");

my $Structures = $MoleculeFileHandler->BuildObjectsFromFile;

my @MoleculeChains = grep { !$_->getMolecule(0)->getIsPatch }  @$Structures;

my ($ProChain) = grep { $_->getMolecule(0)->getName eq 'PRO'} @MoleculeChains;

my $CopiedChain = $ProChain->Copy;

$CopiedChain->Initiate;

my ($CopiedNatom) = grep { $_->getType eq 'N' } @{$CopiedChain->getMolecule(0)->getAtoms}; 

my @CopiedNodes = grep { $_->getAtomName eq "CA" || $_->getAtomName eq "CD"} @{$CopiedNatom->getLookUpTable->getNodes};

$ProChain->Initiate;

my ($Natom) = grep { $_->getType eq 'N' } @{$ProChain->getMolecule(0)->getAtoms}; 

print "START!!!\n";

my $Consensus = $Natom->getLookUpTable->DetermineConsensusLookUpTable([]);

my @RealNodes = grep { $_->getAtomName eq "CA" || $_->getAtomName eq "CD"} @{$Natom->getLookUpTable->getNodes};

my @ConsensusNodes = grep { $_->getAtomName eq "CA" || $_->getAtomName eq "CD"} @{$Consensus->getNodes};

#foreach ( 0 .. 10 ){

  TESTDetermineBestMatchUps(2, \@RealNodes, \@ConsensusNodes);

#}
