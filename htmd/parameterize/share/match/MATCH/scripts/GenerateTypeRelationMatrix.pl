#!/usr/bin/env perl

use strict;
use warnings;
use Carp;

#Get MATCH libraries 
use lib $ENV{'MATCH'} . "/lib"; 

use MATCHBaseObject;
use BaseObject ':vars';

use MATCHParameters;
use MoleculeFileHandler;

use MATCHer;

#Setup Default Parameters

my $Forcefield = $ARGV[0];

my $DefaultParameters = MATCHParameters->New;

$DefaultParameters->Initiate("../resources/" . $Forcefield . ".par");

$Parameters = $DefaultParameters;

my $TypesFilePath = $Parameters->getTypesFilePath;

open(FILE, $TypesFilePath);

my @TypesFileContents = grep { $_ !~ /^(?:\#| )/} <FILE>;

close(FILE);

my @Types;

my %TypeHash; 

foreach my $Line (@TypesFileContents) {
			
  my @ArrayedLine = split (/\s+/, $Line);

  next if @ArrayedLine != 3;

  my $Type = Type->New(@ArrayedLine);

  $Type->Initiate;

  $TypeHash{$Type->getName} = [];

  push @Types, $Type;
	
}

foreach my $Type (@Types) {
	
	push @{$TypeHash{$Type->getName}}, $Type;
	
}

foreach my $TypeName1 (keys %TypeHash)  {
	
	my @List;
	
	print $TypeName1 . " " . $TypeName1 . " 1 ";
	
	foreach my $TypeName2 (keys %TypeHash) { 
		
	  next if $TypeName1 eq $TypeName2;
	
	  next if substr($TypeName1,0,1) ne substr($TypeName2, 0,1);
		
		push @List, [$TypeName2, CalculateRelated($TypeName1,$TypeName2)];
		
	}
	
	my @SortedList = sort { $b->[1] <=> $a->[1] } @List;
	
	foreach my $Element (@SortedList) {
		
		print $Element->[0] . " " . $Element->[1] . " ";
		
	}
	
	print "\n";
	
	#while (my ($TypeName2, $Types2) = each %TypeHash) {
		

		
	#}
	
}

sub CalculateRelated {
	
	my ($Type1, $Type2) = @_;
	
	my $TypeObjects1 = $TypeHash{$Type1};
	my $TypeObjects2 = $TypeHash{$Type2};
	
	my $Score = 0;
	my $Count = 0;
	
	#print $Type2 . "\n";
	
	foreach my $TypeObject1 (@$TypeObjects1) {
		
		foreach my $TypeObject2 (@$TypeObjects2) {
					
		  my $Consensus = $TypeObject1->getLookUpTable->CompareTypeLookUpTables([$TypeObject2->getLookUpTable]);
			
			#print $Consensus->Stringify . "\n" if $Consensus != 0;
			
			$Score += $Consensus*2 / (length($TypeObject1->getString) + length($TypeObject2->getString));
			$Count++;

	  }
	
  }

  #return 1;
  return $Score / $Count;
	
}


