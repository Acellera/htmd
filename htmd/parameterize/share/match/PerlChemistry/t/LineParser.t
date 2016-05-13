#!/usr/bin/perl -w

use strict;
use warnings;

use lib '../lib';
use Test::More tests => 4;

BEGIN { use_ok('LineParser') };

use LineParser ':all';

#Test Parsing lines from PDB files, Each LineParser can only parse a certain type of line
#For PDBs there are 3
#Atom: Parses lines that declare atoms
#Bond: Parses CONECT lines for bonding information
#End:  Parsers when the Complex is finished, if it is an NMR structure and has more then one complex per file
subtest "PDB Line Parse Tests" => sub {
	
	plan tests => 3;

  #Test the functionality of the Atom line parser, keep in mind this is just a test that what is in the file is collected correctly
  #there is nothing here to make sure that the Atom objects are correctly built
  subtest "Atom Representation Test" => sub {

    #Regular Expression string that will be used to distinguish if a line is an atom declartion
		my $Atom = qr/^(?:ATOM|HETATM)/;
		
		#Variables that will be collected from the atom line
	  my @AtomVariables = qw(Num Name ResidueName ChainName ResidueNum X Y Z Occupancy TempFactor SegmentId);
	
	  #The length of characters that each variable uses, this information is universal for PDB files
	  #For example the Number of an atom occupies the first 5 characters of an Atom line declaration in a PDB file 
		my @AtomSubStrs = (5, 6, 4, 1, 6, 10, 8, 8, 7, 6, -1);
	  
    #Declartion of the Atom line parser, see LineParsers.pm to understand more about how LineParsers work
    my $AtomLineParser =	    LineParser->New('Atom',
																							
																							#Parsing Function
																				    	sub { 
																					
																					      #If doesn't match Atom return nothing since no variables were collected
																								if($_[0] !~ $Atom) { return ( ) }

																								#Remove ATOM or HETATM from beginning of line
																								my $string = substr($_[0], 6); 

																								#Using the length specified in AtomSubStrs, collect the substrings for there respective variables
																								return map { my $substr = substr($string, 0, $_);
																														 $string = substr($string, $_);
																														 $substr } @AtomSubStrs 
																						  }, 
																						
																							#ReturnResults Function
																					  	sub {
																						 
																						    my $Hash;
																						
																								foreach (0 .. $#AtomVariables) {
																									
																									#Remove Whitespace from both ends of Variable string
																							    $_[0]->[$_] =~ s/^\s*(.*?)\s*$/$1/;
																									
																									#Looks crazy but just assign the value of a given variable so it can be looked up by name
																									# Hash->{Name} will return the name of the atom, etc																									
																									$Hash->{$AtomVariables[$_]} = $_[0]->[$_]  if defined $_[0]->[$_];
																								}
																								
																								return $Hash  
																						});

    #Go through examples in File
    TestParsingExamples($AtomLineParser, $ENV{'PerlChemistry'} . "t/LineParserResources/AtomPdbLines.txt");
    
	};

	#Test the functionality of the Bond line parser
	subtest "Bond Representation Test" => sub {

		plan tests => 3;
		
		my $BondPattern = qr/^CONECT\s/;

		my $LineParser = LineParser->New('Bond',
																		 sub { if($_[0] =~ $BondPattern) { my @spl = split /\s+/, $_[0]; shift @spl; @spl }},
																		 sub { return { StartAtomNum => shift @{$_[0]}, BondAtomNums => $_[0] }});

	  $LineParser->Parse("CONECT 1 2 3 4");

	  my $DataHash = $LineParser->ReturnResults;

	  is($DataHash->{'StartAtomNum'}, 1, 'Worked');

	  is_deeply($DataHash->{'BondAtomNums'}, [2,3,4], 'Worked');

	  is($LineParser->Parse("END"), 0, 'Not Matched');

	};
	
	#Test the functionality of the End line parser
	subtest "End Representation Test" => sub {

		plan tests => 3;
		
		my $EndPattern = qr/^(TER)/;
		my $EndComplexPattern = qr/^(ENDMDL)/;

		my $LineParser = LineParser->New('End', sub { $_[0] =~ $EndPattern }, sub { } );

		is($LineParser->Parse("HETATM  869  C11 PAR A  28       7.827   5.531  25.866  1.00 20.00           C "),
			 0,
			 'Not Matched');

		is($LineParser->Parse("TER  890"), 1, 'Matched Correctly');

		is($LineParser->ReturnResults, undef, 'Returns Nothing');

	};
	
};

#Test Parsing lines for MOL2 files, Each LineParser can only parse a certain type of line
#For MOL2 there are 3
#Atom: Parses lines that declare atoms
#Bond: Parses lines that declare bonds
#End:
subtest "MOL2 Line Parse Tests" => sub {
	
	plan tests => 2;
	
	subtest "Atom Representation Test" => sub {
		
		#Regular Expression string that will be used to distinguish if a line is an atom declartion		
	  my $AtomWithCharge = qr/^\s*(\d+)\s+(\S+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+\S+\s+(\d+)\s+(\S+)\s+(-?\d+\.\d+)/;
	  
	  #Variables that will be collected from the atom line
	  my @AtomWithChargeVariables = qw(Num Name X Y Z ResidueNum ResidueName Charge);
		
		my $LineParser = LineParser->New('Atom',
																		 #Parse Function
																		 sub { $_[0] =~ $AtomWithCharge },
																		 #ReturnResults Function
																		 sub { return { map {$AtomWithChargeVariables[$_] => $_[0]->[$_] } (0 .. $#AtomWithChargeVariables) } } );
																		
		#Test Examples
	  TestParsingExamples($LineParser, $ENV{'PerlChemistry'} . "t/LineParserResources/AtomMol2Lines.txt");
				
	};

  subtest "Bond Representation Test" => sub {
	
	  plan tests => 4;
	
	  my $Bond = qr/^\s*(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s*\n/;
		my @BondVariables = qw(Num BondAtomNums Type);
	
	  my $LineParser = LineParser->New('Bond',
    																 sub { $_[0] =~ $Bond; return ($1, [$2, $3], $4) },
																		 sub { return { map {$BondVariables[$_] => $_[0]->[$_] } (0 .. $#BondVariables) } } );

	  is($LineParser->Parse("	8     7    1     0    0\n"),
			 0,
			 'Does not Match Numbers under name');

    is($LineParser->Parse("   1 C1       1.0848     0.1201    -0.0975 C.3      1  TMP  0.0102"),
			 0,
			 'Does not Match Atoms');

    is($LineParser->Parse("	 1    1    2 1  \n"),
 			 1,
			 'Matched');
			
	  my $DataHash = $LineParser->ReturnResults;
		
		is_deeply($DataHash->{'BondAtomNums'}, [1,2], 'BondAtomNums correct')
		
	
  };

	
};

#Test Parsing Lines for RTF files, Each LineParser can only parse a certain type of line
#Atom: Parses lines that declare atoms
#Bond:
#Improper:
#Delete: 
subtest "RTF Line Parse Tests" => sub {

	plan tests => 4;
	
	subtest "Atom Representation Test" => sub {
		
    my $Atom = qr/^\s*ATOM\s+(\S+)\s+(\S+)\s+(-?\d*\.\d+)/;
    my @AtomVariables = qw( Name Type Charge);
  
    my $LineParser = LineParser->New('Atom',
																		 sub { uc($_[0]) =~ $Atom },
																		 sub { return { map {$AtomVariables[$_] => $_[0]->[$_] } (0 .. $#AtomVariables) } } );
	
		
		TestParsingExamples($LineParser, $ENV{'PerlChemistry'} . "t/LineParserResources/AtomRtfLines.txt");

  };

	subtest "Bond Representation Test" => sub {
		
		plan tests => 6;
		
		my $Bond = qr/^\s*BOND\s+([^!]+)/;
		my @BondVariables = qw(BondAtomNames Type);
		
		my $BondLineParser =    LineParser->New('Bond',
																				    sub { 
																				      return () if uc($_[0]) !~ $Bond;
																				      my @spl = split /\s+/, $1;
	                                            return map { [ [$spl[$_], $spl[$_ + 1]], 1] } grep { $_ % 2 == 0 } (0 .. $#spl);																																   
																				    },
																			    	sub { 
																					    my @Hashes;
																					    foreach my $Bond (@{$_[0]}) {
																					      push @Hashes, { map {$BondVariables[$_] => $Bond->[$_] } (0 .. $#BondVariables) }
																					    }
																					    return \@Hashes;
																				    } );
																		
		is($BondLineParser->Parse("BOND O  C   CG  CD2   CE1  NE2"), 1, 'Matched');
		
		my $DataHash = $BondLineParser->ReturnResults;
		
		is(@$DataHash, 3, "Collected the correct number of bonds");
		
		is_deeply($DataHash->[0]->{'BondAtomNames'}, ['O', 'C'], 'Parsed Correctly');
		
		#Test if comments break it
		is($BondLineParser->Parse("BOND O  C   CG  CD2   CE1  NE2 ! TEST COMMENT"), 1, 'Matched 2');
		
		is(@$DataHash, 3, "Collected the correct number of bonds");
		
		is_deeply($DataHash->[0]->{'BondAtomNames'}, ['O', 'C'], 'Parsed Correctly');
		
						
	};
	
  subtest "Improper Representation Test" => sub {
	
	  plan tests => 2;
		
	  my $Improper = qr/^\s*IMPR\s+([^!]+)/;

		my $LineParser = LineParser->New('Improper',
																		 sub { 
																		   return () if $_[0] !~ $Improper;
																																				
																		   my @spl = split /\s+/, $1;
																		
																		   return map { [ $spl[$_], $spl[$_+1], $spl[$_+2], $spl[$_+3]] } grep { $_ % 4 == 0 } (0 .. $#spl)																																   
																		 },
																		 sub {
																			 my @Hashes;
																																						 
																		   foreach my $Improper (@{$_[0]}) {
																			   push @Hashes, { Atoms => $Improper };
																		   }
																		
																			 return \@Hashes;
																			
																		 } );
																		
	  is($LineParser->Parse("IMPR  N  C  HC HT  N  C  HT  HC"),
			 1,
			 'Matched');
			
		my $DataHash = $LineParser->ReturnResults;

		is_deeply($DataHash->[0]->{'Atoms'}, ['N', 'C', 'HC', 'HT'], 'Atoms match up');
			
		
	
  };

  subtest "Delete Representation Test" => sub {
	
    plan tests => 3;  

    my $Delete = qr/^\s*(?:DELE|DELETE)\s+([^!]+)/;
		my @DeleteVariables = qw(Type Name);
		
		my $LineParser = LineParser->New('Delete', 
																		 sub {
																		   return () if $_[0] !~ $Delete;
																			 my @spl = split /\s+/, $1;
																			 return map { [ $spl[$_], $spl[$_ + 1], 1] } grep { $_ % 2 == 0 } (0 .. $#spl)
																		 },
																	   sub { 
																		   my @Hashes;
																			 foreach my $Delete (@{$_[0]}) {
																			   push @Hashes, { map {$DeleteVariables[$_] => $Delete->[$_] } (0 .. $#DeleteVariables) }
																			 }
																			 return \@Hashes
																		 }	);
																		
	  is($LineParser->Parse("	DELE ATOM N9   "), 1, 'Matched');
	
	  my @DataHashes = @{$LineParser->ReturnResults};
	
	  is($DataHashes[0]->{'Type'}, 'ATOM', 'Type is Atom');
	
	  is($LineParser->Parse("	DELE ATOM N9 ATOM N12 ATOM C39   "), 1, 'Matched');
	  

  };
  
};

sub TestParsingExamples {
	
	my ($Parser, $FileName) = @_;
	
	#A list of examples used to make sure the Atom line parser is behaving properly
  open(FILE, $FileName);

  #Read in lines from file, remove lines comment lines
  my @ExampleFileContents =  grep { $_ !~ /^#/ } <FILE>;

  close(FILE);
	
	for(my $i=0; $i < scalar(@ExampleFileContents); $i+=2) {

    my $Line = $ExampleFileContents[$i];
    my %Values = map { my @spl = split /\=/, $_; $spl[0] => $spl[1] } split /\s+/, $ExampleFileContents[$i+1];

    my $IsParsed = $Parser->Parse($Line);

    is($IsParsed, $Values{'Success'}, "Testing Success should be " . $Values{'Success'});	

    if($IsParsed) {
	
	    my $DataHash = $Parser->ReturnResults;
	
	    TestParsedResults($DataHash, \%Values);
	
    }
   

  }
	
	
}


sub TestParsedResults {
	
	my ($DataHash, $Values) = @_;
	
	while ( my ($ValueName, $Value) = each %$Values) {
	
	  next if $ValueName eq "Success";
	
	  #Easy way to handle spaces, declare them as underscores in txt file and the replace them here
	  $Value =~ s/\_/ / if $Value =~ /\_/;
	
	  #If Value is a number, need to check if values are same numerically
	  if($Value =~ /^-?\d+/) {
	
	    ok($DataHash->{$ValueName} == $Value, "$ValueName is correct: $Value");
	
	  }
	
	  #If Value is not a number use string comparison
	  else {
	
	    ok($DataHash->{$ValueName} eq $Value, "$ValueName is correct: $Value");
	 
    }
	
  }
	
	
}

