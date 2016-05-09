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
use AtomCharger;
use AtomTypeSubstituter;
use Type;

#Deal With Initial Parameter Setup, This will be handled in main .pl script

my $DefaultParameters = MATCHParameters->New;

$DefaultParameters->Initiate("../resources/" . $ARGV[0] . ".par");
#$DefaultParameters->Initiate;

$Parameters = $DefaultParameters;

my $Chains =  retrieve("StoredTopologys/" . $ARGV[0] . ".dat");

my %SeenCenterAtoms;
my %SeenType;

my $Parameters = LoadAtomicParametersFromFile("../resources/top_all22_prot/par_all22_prot.prm");

my @Impropers = map { $_->{'Types'} } @{$Parameters->{'Improper'}};

my @Atoms =  map { @{$_->getAtoms} } map { @{$_->getMolecules} } @$Chains;

open(FILE, "../resources/top_all22_prot/top_all22_prot.impr");

my @ImprContents = <FILE>;

close(FILE);

my @ImproperTypes;

foreach my $Line (@ImprContents) {
			
  my @ArrayedLine = split (/\s+/, $Line);

  next if @ArrayedLine != 3;

  my $Type = Type->New(@ArrayedLine);

  $Type->Initiate;

 # print $Type->getName . " " . $Type->getBondNumber . " " . $Type->getString . "\n";

  push @ImproperTypes, $Type;
	
}

@ImproperTypes = sort { length($b->getString) <=> length($a->getString) } @ImproperTypes;

my %IncorrectResidues = ( 'BAMI' => 1, 'PNTM' => 1, #seems unlikely by authors and does not exists in RESI ARG from proteins


                        );

foreach my $Chain (@$Chains) {
		
	print $Chain->getName . "\n";	
	
	next if $IncorrectResidues{$Chain->getName};
	
	my @CorrectImpropers;
		
	foreach my $Molecule (@{$Chain->getMolecules}) {
		
	  my $Impropers =  $Molecule->getImpropers;	
	
	  my $Atoms = $Molecule->getAtoms;
	
	  my %AtomNameHash = map { $_->getName => $_ } @$Atoms;
	
	  foreach my $Atom (@$Atoms) {
		
		  $SeenType{$Atom->getType}++;
		 
	  }
	
	  foreach my $Improper (@$Impropers) {
		
      my @AtomNames = @{$Improper->{'Atoms'}};
      my @ImproperAtoms;

      #print join(" ", @AtomNames) . " => ";

      my ($CenterAtomName) =shift @AtomNames;
      my $CenterAtom;

      if(exists $AtomNameHash{$CenterAtomName}) { 
	
	      $CenterAtom = $AtomNameHash{$CenterAtomName}; 
	  
	    }
	
	    elsif($CenterAtomName =~ /^\+|\-/ && exists $AtomNameHash{substr($CenterAtomName, 1) }) {
			  
			  $CenterAtom = $AtomNameHash{substr($CenterAtomName, 1)};		
						
	    }
	
	    else {
		
		    print "Cannot find Atom: " . $CenterAtomName . "\n";
		
	    }
		
      foreach my $AtomName (@AtomNames) {
	
	      if(exists $AtomNameHash{$AtomName}) {
			
	        push @ImproperAtoms, $AtomNameHash{$AtomName};      
	
	      }
	
	      elsif($AtomName =~ /^\+|\-/ && exists $AtomNameHash{substr($AtomName, 1) } ) {
		
	        push @ImproperAtoms, $AtomNameHash{substr($AtomName, 1)};
	
	      }
	
	      else { last; }
	
      }

      next if scalar(@ImproperAtoms) != 3;
     
      $SeenCenterAtoms{$CenterAtom->getType}++;			

      foreach my $Atom ($CenterAtom, @ImproperAtoms) {
	
	  #    print $Atom->getName . " ";
	
      }

     # print " => ";

			foreach my $Atom ($CenterAtom, @ImproperAtoms) {

		 #   print $Atom->getType . " ";

	    }
	
	    #print "\n";
	
	    push(@CorrectImpropers, [$CenterAtom, @ImproperAtoms]);
		
   	}
	
  	#print "Guessed Impropers\n";

    my @GuessedImpropers = FindImproperAtoms($Molecule);

    print "EXISTS : ";

    foreach my $CorrectImproper (@CorrectImpropers) {
	
	    print join(" ", map { $_->getName } @$CorrectImproper) . "   ";
			
	
    }

    print "\n";

    print "GUESSED: ";

		foreach my $GuessedImproper (@GuessedImpropers) {

		  print join(" ", map { $_->getName } @$GuessedImproper) . "   ";

	  }

	  print "\n";

    foreach my $CorrectImproper (@CorrectImpropers) {
	
	    my $Flag = 0;
	    my $StringedNameImproper = join(" ", map { $_->getName } @$CorrectImproper);
	
	    foreach my $GuessedImproper (@GuessedImpropers) {
		
		    my $StringedGuessedNameImproper =  join(" ", map { $_->getName } @$GuessedImproper);
		
		    if($StringedNameImproper eq $StringedGuessedNameImproper) {
			
			    $Flag = 1;
			    last;
			
		    } 
		    
	    }
	
	    if($Flag == 0 ) {
	
	      print "Not FOUND: " . $StringedNameImproper . "\n";
	
	    }
	
    }
		
	}
	
}

exit 1;

while (my ($key, $value) = each %SeenCenterAtoms) {
	
	print $key . " " . $value . " " . $SeenType{$key} . "\n\n";
	
	my @TypeAtoms = grep { $_->getType eq $key } @Atoms;
	
	foreach my $Atom (@TypeAtoms) {
		
		#print $Atom->getType . " " . $Atom->getMolecule->getName . "\n"; 
	  print $Atom->getName . " " . $Atom->getMolecule->getName . "\n";
	
	  
	  my $Matched = 0;
	  my $NumOfImpropers;
	
	  foreach my $ImproperType (@ImproperTypes) {
		
		  next if $ImproperType->getName ne $Atom->getType;
		
		  next unless $ImproperType->getLookUpTable->AreAllNodesSharedBetween([$Atom->getLookUpTable]);
    
	    $Matched = 1; $NumOfImpropers = $ImproperType->getBondNumber; last;
		
	  }
	
	  next if $Matched == 0;
	
	  my @MatchingImpropers;
	
	  my $ImproperCount = 0;
	
	  foreach my $Improper (@Impropers) {
		
	    if($Atom->getType eq $Improper->[0] || $Atom->getType eq $Improper->[3]) {
		
		    my $BondedAtoms = $Atom->getBondedAtoms;
		   
		    my @ToMatch;
		 
		    if($Atom->getType eq $Improper->[0]) {
			
		      @ToMatch = @$Improper;
				
		    }
		
		    else {
			
			    @ToMatch = @$Improper;
		      @ToMatch = reverse(@ToMatch);
			 
		    }
		 
		    shift @ToMatch;
	     
	      my @Matches;	
	      my @PossibleMatches;
		
	      foreach my $ImproperNum (0 .. @ToMatch-1) {
		
		      foreach my $BondedAtom (@$BondedAtoms) {
			
			      if($ToMatch[$ImproperNum] eq $BondedAtom->getType || $ToMatch[$ImproperNum] eq 'X' ) {
			
			        push @Matches, [{$ImproperNum => 1, $BondedAtom->getName => 1}, [[$ImproperNum, $BondedAtom->getName]]];
			        push @PossibleMatches, [$ImproperNum, $BondedAtom->getName ];
			
		        }
			 
		      } 
	
	      }
	
	      my $Done = 0;
	      my @NewMatches;	
	      my %SeenMatches;
		      
	      while(!$Done) {
		
		      @NewMatches = ();
		      my $FoundNew = 0;
				
		      foreach my $Match (@Matches) {
			
			      if(@{$Match->[1]} == 3 ) {
				
			        push @NewMatches, $Match; next;
			
			      }
	
	          foreach my $PossibleMatch (@PossibleMatches) {
		
		          next if exists $Match->[0]->{$PossibleMatch->[0]} || exists $Match->[0]->{$PossibleMatch->[1]};
				
							my $NewMatch;
							my %NewHash;
							my @NewList = map { $_ } @{$Match->[1]};
							
							foreach my $key (keys %{$Match->[0]}) { $NewHash{$key} = 1 }
							
							$NewHash{$PossibleMatch->[0]} = 1; $NewHash{$PossibleMatch->[1]} = 1;
							
							push @NewList, $PossibleMatch;
											
							next if exists $SeenMatches{ join(" ", sort map { join("-", @$_)} @NewList ) };
														
							push @NewMatches, [\%NewHash, \@NewList];
		
		          $SeenMatches{ join(" ", sort map { join("-", @$_) } @NewList ) } = 1;
		
		          $FoundNew = 1;
				
	          }
	
          }

          @Matches = @NewMatches;

          $Done = 1 unless $FoundNew;
	
        }


        foreach my $Match (@Matches) {
	
          if(@{$Match->[1]} == 3 && $ImproperCount <= $NumOfImpropers) {
	
	          push(@MatchingImpropers, $Improper);
			
          }

        }
		
	      
		
		
	    } 
		
	  }
	
    foreach my $MatchingImproper (@MatchingImpropers) {
	
	    print join(" ", @$MatchingImproper) . "\n" if $ImproperCount < $NumOfImpropers;
	
	    $ImproperCount++;
	
    }


  }

  print "\n";
	
}

sub FindImproperAtoms {
	
	my $Molecule = shift;
	
	my $Atoms = $Molecule->getAtoms;
	
	my @GuessedImpropers;
	
	foreach my $Atom (@$Atoms) {
		
		my $Matched = 0;
		my $ImproperNum;
	
	  foreach my $ImproperType (@ImproperTypes) {
		
      next if $ImproperType->getName ne $Atom->getType;
		
		  next unless $ImproperType->getLookUpTable->AreAllNodesSharedBetween([$Atom->getLookUpTable]);
    
	    $Matched = 1; $ImproperNum = $ImproperType->getBondNumber; last;
		
	  }
	
	  next if $Matched == 0;
	 
	  next if $ImproperNum == 0;
	
	  #print $Atom->getName . "\n";
	
	  my $BondedAtoms = $Atom->getBondedAtoms;
	
	  my %InAtomRing;
	  
	  foreach my $Atom (@$BondedAtoms) {
		
		  $InAtomRing{$Atom->getName} = 0;
		
	  }
	
	  if($Atom->getRingAtom) {
		
		  my $Rings = $Molecule->getRings;
				
		  foreach my $Ring (@$Rings) {
			
			  my $IsAtomInTheRing = 0;
			
			  foreach my $RingAtom (@$Ring) {
				
				  if($RingAtom eq $Atom) {
					
					  $IsAtomInTheRing = 1;
					  last; 
					
				  }
				
			  }
			
			  if($IsAtomInTheRing) {
				
				  foreach my $RingAtom (@$Ring) {
					
					  $InAtomRing{$RingAtom->getName} = 1;
					
				  }
				
			  }
			
		  }
			
	  }
	
	  my @SortedAtoms = sort { $InAtomRing{$b->getName} <=> $InAtomRing{$a->getName} || $b->getSumofBondTypes <=> $a->getSumofBondTypes } @$BondedAtoms;
		
		unshift (@SortedAtoms, $Atom);
		
	  my $Improper = \@SortedAtoms;
	
	  push @GuessedImpropers, $Improper;
	
	  if($ImproperNum == 2) {
		
		  my $NewImproper = [$Improper->[0], $Improper->[2], $Improper->[1], $Improper->[3]];
		
		  push @GuessedImpropers, $NewImproper;
		
	  }
	
	  elsif($ImproperNum == 3) {
		
		  my $NewImproper = [$Improper->[0], $Improper->[1], $Improper->[3], $Improper->[2]];

			push @GuessedImpropers, $NewImproper;
		
	  }
	
	  #print join(" ", ($Atom->getName, map { $_->getName} @SortedAtoms )) . "\n";
	
	}
	
	return @GuessedImpropers;
	
}

sub FindImpropers {

  my $Atom = shift;
	
  my @MatchingImpropers;
		
	foreach my $Improper (@Impropers) {
		
	 if($Atom->getType eq $Improper->[0] || $Atom->getType eq $Improper->[3]) {
		
	   my $BondedAtoms = $Atom->getBondedAtoms;
		   
	   my @ToMatch;
		 
	 	 if($Atom->getType eq $Improper->[0]) {
			
		   @ToMatch = @$Improper;
				
		 }
		
		 else {
			
		   @ToMatch = @$Improper;
		   @ToMatch = reverse(@ToMatch);
			 
		 }
		 
		 shift @ToMatch;
	     
	   my @Matches;	
	   my @PossibleMatches;
		
	   foreach my $ImproperNum (0 .. @ToMatch-1) {
		
		   foreach my $BondedAtom (@$BondedAtoms) {
			
			   if($ToMatch[$ImproperNum] eq $BondedAtom->getType || $ToMatch[$ImproperNum] eq 'X' ) {
			
			     push @Matches, [{$ImproperNum => 1, $BondedAtom->getName => 1}, [[$ImproperNum, $BondedAtom->getName]]];
		       push @PossibleMatches, [$ImproperNum, $BondedAtom->getName ];
			
		     }
			 
		   } 
	
  	 } 
	
	   my $Done = 0;
	   my @NewMatches;	
	   my %SeenMatches;
		      
	   while(!$Done) {
		
	 	   @NewMatches = ();
		   my $FoundNew = 0;
				
		   foreach my $Match (@Matches) {
			
			   if(@{$Match->[1]} == 3 ) {
				
			     push @NewMatches, $Match; next;
			
			 }
	
	     foreach my $PossibleMatch (@PossibleMatches) {
		
		     next if exists $Match->[0]->{$PossibleMatch->[0]} || exists $Match->[0]->{$PossibleMatch->[1]};
				
			   my $NewMatch;
			   my %NewHash;
				 my @NewList = map { $_ } @{$Match->[1]};
							
				 foreach my $key (keys %{$Match->[0]}) { $NewHash{$key} = 1 }
							
				   $NewHash{$PossibleMatch->[0]} = 1; $NewHash{$PossibleMatch->[1]} = 1;
							
					 push @NewList, $PossibleMatch;
											
					 next if exists $SeenMatches{ join(" ", sort map { join("-", @$_)} @NewList ) };
														
					 push @NewMatches, [\%NewHash, \@NewList];
		
		       $SeenMatches{ join(" ", sort map { join("-", @$_) } @NewList ) } = 1;
		
		       $FoundNew = 1;
				
	       }
	
       }

       @Matches = @NewMatches;

       $Done = 1 unless $FoundNew;
  	
     }   

     foreach my $Match (@Matches) {
	
	     if (@{$Match->[1]} == 3) {

	       push(@MatchingImpropers, [$Improper, [$Atom->getName, map { $_->[1] } @{ $Match->[1] } ] ] );
	
       }

		 } 
		
		}

   }
	
	return @MatchingImpropers;
  	
	
}
