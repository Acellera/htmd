package MoleculeFileHandler::RtfFileHandler;

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

#use 5.010000;
use strict;
use warnings;
use Carp;
use MoleculeFileHandler;
use Data::Dumper;

require Exporter;

our @ISA = qw(MoleculeFileHandler Exporter);

# Items to export into callers namespace by default. Note: do not export
# names by default without a very good reason. Use EXPORT_OK instead.
# Do not simply export all your public functions/methods/constants.

# If you do not need this, moving things directly into @EXPORT or @EXPORT_OK
# will save memory.
our %EXPORT_TAGS = ( 'all'  => [ qw( )],
											 
										 'vars' => [ qw () ],
										
										 'func' => [ qw ()]);
							

our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

our @EXPORT = qw(
	
);

our $VERSION = '0.01';

#Exported Variables


#Non Exported Variables


# Preloaded methods go here.

=head2 Initiate;

Usage: $RtfFileHandlerInstance->Initiate;

Arguments:

Synopsis:
  Setups $RtfFileHandlerInstance must be called

=cut

sub Initiate {
	
	my $Self  = shift;

  my $Atom = qr/^\s*ATOM\s+(\S+)\s+(\S+)\s+(-?\d*\.\d+)/;
  my @AtomVariables = qw( Name Type Charge);

  my $Bond = qr/^\s*BOND\s+([^!]+)/;
	my @BondVariables = qw(BondAtomNames Type);
	
	my $Double = qr/^\s*(?:DOUB|DOUBLE)\s+([^!]+)/;
	
	my $Triple = qr/^\s*(?:TRIP|TRIPLE)\s+([^!]+)/;
	
	
	my $EndChain = qr/^(?:IC|BILD)/;
	
	my $Resi = qr/^\s*(RESI|PRES)\s*(\S+)\s+(-?\d+\.\d+)/;
	my @ResiVariables = qw(IsPatch ResidueName ResidueCharge);
	
	my $Improper = qr/^\s*IMPR\s+([^!]+)/;
	
	my $Delete = qr/^\s*(?:DELE|DELETE)\s+([^!]+)/;
	my @DeleteVariables = qw(Type Name);
  
  my $AtomLineParser =    LineParser->New('Atom',
																		      sub { uc($_[0]) =~ $Atom },
																		      sub { return { map {$AtomVariables[$_] => $_[0]->[$_] } (0 .. $#AtomVariables) } } );														

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

	my $DoubleLineParser =  LineParser->New('Bond',
																			    sub { 
																			      return () if uc($_[0]) !~ $Double;
																				    my @spl = split /\s+/, $1;
																		        return map { [ [$spl[$_], $spl[$_ + 1]], 2] } grep { $_ % 2 == 0 } (0 .. $#spl);																																   
																			    },
																			    sub { 
																			      my @Hashes;
																				    foreach my $Bond (@{$_[0]}) {
																				      push @Hashes, { map {$BondVariables[$_] => $Bond->[$_] } (0 .. $#BondVariables) }
																				    }
																				    return \@Hashes;
																			    } );

	my $TripleLineParser =  LineParser->New('Bond',
																					sub { 
																					  return () if uc($_[0]) !~ $Triple;
																						my @spl = split /\s+/, $1;
																						return map { [ [$spl[$_], $spl[$_ + 1]], 3] } grep { $_ % 2 == 0 } (0 .. $#spl);																																   
																					},
																					sub { 
																					  my @Hashes;
																						foreach my $Bond (@{$_[0]}) {
																						  push @Hashes, { map {$BondVariables[$_] => $Bond->[$_] } (0 .. $#BondVariables) }
																						}
																						return \@Hashes;
																					} );


 	my $ResiLineParser =    LineParser->New('Resi',
																			    sub { uc($_[0]) =~ $Resi },
																			    sub { return { map {$ResiVariables[$_] => $_[0]->[$_] } (0 .. $#ResiVariables) } } );

  my $EndChainLineParser= LineParser->New('EndChain',
																			    sub { uc($_[0]) =~ $EndChain },
																			    sub { } );
																			
	my $BlankLineParser   = LineParser->New('BlankLine',
																					sub { $_[0] =~ qr/^\s*\n/ },
																					sub { } );
	
	my $PatchLineParser   = LineParser->New('EndChain',
																					sub { $_[0] =~ qr/^\s*PATC/ },
																					sub { } );
																					

	my $ImprLineParser = LineParser->New('Improper',
																	 sub { 
																	   return () if uc($_[0]) !~ $Improper;
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
																			
	my $DeleteLineParser = LineParser->New('Delete', 
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
	
	$Self->{_Patterns} = [ $AtomLineParser,
												 $BondLineParser,
												 $DoubleLineParser,
												 $TripleLineParser,
												 $ImprLineParser,
												 $PatchLineParser,
												 $DeleteLineParser,
												 $ResiLineParser,
												 $EndChainLineParser,
												 $BlankLineParser ];
												

}

=head2 ReadFile;

Usage: $PdbFileHandlerInstance->ReadFile($File);

Arguments:

Synopsis:
  This is where all the data is collected from File

=cut

sub ReadFile {
			
	my ($Self, $FilePath) = @_;
	
	open(FILE, $FilePath);
	
	my @Chains;
	my $Finished = 0;
	
	while(!$Finished) {
		
	  push @Chains, $Self->ReadChainFromFile(*FILE); 
	
	  $Finished = 1 if($Chains[-1]->[0] eq 'EOF');
	
  }
	
	close(FILE);
	
	my @FinishedChains;

	#Specific Post Processing for RTFs, they are strange and have a lot of extra useful info	
	foreach my $Chain (@Chains) {
		
		next unless defined $Chain->[1]->{'Resi'}[0]->{'ResidueName'};
				
	  my $ResidueName = $Chain->[1]->{'Resi'}[0]->{'ResidueName'};
	
	  my $ResidueIndentifier = $ResidueName . " " . 1;
		
		$Chain->[1]->{'Molecules'}{$ResidueIndentifier}{'Structure'} = 'Molecule';
		$Chain->[1]->{'Molecules'}{$ResidueIndentifier}{'IsPatch'} =  $Chain->[1]->{'Resi'}[0]->{'IsPatch'} eq 'RESI' ? 0 : 1;
	  $Chain->[1]->{'Molecules'}{$ResidueIndentifier}{'Bonds'} = delete $Chain->[1]->{'Bond'};
	  $Chain->[1]->{'Molecules'}{$ResidueIndentifier}{'Atoms'} = delete $Chain->[1]->{'Atom'};
	  $Chain->[1]->{'Molecules'}{$ResidueIndentifier}{'Charge'} = $Chain->[1]->{'Resi'}[0]->{'ResidueCharge'};
	  $Chain->[1]->{'Molecules'}{$ResidueIndentifier}{'Name'} =  $ResidueName; 
	  $Chain->[1]->{'Molecules'}{$ResidueIndentifier}{'Num'} = 1;
	  $Chain->[1]->{'Molecules'}{$ResidueIndentifier}{'Deletes'} = delete $Chain->[1]->{'Delete'};
	  $Chain->[1]->{'Molecules'}{$ResidueIndentifier}{'Impropers'} = delete $Chain->[1]->{'Improper'};
	
	
	  push @FinishedChains, $Chain->[1];
    	
	}
	
	#print "Finished: " . scalar(@Chains) . "\n";
	
	$Self->setChains(\@FinishedChains);
  return 1;

}


=head2 ReadChainFromFile;

Usage: $PdbFileHandlerInstance->ReadChainFromFile($File);

Arguments:

Synopsis:
  Collects the data on a specific Chain, this allows the break up of the data

=cut

sub ReadChainFromFile {
	  
  my ($Self, $FH) = @_;
		
	my $Patterns = $Self->getPatterns;

	my $Sections;
				
	my $CurrentRepresentation = "";
	my $LastRepresentation = "";
		
	my $DoneReading = 0;
		
	while(my $Line = <$FH>) {
			
	  if(eof($FH)) { $CurrentRepresentation = 'EOF'; last; }
													
		foreach my $LineParser (@$Patterns) {
				
		  next if ! $LineParser->Parse($Line);
		
		  if(exists $Sections->{'Resi'} and $LineParser->getRepresentation eq 'Resi') {
			
			  seek($FH, -length($Line), 1);
			
			  $DoneReading = 1; last;
			
		  }		  
					
		  my $DataHash = $LineParser->ReturnResults;
								
		  if($DataHash) {				
			 			
		    push @{$Sections->{$LineParser->getRepresentation}}, $DataHash  if ref($DataHash) eq 'HASH';
			  push @{$Sections->{$LineParser->getRepresentation}}, @$DataHash if ref($DataHash) eq 'ARRAY';
		  }
			
			$LastRepresentation = $CurrentRepresentation;								
		  $CurrentRepresentation = $LineParser->getRepresentation;
			  			
		  last;
			  					
	  }		
	
	  last if $DoneReading;
										
		
	}		
			
	return [$CurrentRepresentation, $Sections];
}


1;
__END__
