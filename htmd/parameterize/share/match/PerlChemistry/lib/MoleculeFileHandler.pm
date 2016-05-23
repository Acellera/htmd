package MoleculeFileHandler;
use Data::Dumper;

=head1 NAME

MoleculeFileHandler - Abstract class that lays out the frame work for all MoleculeFileHandler classes 

=head1 SYNOPSIS

use MoleculeFileHandler;
@ISA = qw(MoleculeFileHandler);

=head1 DESCRIPTION

MoleculeFileHandler is the abstract class that is inherited by all MoleculeFileHandler classes

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
use LineParser;
use Complex;

require Exporter;

our @ISA = qw(BaseObject Exporter);

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


=head2 New

Usage: MoleculeFileHandler->New;

Arguments:
  $Class: should be 'MoleculeFileHandler'
  $FilePath: the location of the File trying to be loaded

Synopsis:
  Creates MoleculeFileHandler object, this object is just an abstract class, New actually allows for polymorphism and creates a class that handles that
specific MoleculeFile type 

=cut

sub New {
	
	my $Class = shift;
	my $FilePath = shift;
	
	croak("$FilePath cannot be loaded as it does not exist! ") unless -e $FilePath;
	
	my $Self = {
		
	  _Chains								=> [],
	  _Complexes						=> [],
	  _FileName							=> undef,
	  _Patterns							=> [],
		
	};
	
	my $FileHandlerType = DetermineWhatFileHandlerIsNeccessary($FilePath);

	eval "require $FileHandlerType; ";
	
  croak("$FileHandlerType does not exist!") if $@;
	
	bless $Self, $FileHandlerType;
	
	my @FileSplitOverSlash = split("/", $FilePath); 
	
	my $LastPart = pop @FileSplitOverSlash;
			
	my ($FileName) = ($LastPart =~/(\S+)\.\S+$/ );
	
  $Self->setFileName($FileName);
	
	$Self->Initiate;
	$Self->ReadFile($FilePath);
	
	return $Self;
	
	
}

=head2 DetermineWhatFileHandlerIsNeccessary($)

Usage: MoleculeFileHandler->DetermineWhatFileHandlerIsNeccessary($FilePath)

Arguments:
  $FilePath: The file path of the file trying to be loaded in

Synopsis:
  Determines what class needs to be used to load in this file 

=cut

sub DetermineWhatFileHandlerIsNeccessary($) {
	
	my $FilePath = shift;
	
	my @FileSplitOverPeriod = split(/\./, $FilePath);
	
	my $FileExtension = pop @FileSplitOverPeriod;
		
	my $Class = "MoleculeFileHandler";
	
	if(uc($FileExtension) eq "MOL2")   { $Class .= "::Mol2FileHandler"}
	elsif(uc($FileExtension) eq "SDF") { $Class .= "::SdfFileHandler"}
	elsif(uc($FileExtension) eq "MOL") { $Class .= "::SdfFileHandler"}
	elsif(uc($FileExtension) eq "PDB") { $Class .= "::PdbFileHandler"}
	elsif(uc($FileExtension) eq "PRM") { $Class .= "::PrmFileHandler"}
	elsif(uc($FileExtension) eq "RTF") { $Class .= "::RtfFileHandler"}
	elsif(uc($FileExtension) eq "STR") { $Class .= "::RtfFileHandler"}
	elsif(uc($FileExtension) eq "INP") { $Class .= "::RtfFileHandler"}
	else															 { croak "FileType " . uc($FileExtension) . " is not supported!"}
	
  return $Class;
 
}

=head2 BuildFromFile;

Usage: $MoleculeFileHandlerInstance->BuildFromFile;

Arguments:

Synopsis:
  Creates the neccessary Objects to represent the Data in the file, i.e. if its only one molecule, this function will return a new molecule object 
with all the data stored in it, see MoleculeFileHandler.t if interested

=cut

sub BuildObjectsFromFile {
	
	my $Self = shift;
	
	my $DataStructuresFromFile = $Self->ProcessDataFromFile;
	
	my %LoadedOverview; 
		
	my @StructuralObjects;
	
	PrintMessage("Warning Nothing is Loaded!!, from File: " . $Self->getFileName, 1) if !@$DataStructuresFromFile;
	
	my $Count =  0;			
				
	foreach (@$DataStructuresFromFile) {
		
		my $StructualObject;
		
		$LoadedOverview{$_->{'Structure'}}++;
		
		$_->{'FileName'} = $Self->getFileName;
		
		if($_->{'Structure'} eq 'Complex')  { $StructualObject = Complex->New($_); }
		elsif($_->{'Structure'} eq 'Chain') { $StructualObject = Chain->New($_);   }
									
		push @StructuralObjects, $StructualObject;
		
	}
		
	PrintMessage("Loaded from File: ", 6);
	
	while( my ($StructureType, $Count) = each(%LoadedOverview)) {
		PrintMessage("$StructureType => $Count", 6);
	}
		
	return \@StructuralObjects;
}

=head2 ProcessDataFromFile;

Usage: $MoleculeFileHandlerInstance->ProcessDataFromFile;

Arguments:

Synopsis:
  Gathers all data from File but does not create no objects to hold it, this is used if you want to add some new data to objects that are already created

=cut

sub ProcessDataFromFile {
	
	my $Self = shift;
	
	my $Chains = $Self->getChains;
	
	my @FinishedChains;  
	
  foreach (@$Chains) {	  
	
	  my @ChainHashes = NewProcessChainFromHashedFileData($_);
	
	  push @FinishedChains, @ChainHashes;

  }

  my $Complexes = $Self->getComplexes;

  my @FinishedComplexes;

  foreach my $Complex (@$Complexes) {
	
	  my $ComplexHash;
	  $ComplexHash->{'Structure'} = 'Complex';
		
	  foreach my $ChainNum (0 .. @$Complex - 1) {
		  
		  my @ChainHashes = NewProcessChainFromHashedFileData($Complex->[$ChainNum]);
		  my $Count = 0;
		
		  foreach (@ChainHashes) {
					
		  	if($Count == 0)  { $ComplexHash->{'Chains'}{$_} = $_; }
		    else						 { $ComplexHash->{'Chains'}{$_ . $Count} = $_; }
		
		    $Count++;
				
			}				
	  }
	
	
	  push @FinishedComplexes, $ComplexHash;
	
	  last; #!!!!!!only care about first model add flag!!!
	}
	
	if(@FinishedComplexes) {
		
	  if(! exists $FinishedComplexes[-1]->{'Chains'}) {
	  
	    splice(@FinishedComplexes, -1, 1);
	
     #my @Bonds = @{$FinishedComplexes[-1]->{'Chains'}{0}{'Bonds'}};
		
	  }
	
  }
	
	return [@FinishedChains, @FinishedComplexes];	
	
	
}


sub ExtractCoordinatesHash {
	
	my $Self = shift;
	
	my $Chains = $Self->getChains;
  
  my $Atoms;
	
  foreach my $Chain (@$Chains) {

  	if(exists $Chain->{'Molecules'}) {
	
			#$ChainHashes->{""}{'Molecules'} = delete $Chain->{'Molecules'};
			#$ChainHashes->{""}{'Structure'} = 'Chain';
		
		}
		
		while( my ($Representation, $Containers) = each(%$Chain)) {

	  	if ($Representation eq 'Atom') {
		
		  	foreach (@$Containers) {
		    
      		my $ChainName = exists $_->{'ChainName'} ? delete $_->{'ChainName'} : "";
	    		my $ResidueName = exists $_->{'ResidueName'} ? delete $_->{'ResidueName'} : "UNK";
	    		my $SegmentId = exists $_->{'SegmentId'} ? delete $_->{'SegmentId'} : "";

      	  $Atoms->{"$ChainName-$ResidueName-$SegmentId-" . $_->{'Name'}} = [$_->{'X'}, $_->{'Y'}, $_->{'Z'}];


			  }

    	}

  	}

  }
	
	return $Atoms;
}



=head2 ProcessChainFromHashedFileData($);

Usage: ProcessChainFromHashedFileData($Chain);

Arguments:
  $Chain: The results from ReadFile that is thought to be organized as a chain

Synopsis:
  Organizes Hashed data from ReadFile into A Hash that can be loaded directly into a Chain objects

=cut

sub NewProcessChainFromHashedFileData($) {

  my $Chain = shift;

  my $TempHash;

  my $ChainHashes;

	#A hack to move preprocessed Molecule Hashes directly into ChainHashes <= This is ugly
	if(exists $Chain->{'Molecules'}) {
				
		$ChainHashes->{""}{'Molecules'} = delete $Chain->{'Molecules'};
		$ChainHashes->{""}{'Structure'} = 'Chain';
	
	}
		
  while( my ($Representation, $Containers) = each(%$Chain)) {
	
	  if ($Representation eq 'Atom') {
	
	    foreach (@$Containers) {
		
	      my $ChainName = exists $_->{'ChainName'} ? delete $_->{'ChainName'} : "";
	      my $ResidueName = exists $_->{'ResidueName'} ? delete $_->{'ResidueName'} : "UNK";
	      my $SegmentId = exists $_->{'SegmentId'} ? delete $_->{'SegmentId'} : "";
	      
				#Somtimes atomic elements are stored here for some reason
				#+0 is included by some zeros I would guess for formal charge, this sometimes makes it length 4
        if(length($SegmentId) != 4 || $SegmentId =~ /[+,-]/) {
	
				  $SegmentId = "";
	
				}
	
	      confess("Atom does not have ResidueNum value!") if ! exists $_->{'ResidueNum'};
	       
	      my $ResidueNum = delete $_->{'ResidueNum'};
	
	      my $ChainIndentifier = $ChainName . " " . $SegmentId;
	
	      my $ResidueIndentifier = $ResidueName . " " . $ResidueNum;
	
	      if(! exists $ChainHashes->{$ChainIndentifier}) {
		
	        $ChainHashes->{$ChainIndentifier}{'Structure'} = 'Chain';
	        $ChainHashes->{$ChainIndentifier}{'Name'} = $ChainName;  
	        $ChainHashes->{$ChainIndentifier}{'SegmentId'} = $SegmentId;		      
	
	      }
	
	      if(! exists $ChainHashes->{$ChainIndentifier}{'Molecules'}{$ResidueIndentifier}) {
				
		      $ChainHashes->{$ChainIndentifier}{'Molecules'}{$ResidueIndentifier}{'Name'} = $ResidueName;
		      $ChainHashes->{$ChainIndentifier}{'Molecules'}{$ResidueIndentifier}{'Num'} = $ResidueNum;
		      $ChainHashes->{$ChainIndentifier}{'Molecules'}{$ResidueIndentifier}{'Structure'} = 'Molecule';
			      	
	      }
	
	      push @{$ChainHashes->{$ChainIndentifier}{'Molecules'}{$ResidueIndentifier}{'Atoms'}}, $_;        
	
	    }
	
    }
	
	}
	
	return values (%$ChainHashes);
	

}

sub ProcessChainFromHashedFileData($) {
	
	my $Chain = shift;
  my $ChainHash;
	  
	$ChainHash->{'Structure'} = 'Chain';
			  
	while( my ($Representation, $Containers) = each(%$Chain)) {
				
		if($Representation eq 'Atom') {
		    
      foreach (@$Containers) {				
               
				#Set Name for Chain if ChainName is defined for Atom 
        if(exists $_->{'ChainName'} && ! exists $ChainHash->{'Name'}) { $ChainHash->{'Name'} = delete $_->{'ChainName'};  }

				#Chain Name randomly Changes without using TER to break it throw error
        elsif(exists $_->{'ChainName'} && exists $ChainHash->{'Name'}) {
	        delete $_->{'ChainName'};
	        #confess("Inconsistent Chain Naming: "  . $_->{'ChainName'} . " " . $_->{'Num'} . "\n") if $ChainHash->{'Name'} ne delete $_->{'ChainName'};
        }

				#Check to see if SegmentId stay consistent throughout Chain
        if(exists $_->{'SegmentId'} && ! exists $ChainHash->{'SegmentId'}) { $ChainHash->{'SegmentId'} = delete $_->{'SegmentId'};  }

	      elsif(exists $_->{'SegmentId'} && exists $ChainHash->{'SegmentId'}) {
		      confess("Inconsistent Segment Ids") if $ChainHash->{'SegmentId'} ne delete $_->{'SegmentId'};
	      }

        my $ResidueNum;

        if(exists $_->{'ResidueNum'} && ! exists $ChainHash->{'Molecules'}{$_->{'ResidueNum'}} ) {
	        $ResidueNum = delete $_->{'ResidueNum'};
	        $ChainHash->{'Molecules'}{$ResidueNum}{'Name'} = exists $_->{'ResidueName'} ? delete $_->{'ResidueName'} : "UNK" . $ResidueNum;
	        $ChainHash->{'Molecules'}{$ResidueNum}{'Structure'} = 'Molecule';
		    }
				
		    elsif(exists $_->{'ResidueNum'}) { 
			    $ResidueNum = delete $_->{'ResidueNum'};  
			  
			
			  }
			
	
		
		    else { print $_->{'Num'} . "\n"; confess("Atom does not have ResidueNum value!") }
								
		    push @{$ChainHash->{'Molecules'}{$ResidueNum}{'Atoms'}}, $_;
	      
		
	    }	
		} 		  

	}
	
	return $ChainHash;
}

=head2 Initiate;

Usage: $MoleculeFileHandlerInstance->Initiate;

Arguments:

Synopsis:
  Setups $MoleculeFileHandlerInstance must be called

=cut

sub Initiate {
	
}

=head2 ReadFile;

Usage: $MoleculeFileHandlerInstance->ReadFile($File);

Arguments:

Synopsis:
  This is where all the data is collected from File

=cut

sub ReadFile {
	
}




1;
__END__
