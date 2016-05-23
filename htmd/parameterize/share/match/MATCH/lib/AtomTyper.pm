package AtomTyper;

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
use MATCHBaseObject;
use MATCHFunctions ':all';
use BaseObject ':vars';
use Type;

require Exporter;

our @ISA = qw(MATCHBaseObject Exporter);

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

our $DistanceBondingRules = { };


# Preloaded methods go here.

=head2 New

Usage: AtomTyper->New;

Arguments:
  $Class: should be 'AtomTyper'

Synopsis:
  creates a new AtomTyper object

=cut

sub New {
	
	my $Class = shift;
	
	my $Self = {
		
		_Types 		=>   [],
		
	};
	
	
	bless $Self, $Class;
	return $Self;
	
}

=head2 Initiate

Usage: $AtomTyperInstance->Initiate;

Arguments:

Synopsis:
  Sets up new AtomTyperInstance must be called!

=cut

sub Initiate {
	
	my $Self = shift;
	
	my $TypesFilePath = $Parameters->getTypesFilePath;
	
	open(FILE, $TypesFilePath);
	
	my @TypesFileContents = grep { $_ !~ /^(?:\#| )/} <FILE>;
	
	close(FILE);
	
	my @Types;
	
	my @Subs;
	my @Lines;
	
	foreach my $Line (@TypesFileContents) {
		
		if($Line =~ /\s*\S+\s+=\s+\S+\s+/) {
			
			push @Subs, $Line;
			
		}
		
		else { push @Lines, $Line }
		
	}
	
	my $Lines = ApplyStringNamingShortCuts(\@Subs, \@Lines);
	
	foreach my $Line (@$Lines) {
				
    my @ArrayedLine = split (/\s+/, $Line);

    next if @ArrayedLine != 3;

    my $Type = Type->New(@ArrayedLine);

    $Type->Initiate;

    push @Types, $Type;
		
	}
	
	my @SortedTypes = sort { length($b->getString) <=> length($a->getString) } @Types;
	
	$Self->setTypes(\@SortedTypes);
	
	
}


sub TypeAtomsInChain {
	
	my ($Self, $Chain) = @_;
	
	my @Atoms = map { @{$_->getAtoms} } @{$Chain->getMolecules};
	
	return $Self->TypeAtoms(\@Atoms);
	
}

=head2 TypeAtomsInMolecule

Usage: $AtomTyperInstance->TypeAtomsInMolecule;

Arguments:
  $Molecule: The molecule in which you to get the Atom types for

Synopsis:
  Assigns the Atom type for all atoms in Molecule object

=cut

sub TypeAtomsInMolecule {
	
	my ($Self, $Molecule) = @_;
	
	my $Atoms = $Molecule->getAtoms;
	
  return $Self->TypeAtoms($Atoms);
	
}


sub TypeAtoms {
	
	my ($Self, $Atoms) = @_;
	
	my $Types = $Self->getTypes;
	
	my $ExitifNotTyped = $Parameters->getExitifNotTyped;
	
	my @CannotType;
	
	foreach my $Atom (@$Atoms) {
		
		my $TypeMatchToAtom = 0;
				
		foreach my $Type (@$Types) {
		
      next unless $Type->getLookUpTable->AreAllNodesSharedBetween([$Atom->getLookUpTable]);

      $Atom->setType($Type->getName);

      $TypeMatchToAtom = 1; 

      if($Verbosity > 3 ) {

        print sprintf("%-6s",$Atom->getName) . " ";
	      print sprintf("%-6s",$Atom->Stringify) . " " . sprintf("%-6s", $Type->getName) . " " . sprintf("%-6s", $Type->getString);
		  	print "\n";
			
		  }

      last;

    }

    

    croak("Could Not Type Atom " . $Atom->getName . " With State " . $Atom->Stringify ) if ! $TypeMatchToAtom && $ExitifNotTyped;
		
		push @CannotType, $Atom if ! $TypeMatchToAtom && !$ExitifNotTyped;
				
	}
	
	if(!$ExitifNotTyped && scalar(@CannotType) > 0) {
		
		return -1;
		
	} 
	
	$Self->PerformForceFieldSpecificNonFragmentTyping($Atoms);
	
	return 1;

}


=head2 PerformForceFieldSpecificNonFragmentTyping

Usage: $AtomTyperInstance->PerformForceFieldSpecificNonFragmentTyping($Atoms);

Arguments:
  $Atoms: an array of the atoms within a chain/molecule 

Synopsis:
  Handles all typing procedures that are not done via straight up fragment matching

=cut

sub PerformForceFieldSpecificNonFragmentTyping {
	
	my ($Self, $Atoms) = @_;
	
	if(uc($Parameters->getForceField) eq "CGENFF") {
		
		#Alternates CG2DC1 / CG2DC2 for conjugated alkenes
		HandleCG2DCTypes($Atoms);
		
	}
 	
}


sub HandleCG2DCTypes {
	
	my $Atoms = shift;
	
  my %SeenCG2DC;

	my %AcceptableEnd = map { $_ => 1 } qw(CG2D1O CG2D2O CG2DC3);
	
	foreach my $Atom (@$Atoms) {
		
		#Only need to consider atoms with type CG2DC? and that have not been encountered yet 
		next if $Atom->getType ne "CG2DC?" || exists $SeenCG2DC{$Atom};

		$SeenCG2DC{$Atom} = 1;

    #Start the chain
		my @Chain = ($Atom);

		my $Finished = 0;

		while(!$Finished) {

		  $Finished = 1;

		  my $StartBonds = $Chain[0]->getBondedAtoms;

      #Look to see if any of the bonded atoms at the start and the end of the current chain are also type "CG2DC?" to extend it
		  foreach my $Bonded (@$StartBonds) {

			  next if $Bonded->getType ne "CG2DC?" || exists $SeenCG2DC{$Bonded};

			  unshift @Chain, $Bonded; $Finished = 0;

			  $SeenCG2DC{$Bonded} = 1;

			  last;

		  }

		  my $EndBonds = $Chain[-1]->getBondedAtoms;

		  foreach my $Bonded (@$EndBonds) {

			  next if $Bonded->getType ne "CG2DC?" || exists $SeenCG2DC{$Bonded};

			  push @Chain, $Bonded; $Finished = 0;

			  $SeenCG2DC{$Bonded} = 1;

			  last;

		  }  	

		}
		
    #Types that are bonded to the chain of CG2DC1/CG2DC2 are strong indicators of whether to use 1 or 2
    my $Start = ""; my $End = "";

		foreach my $Bonded (@{$Chain[0]->getBondedAtoms}) {

	  	$Start = $Bonded->getType if exists $AcceptableEnd{$Bonded->getType};

		}

		foreach my $Bonded (@{$Chain[-1]->getBondedAtoms}) {

			$End = $Bonded->getType if(exists $AcceptableEnd{$Bonded->getType});

		}
		
		#Currently only have a ragtag series of conditions, probably should work on getting some real rules!
		if (@Chain == 1) {

			$Chain[0]->setType("CG2DC1");

		}

    #CG2DC3 is always bonded to CG2DC2, except in RESI MECH
		elsif($Start eq "CG2DC3" && $End eq "CG2DC3") {

      #Organize by name, so C1 beats C2
			if($Chain[0]->getName ge $Chain[-1]->getName) {

				@Chain = reverse(@Chain);

			}

			$Chain[0]->setType("CG2DC2");

		}

		elsif($End eq "CG2DC3") {

		  @Chain = reverse(@Chain);

			$Chain[0]->setType("CG2DC2");

		}

		elsif($Start eq "CG2D1O") {

			$Chain[0]->setType("CG2DC1");

		}

		elsif($End eq "CG2D1O") {

			@Chain = reverse(@Chain);

			$Chain[0]->setType("CG2DC1");

		}

    $Chain[0]->setType("CG2DC2") if $Chain[0]->getType eq "CG2DC?";
		
		foreach my $i (1 .. @Chain-1) {

	    my $Same; my $Alt;

			if($Chain[$i-1]->getType eq "CG2DC1") {

				$Same = "CG2DC1"; $Alt = "CG2DC2";

			} 

			else {

				$Same = "CG2DC2"; $Alt = "CG2DC1";

			}

	    my $CurrentBondType = $Chain[$i-1]->getBondwithAtom($Chain[$i])->getType;

			if($CurrentBondType == 2) {

				$Chain[$i]->setType($Same);

			}

			else {

				$Chain[$i]->setType($Alt);

			}

		}

	}
	
}

1;
__END__

