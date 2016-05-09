package ChemicalSpaceMatch;

use strict;
use warnings;

use lib $ENV{'PerlChemistry'} . "/lib";

use BaseObject ':all';

require Exporter;

our @ISA = qw(BaseObject Exporter);

sub New { 
	
	my ($Class, $PAtom, $SAtom, $Level) = @_;
	
	my $Self = {
		
    _AtomDistance          => $Level,
    _IsAPrerequisite       => 0,
		_Prerequisites         => [ ],
		_APrerequisitesOf      => [ ],
		_PAtom                 => $PAtom,
		_SAtom                 => $SAtom
		
	};
	
	bless $Self, $Class;
	
	return $Self;
	
}

sub Stringify {
	
	my $Self = shift;
	
	return $Self->getPAtom->getName . " " . $Self->getSAtom->getName;
	
}

sub AddAPrerequisiteOf {
	
	my ($Self, $ChemicalSpaceMatch) = @_;
	
	push @{$Self->getAPrerequisitesOf}, $ChemicalSpaceMatch;
	
}

sub CheckForPrerequisite {
	
	my ($Self, $ChemicalSpaceMatch) = @_;
	
	if($Self->{_AtomDistance} > $ChemicalSpaceMatch->{_AtomDistance} && $Self->{_PAtom}->IsBondedTo($ChemicalSpaceMatch->{_PAtom}) && $Self->{_SAtom}->IsBondedTo($ChemicalSpaceMatch->{_SAtom}) ) {
		
    push @{$Self->getPrerequisites}, ($ChemicalSpaceMatch, @{$ChemicalSpaceMatch->getPrerequisites});

    $ChemicalSpaceMatch->setIsAPrerequisite(1);
		
	}
	
	elsif($Self->{_AtomDistance} > $ChemicalSpaceMatch->{_AtomDistance} && $Self->{_SAtom}->IsBondedTo($ChemicalSpaceMatch->{_PAtom}) && $Self->{_PAtom}->IsBondedTo($ChemicalSpaceMatch->{_SAtom}) ) {
	
    push @{$Self->getPrerequisites}, ($ChemicalSpaceMatch, @{$ChemicalSpaceMatch->getPrerequisites});

    $ChemicalSpaceMatch->setIsAPrerequisite(1);

  }  
	
}

1;