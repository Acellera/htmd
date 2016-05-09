#!/usr/bin/env perl

=head1 NAME


=head1 SYNOPSIS


=head1 DESCRIPTION


=head2 EXPORT


=head1 AUTHOR

Joseph Yesselman, E<lt>jyesselm@umich.eduE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2010 by Joseph Yesselman
c
This library is free software; you can redistribute it and/or modifycd
it under the same terms as Perl itself, either Perl version 5.10.0 or,
at your option, any later version of Perl 5 you may have available.

=head1 FUNCTIONS

=cut

#use 5.010000;
use strict;
use warnings;
use Carp;

#Load MATCH Libraries 
use lib $ENV{"MATCH"} . "/lib";

use MATCHBaseObject;
use MATCHFunctions ':func';
use LookUpTable;
use Storable;
use Parameters;
use MATCHParameters;
use BaseObject ':vars';

our $VERSION = '0.01';
#Exported Variables


#Non Exported Variables

my $DefaultParameters = MATCHParameters->New;

$DefaultParameters->Initiate;

$Parameters = $DefaultParameters;


# Preloaded methods go here.

srand(time ^ ($$ + ($$ << 15)));

my $FileName = $ARGV[0];

my $TopologyFilePath = $ENV{"MATCH"} . "/resources/$FileName/$FileName.rtf";
my $TopologyPatchingFile = $ENV{"MATCH"} . "/resources/$FileName/$FileName.patches";

my $Chains = SetupMoleculesFromTopology($TopologyFilePath, $TopologyPatchingFile);

open(FILE, $ENV{"MATCH"} . "/resources/Incorrect_Residue_List.txt");
my %IncorrectResidues = map { chomp($_); $_ => 1 } <FILE>;
close(FILE);

print scalar(@$Chains) . "\n";
my @CorrectChains = grep { ! exists $IncorrectResidues{$_->getMolecule(0)->getName } } @$Chains;
print scalar(@CorrectChains) . "\n";
#my $Chains = retrieve("Chains.dat");

store(\@CorrectChains, $ENV{"MATCH"} . "/t/StoredTopologys/$FileName.dat");

exit 1;

1;
__END__
