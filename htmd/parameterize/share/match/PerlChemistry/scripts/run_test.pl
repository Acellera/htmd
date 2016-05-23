#!/usr/bin/perl

use strict;
use warnings;

my @Files = map { "~/Work/CDOCKER/COM_$_\_0.pdb" } (0 .. 100);

my $String = join(" ", @Files);


system("./Test3.pl ~/Work/CDOCKER/COM_0_0.pdb $String");
	
