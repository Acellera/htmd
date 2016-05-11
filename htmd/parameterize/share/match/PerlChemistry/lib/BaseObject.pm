package BaseObject;

=head1 NAME

BaseObject - The Abstract Class that all object in PerlChemistry are derived from 

=head1 SYNOPSIS

use BaseObject;
@ISA = qw(BaseObject);

=head1 DESCRIPTION

BaseObject is to be inherited by all other objects in the Module

=head2 EXPORT

:all => $Parameters, $Verbosity, PrintMessage($$)
:vars => $Parameters, $Verbosity
:func => PrinMessage($$)

=head1 AUTHOR

Joseph Yesselman, E<lt>jyesselm@umich.eduE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2010 by Joseph Yesselman

This library is free software; you can redistribute it and/or modified
it under the same terms as Perl itSelf, either Perl version 5.10.0 or,
at your option, any later version of Perl 5 you may have available.

=head1 FUNCTIONS

=cut

#use 5.010000;
use strict;
use warnings;
use Carp;

require Exporter;

our @ISA = qw(Exporter);

# Items to export into callers namespace by default. Note: do not export
# names by default without a very good reason. Use EXPORT_OK instead.
# Do not simply export all your public functions/methods/constants.

# This allows declaration	use BaseObject ':all';
# If you do not need this, moving things directly into @EXPORT or @EXPORT_OK
# will save memory.
our %EXPORT_TAGS = ( 'all'  => [ qw( $Parameters
	                                   $Verbosity
															  		 PrintMessage)],
											 
										 'vars' => [ qw ($Parameters
																		 $Verbosity) ],
										
										 'func' => [ qw (PrintMessage)]);
							

our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

our @EXPORT = qw(
);

our $VERSION = '0.01';

#Exported Variables
our $Parameters;  #Contains all program parameters in Parameters object can be altered by reading a .par file
our $Verbosity  = 5;

#Non Exported Variables
our $AUTOLOAD;

# Preloaded methods go here.

=head2 PrintMessage($$)

Arguments:
  $Message: The String on which to print
  $Level  : The Verbosity Level, if $Level is less than $Verbosity then it will not be printed

Synopsis:
  Print outs a message to the screen at a given verbosity level

=cut 

sub PrintMessage($$) {

  my ($Message, $Level) = @_;

  confess("Level must be a number!\n") if $Level !~ /^\d/;

  print $Message . "\n" if $Level <= $Verbosity;

}

sub AUTOLOAD {
  my ($Self, $newvalue) = @_;
  my $type = ref($Self) 
             or print "$Self is not an object\n";


  if($AUTOLOAD =~ /DESTROY/ ) { return; } 

  my ($operation, $attribute) = ($AUTOLOAD =~ /(get|set)(\w+)$/);
  
  confess("Method name $AUTOLOAD in " . ref($Self) . " is not in the recognized form (get|set)_attribute") unless ($operation && $attribute);
	
	confess("Parameter $attribute does not exist in " . ref($Self)) unless (exists $Self->{"_" . $attribute});
	
  no strict "refs";
 
  if($operation eq 'set')  { 
	  $Self->{"_" . $attribute} = $newvalue;

	  *{$AUTOLOAD} = sub { 
		  						   my ($Self, $newvalue) = @_;
		  							 $Self->{"_" . $attribute} = $newvalue; 	
								   };
  }

  elsif($operation eq 'get')    { 
	  *{$AUTOLOAD} = sub { 
		                 my ($Self) = @_;
									   $Self->{"_" . $attribute};
									 }; 
	}

  use strict "refs";

  #confess("Cannot retrieve parameter: " . $attribute . " because it is undefined, add its value to your .par file\n") if(! defined $Self->{"_" . $attribute});

  return $Self->{"_" . $attribute};
}





1;
__END__


