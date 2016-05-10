package MATCHBaseObject;


#use 5.010000;
use strict;
use warnings;
use Carp;

confess("PerlChemistry Library path is not set") unless exists $ENV{'PerlChemistry'};

use lib $ENV{'PerlChemistry'} . "/lib";

use BaseObject;

require Exporter;

our @ISA = qw(BaseObject Exporter);

# Items to export into callers namespace by default. Note: do not export
# names by default without a very good reason. Use EXPORT_OK instead.
# Do not simply export all your public functions/methods/constants.

# This allows declaration	use MATCH ':all';
# If you do not need this, moving things directly into @EXPORT or @EXPORT_OK
# will save memory.
our %EXPORT_TAGS = ( 'all'  => [ qw( $Parameters
	                                   $Verbosity
															  		 PrintMessage)],
											 
										 'vars' => [ qw ($Parameters
																		 $Verbosity) ],
										
										 'func' => [ qw (PrintMessage)]);
							

our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';


# Preloaded methods go here.

1;
__END__
# Below is stub documentation for your module. You'd better edit it!

=head1 NAME

MATCH - Perl extension for blah blah blah

=head1 SYNOPSIS

  use MATCH;
  blah blah blah

=head1 DESCRIPTION

Stub documentation for MATCH, created by h2xs. It looks like the
author of the extension was negligent enough to leave the stub
unedited.

Blah blah blah.

=head2 EXPORT

None by default.



=head1 SEE ALSO

Mention other useful documentation such as the documentation of
related modules or operating system documentation (such as man pages
in UNIX), or any relevant external documentation such as RFCs or
standards.

If you have a mailing list set up for your module, mention it here.

If you have a web site set up for your module, mention it here.

=head1 AUTHOR

Joseph Yesselman, E<lt>skullnite@apple.comE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2010 by Joseph Yesselman

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.10.0 or,
at your option, any later version of Perl 5 you may have available.


=cut
