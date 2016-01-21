#-*- buffer-read-only: t -*-

# This Makefile.PL was written by make_ext.pl.
# It will be deleted automatically by make realclean

use strict;
use ExtUtils::MakeMaker;

WriteMakefile(
    NAME          => 'PerlIO::scalar',
    VERSION_FROM  => 'scalar.pm',
    ABSTRACT_FROM => 'scalar.pm',
    realclean     => {FILES => 'Makefile.PL'},
);

# ex: set ro:
