#-*- buffer-read-only: t -*-

# This Makefile.PL was written by make_ext.pl.
# It will be deleted automatically by make realclean

use strict;
use ExtUtils::MakeMaker;

WriteMakefile(
    NAME          => 'PerlIO::via',
    VERSION_FROM  => 'via.pm',
    ABSTRACT_FROM => 'via.pm',
    realclean     => {FILES => 'Makefile.PL'},
);

# ex: set ro:
