#-*- buffer-read-only: t -*-

# This Makefile.PL was written by make_ext.pl.
# It will be deleted automatically by make realclean

use strict;
use ExtUtils::MakeMaker;

WriteMakefile(
    NAME          => 'Safe',
    VERSION_FROM  => 'Safe.pm',
    ABSTRACT_FROM => 'Safe.pm',
    realclean     => {FILES => 'Makefile.PL'},
);

# ex: set ro:
