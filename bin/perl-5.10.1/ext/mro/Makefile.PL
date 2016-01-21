#-*- buffer-read-only: t -*-

# This Makefile.PL was written by make_ext.pl.
# It will be deleted automatically by make realclean

use strict;
use ExtUtils::MakeMaker;

WriteMakefile(
    NAME          => 'mro',
    VERSION_FROM  => 'mro.pm',
    ABSTRACT_FROM => 'mro.pm',
    realclean     => {FILES => 'Makefile.PL'},
);

# ex: set ro:
