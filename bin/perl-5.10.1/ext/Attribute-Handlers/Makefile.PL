#-*- buffer-read-only: t -*-

# This Makefile.PL was written by make_ext.pl.
# It will be deleted automatically by make realclean

use strict;
use ExtUtils::MakeMaker;

WriteMakefile(
    NAME          => 'Attribute::Handlers',
    VERSION_FROM  => 'lib/Attribute/Handlers.pm',
    ABSTRACT_FROM => 'lib/Attribute/Handlers.pm',
    realclean     => {FILES => 'Makefile.PL'},
);

# ex: set ro:
