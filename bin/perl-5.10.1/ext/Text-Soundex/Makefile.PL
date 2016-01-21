#-*- buffer-read-only: t -*-

# This Makefile.PL was written by make_ext.pl.
# It will be deleted automatically by make realclean

use strict;
use ExtUtils::MakeMaker;

WriteMakefile(
    NAME          => 'Text::Soundex',
    VERSION_FROM  => 'Soundex.pm',
    ABSTRACT_FROM => 'Soundex.pm',
    realclean     => {FILES => 'Makefile.PL'},
);

# ex: set ro:
