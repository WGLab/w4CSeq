#-*- buffer-read-only: t -*-

# This Makefile.PL was written by make_ext.pl.
# It will be deleted automatically by make realclean

use strict;
use ExtUtils::MakeMaker;

WriteMakefile(
    NAME          => 'Filter::Util::Call',
    VERSION_FROM  => 'Call.pm',
    ABSTRACT_FROM => 'Call.pm',
    realclean     => {FILES => 'Makefile.PL'},
);

# ex: set ro:
