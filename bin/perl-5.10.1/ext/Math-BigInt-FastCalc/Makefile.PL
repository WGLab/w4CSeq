#-*- buffer-read-only: t -*-

# This Makefile.PL was written by make_ext.pl.
# It will be deleted automatically by make realclean

use strict;
use ExtUtils::MakeMaker;

WriteMakefile(
    NAME          => 'Math::BigInt::FastCalc',
    VERSION_FROM  => 'FastCalc.pm',
    ABSTRACT_FROM => 'FastCalc.pm',
    realclean     => {FILES => 'Makefile.PL'},
);

# ex: set ro:
