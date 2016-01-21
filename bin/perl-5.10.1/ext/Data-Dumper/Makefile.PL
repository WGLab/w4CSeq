#-*- buffer-read-only: t -*-

# This Makefile.PL was written by make_ext.pl.
# It will be deleted automatically by make realclean

use strict;
use ExtUtils::MakeMaker;

WriteMakefile(
    NAME          => 'Data::Dumper',
    VERSION_FROM  => 'Dumper.pm',
    ABSTRACT_FROM => 'Dumper.pm',
    realclean     => {FILES => 'Makefile.PL'},
);

# ex: set ro:
