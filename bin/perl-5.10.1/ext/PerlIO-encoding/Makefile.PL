#-*- buffer-read-only: t -*-

# This Makefile.PL was written by make_ext.pl.
# It will be deleted automatically by make realclean

use strict;
use ExtUtils::MakeMaker;

WriteMakefile(
    NAME          => 'PerlIO::encoding',
    VERSION_FROM  => 'encoding.pm',
    ABSTRACT_FROM => 'encoding.pm',
    realclean     => {FILES => 'Makefile.PL'},
);

# ex: set ro:
