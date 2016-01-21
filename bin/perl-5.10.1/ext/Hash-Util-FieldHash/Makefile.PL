#-*- buffer-read-only: t -*-

# This Makefile.PL was written by make_ext.pl.
# It will be deleted automatically by make realclean

use strict;
use ExtUtils::MakeMaker;

WriteMakefile(
    NAME          => 'Hash::Util::FieldHash',
    VERSION_FROM  => 'lib/Hash/Util/FieldHash.pm',
    ABSTRACT_FROM => 'lib/Hash/Util/FieldHash.pm',
    realclean     => {FILES => 'Makefile.PL'},
);

# ex: set ro:
