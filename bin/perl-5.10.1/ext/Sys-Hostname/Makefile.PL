#-*- buffer-read-only: t -*-

# This Makefile.PL was written by make_ext.pl.
# It will be deleted automatically by make realclean

use strict;
use ExtUtils::MakeMaker;

WriteMakefile(
    NAME          => 'Sys::Hostname',
    VERSION_FROM  => 'Hostname.pm',
    ABSTRACT_FROM => 'Hostname.pm',
    realclean     => {FILES => 'Makefile.PL'},
);

# ex: set ro:
