cd ext/DynaLoader
if test ! -f Makefile -a -f Makefile.old; then
    echo "Note: Using Makefile.old"
    make -f Makefile.old realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=static CCCDLFLAGS=
else
    if test ! -f Makefile ; then
	echo "Warning: No Makefile!"
    fi
    make realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=static CCCDLFLAGS=
fi
cd ../..
cd ext/attrs
if test ! -f Makefile -a -f Makefile.old; then
    echo "Note: Using Makefile.old"
    make -f Makefile.old realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
else
    if test ! -f Makefile ; then
	echo "Warning: No Makefile!"
    fi
    make realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
fi
cd ../..
cd ext/B
if test ! -f Makefile -a -f Makefile.old; then
    echo "Note: Using Makefile.old"
    make -f Makefile.old realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
else
    if test ! -f Makefile ; then
	echo "Warning: No Makefile!"
    fi
    make realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
fi
cd ../..
cd ext/Compress-Raw-Bzip2
if test ! -f Makefile -a -f Makefile.old; then
    echo "Note: Using Makefile.old"
    make -f Makefile.old realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
else
    if test ! -f Makefile ; then
	echo "Warning: No Makefile!"
    fi
    make realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
fi
cd ../..
cd ext/Compress-Raw-Zlib
if test ! -f Makefile -a -f Makefile.old; then
    echo "Note: Using Makefile.old"
    make -f Makefile.old realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
else
    if test ! -f Makefile ; then
	echo "Warning: No Makefile!"
    fi
    make realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
fi
cd ../..
cd ext/Cwd
if test ! -f Makefile -a -f Makefile.old; then
    echo "Note: Using Makefile.old"
    make -f Makefile.old realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
else
    if test ! -f Makefile ; then
	echo "Warning: No Makefile!"
    fi
    make realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
fi
cd ../..
cd ext/Data-Dumper
if test ! -f Makefile -a -f Makefile.old; then
    echo "Note: Using Makefile.old"
    make -f Makefile.old realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
else
    if test ! -f Makefile ; then
	echo "Warning: No Makefile!"
    fi
    make realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
fi
cd ../..
cd ext/DB_File
if test ! -f Makefile -a -f Makefile.old; then
    echo "Note: Using Makefile.old"
    make -f Makefile.old realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
else
    if test ! -f Makefile ; then
	echo "Warning: No Makefile!"
    fi
    make realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
fi
cd ../..
cd ext/Devel-DProf
if test ! -f Makefile -a -f Makefile.old; then
    echo "Note: Using Makefile.old"
    make -f Makefile.old realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
else
    if test ! -f Makefile ; then
	echo "Warning: No Makefile!"
    fi
    make realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
fi
cd ../..
cd ext/Devel-Peek
if test ! -f Makefile -a -f Makefile.old; then
    echo "Note: Using Makefile.old"
    make -f Makefile.old realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
else
    if test ! -f Makefile ; then
	echo "Warning: No Makefile!"
    fi
    make realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
fi
cd ../..
cd ext/Devel-PPPort
if test ! -f Makefile -a -f Makefile.old; then
    echo "Note: Using Makefile.old"
    make -f Makefile.old realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
else
    if test ! -f Makefile ; then
	echo "Warning: No Makefile!"
    fi
    make realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
fi
cd ../..
cd ext/Digest-MD5
if test ! -f Makefile -a -f Makefile.old; then
    echo "Note: Using Makefile.old"
    make -f Makefile.old realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
else
    if test ! -f Makefile ; then
	echo "Warning: No Makefile!"
    fi
    make realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
fi
cd ../..
cd ext/Digest-SHA
if test ! -f Makefile -a -f Makefile.old; then
    echo "Note: Using Makefile.old"
    make -f Makefile.old realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
else
    if test ! -f Makefile ; then
	echo "Warning: No Makefile!"
    fi
    make realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
fi
cd ../..
cd ext/Encode
if test ! -f Makefile -a -f Makefile.old; then
    echo "Note: Using Makefile.old"
    make -f Makefile.old realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
else
    if test ! -f Makefile ; then
	echo "Warning: No Makefile!"
    fi
    make realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
fi
cd ../..
cd ext/Fcntl
if test ! -f Makefile -a -f Makefile.old; then
    echo "Note: Using Makefile.old"
    make -f Makefile.old realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
else
    if test ! -f Makefile ; then
	echo "Warning: No Makefile!"
    fi
    make realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
fi
cd ../..
cd ext/File-Glob
if test ! -f Makefile -a -f Makefile.old; then
    echo "Note: Using Makefile.old"
    make -f Makefile.old realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
else
    if test ! -f Makefile ; then
	echo "Warning: No Makefile!"
    fi
    make realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
fi
cd ../..
cd ext/Filter-Util-Call
if test ! -f Makefile -a -f Makefile.old; then
    echo "Note: Using Makefile.old"
    make -f Makefile.old realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
else
    if test ! -f Makefile ; then
	echo "Warning: No Makefile!"
    fi
    make realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
fi
cd ../..
cd ext/GDBM_File
if test ! -f Makefile -a -f Makefile.old; then
    echo "Note: Using Makefile.old"
    make -f Makefile.old realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
else
    if test ! -f Makefile ; then
	echo "Warning: No Makefile!"
    fi
    make realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
fi
cd ../..
cd ext/Hash-Util
if test ! -f Makefile -a -f Makefile.old; then
    echo "Note: Using Makefile.old"
    make -f Makefile.old realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
else
    if test ! -f Makefile ; then
	echo "Warning: No Makefile!"
    fi
    make realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
fi
cd ../..
cd ext/Hash-Util-FieldHash
if test ! -f Makefile -a -f Makefile.old; then
    echo "Note: Using Makefile.old"
    make -f Makefile.old realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
else
    if test ! -f Makefile ; then
	echo "Warning: No Makefile!"
    fi
    make realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
fi
cd ../..
cd ext/I18N-Langinfo
if test ! -f Makefile -a -f Makefile.old; then
    echo "Note: Using Makefile.old"
    make -f Makefile.old realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
else
    if test ! -f Makefile ; then
	echo "Warning: No Makefile!"
    fi
    make realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
fi
cd ../..
cd ext/IO
if test ! -f Makefile -a -f Makefile.old; then
    echo "Note: Using Makefile.old"
    make -f Makefile.old realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
else
    if test ! -f Makefile ; then
	echo "Warning: No Makefile!"
    fi
    make realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
fi
cd ../..
cd ext/IO-Compress
if test ! -f Makefile -a -f Makefile.old; then
    echo "Note: Using Makefile.old"
    make -f Makefile.old realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
else
    if test ! -f Makefile ; then
	echo "Warning: No Makefile!"
    fi
    make realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
fi
cd ../..
cd ext/IPC-SysV
if test ! -f Makefile -a -f Makefile.old; then
    echo "Note: Using Makefile.old"
    make -f Makefile.old realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
else
    if test ! -f Makefile ; then
	echo "Warning: No Makefile!"
    fi
    make realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
fi
cd ../..
cd ext/List-Util
if test ! -f Makefile -a -f Makefile.old; then
    echo "Note: Using Makefile.old"
    make -f Makefile.old realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
else
    if test ! -f Makefile ; then
	echo "Warning: No Makefile!"
    fi
    make realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
fi
cd ../..
cd ext/Math-BigInt-FastCalc
if test ! -f Makefile -a -f Makefile.old; then
    echo "Note: Using Makefile.old"
    make -f Makefile.old realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
else
    if test ! -f Makefile ; then
	echo "Warning: No Makefile!"
    fi
    make realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
fi
cd ../..
cd ext/MIME-Base64
if test ! -f Makefile -a -f Makefile.old; then
    echo "Note: Using Makefile.old"
    make -f Makefile.old realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
else
    if test ! -f Makefile ; then
	echo "Warning: No Makefile!"
    fi
    make realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
fi
cd ../..
cd ext/mro
if test ! -f Makefile -a -f Makefile.old; then
    echo "Note: Using Makefile.old"
    make -f Makefile.old realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
else
    if test ! -f Makefile ; then
	echo "Warning: No Makefile!"
    fi
    make realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
fi
cd ../..
cd ext/NDBM_File
if test ! -f Makefile -a -f Makefile.old; then
    echo "Note: Using Makefile.old"
    make -f Makefile.old realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
else
    if test ! -f Makefile ; then
	echo "Warning: No Makefile!"
    fi
    make realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
fi
cd ../..
cd ext/Opcode
if test ! -f Makefile -a -f Makefile.old; then
    echo "Note: Using Makefile.old"
    make -f Makefile.old realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
else
    if test ! -f Makefile ; then
	echo "Warning: No Makefile!"
    fi
    make realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
fi
cd ../..
cd ext/PerlIO-encoding
if test ! -f Makefile -a -f Makefile.old; then
    echo "Note: Using Makefile.old"
    make -f Makefile.old realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
else
    if test ! -f Makefile ; then
	echo "Warning: No Makefile!"
    fi
    make realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
fi
cd ../..
cd ext/PerlIO-scalar
if test ! -f Makefile -a -f Makefile.old; then
    echo "Note: Using Makefile.old"
    make -f Makefile.old realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
else
    if test ! -f Makefile ; then
	echo "Warning: No Makefile!"
    fi
    make realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
fi
cd ../..
cd ext/PerlIO-via
if test ! -f Makefile -a -f Makefile.old; then
    echo "Note: Using Makefile.old"
    make -f Makefile.old realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
else
    if test ! -f Makefile ; then
	echo "Warning: No Makefile!"
    fi
    make realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
fi
cd ../..
cd ext/POSIX
if test ! -f Makefile -a -f Makefile.old; then
    echo "Note: Using Makefile.old"
    make -f Makefile.old realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
else
    if test ! -f Makefile ; then
	echo "Warning: No Makefile!"
    fi
    make realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
fi
cd ../..
cd ext/re
if test ! -f Makefile -a -f Makefile.old; then
    echo "Note: Using Makefile.old"
    make -f Makefile.old realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
else
    if test ! -f Makefile ; then
	echo "Warning: No Makefile!"
    fi
    make realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
fi
cd ../..
cd ext/SDBM_File
if test ! -f Makefile -a -f Makefile.old; then
    echo "Note: Using Makefile.old"
    make -f Makefile.old realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
else
    if test ! -f Makefile ; then
	echo "Warning: No Makefile!"
    fi
    make realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
fi
cd ../..
cd ext/Socket
if test ! -f Makefile -a -f Makefile.old; then
    echo "Note: Using Makefile.old"
    make -f Makefile.old realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
else
    if test ! -f Makefile ; then
	echo "Warning: No Makefile!"
    fi
    make realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
fi
cd ../..
cd ext/Storable
if test ! -f Makefile -a -f Makefile.old; then
    echo "Note: Using Makefile.old"
    make -f Makefile.old realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
else
    if test ! -f Makefile ; then
	echo "Warning: No Makefile!"
    fi
    make realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
fi
cd ../..
cd ext/Sys-Hostname
if test ! -f Makefile -a -f Makefile.old; then
    echo "Note: Using Makefile.old"
    make -f Makefile.old realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
else
    if test ! -f Makefile ; then
	echo "Warning: No Makefile!"
    fi
    make realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
fi
cd ../..
cd ext/Sys-Syslog
if test ! -f Makefile -a -f Makefile.old; then
    echo "Note: Using Makefile.old"
    make -f Makefile.old realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
else
    if test ! -f Makefile ; then
	echo "Warning: No Makefile!"
    fi
    make realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
fi
cd ../..
cd ext/Text-Soundex
if test ! -f Makefile -a -f Makefile.old; then
    echo "Note: Using Makefile.old"
    make -f Makefile.old realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
else
    if test ! -f Makefile ; then
	echo "Warning: No Makefile!"
    fi
    make realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
fi
cd ../..
cd ext/threads
if test ! -f Makefile -a -f Makefile.old; then
    echo "Note: Using Makefile.old"
    make -f Makefile.old realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
else
    if test ! -f Makefile ; then
	echo "Warning: No Makefile!"
    fi
    make realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
fi
cd ../..
cd ext/threads-shared
if test ! -f Makefile -a -f Makefile.old; then
    echo "Note: Using Makefile.old"
    make -f Makefile.old realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
else
    if test ! -f Makefile ; then
	echo "Warning: No Makefile!"
    fi
    make realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
fi
cd ../..
cd ext/Time-HiRes
if test ! -f Makefile -a -f Makefile.old; then
    echo "Note: Using Makefile.old"
    make -f Makefile.old realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
else
    if test ! -f Makefile ; then
	echo "Warning: No Makefile!"
    fi
    make realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
fi
cd ../..
cd ext/Time-Piece
if test ! -f Makefile -a -f Makefile.old; then
    echo "Note: Using Makefile.old"
    make -f Makefile.old realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
else
    if test ! -f Makefile ; then
	echo "Warning: No Makefile!"
    fi
    make realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
fi
cd ../..
cd ext/Unicode-Normalize
if test ! -f Makefile -a -f Makefile.old; then
    echo "Note: Using Makefile.old"
    make -f Makefile.old realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
else
    if test ! -f Makefile ; then
	echo "Warning: No Makefile!"
    fi
    make realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
fi
cd ../..
cd ext/XS-APItest
if test ! -f Makefile -a -f Makefile.old; then
    echo "Note: Using Makefile.old"
    make -f Makefile.old realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
else
    if test ! -f Makefile ; then
	echo "Warning: No Makefile!"
    fi
    make realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
fi
cd ../..
cd ext/XS-Typemap
if test ! -f Makefile -a -f Makefile.old; then
    echo "Note: Using Makefile.old"
    make -f Makefile.old realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
else
    if test ! -f Makefile ; then
	echo "Warning: No Makefile!"
    fi
    make realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a LINKTYPE=dynamic
fi
cd ../..
cd ext/Attribute-Handlers
if test ! -f Makefile -a -f Makefile.old; then
    echo "Note: Using Makefile.old"
    make -f Makefile.old realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a
else
    if test ! -f Makefile ; then
	echo "Warning: No Makefile!"
    fi
    make realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a
fi
cd ../..
cd ext/Errno
if test ! -f Makefile -a -f Makefile.old; then
    echo "Note: Using Makefile.old"
    make -f Makefile.old realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a
else
    if test ! -f Makefile ; then
	echo "Warning: No Makefile!"
    fi
    make realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a
fi
cd ../..
cd ext/Module-Pluggable
if test ! -f Makefile -a -f Makefile.old; then
    echo "Note: Using Makefile.old"
    make -f Makefile.old realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a
else
    if test ! -f Makefile ; then
	echo "Warning: No Makefile!"
    fi
    make realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a
fi
cd ../..
cd ext/Safe
if test ! -f Makefile -a -f Makefile.old; then
    echo "Note: Using Makefile.old"
    make -f Makefile.old realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a
else
    if test ! -f Makefile ; then
	echo "Warning: No Makefile!"
    fi
    make realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a
fi
cd ../..
cd ext/Test-Harness
if test ! -f Makefile -a -f Makefile.old; then
    echo "Note: Using Makefile.old"
    make -f Makefile.old realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a
else
    if test ! -f Makefile ; then
	echo "Warning: No Makefile!"
    fi
    make realclean MAKE='make' PERL_CORE=1 LIBPERL_A=libperl.a
fi
cd ../..
