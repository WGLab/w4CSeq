#define PERL_constant_NOTFOUND	1
#define PERL_constant_NOTDEF	2
#define PERL_constant_ISIV	3
#define PERL_constant_ISNO	4
#define PERL_constant_ISNV	5
#define PERL_constant_ISPV	6
#define PERL_constant_ISPVN	7
#define PERL_constant_ISSV	8
#define PERL_constant_ISUNDEF	9
#define PERL_constant_ISUV	10
#define PERL_constant_ISYES	11

#ifndef NVTYPE
typedef double NV; /* 5.6 and later define NVTYPE, and typedef NV to it.  */
#endif
#ifndef aTHX_
#define aTHX_ /* 5.6 or later define this for threading support.  */
#endif
#ifndef pTHX_
#define pTHX_ /* 5.6 or later define this for threading support.  */
#endif

static int
constant (pTHX_ const char *name, STRLEN len, IV *iv_return) {
  /* Initially switch on the length of the name.  */
  /* When generated this function returned values for the list of names given
     in this section of perl code.  Rather than manually editing these functions
     to add or remove constants, which would result in this comment and section
     of code becoming inaccurate, we recommend that you edit this section of
     code, and use it to regenerate a new set of constant functions which you
     then use to replace the originals.

     Regenerate these constant functions by feeding this entire source file to
     perl -x

#!/var/www/html/w4cseq/bin/perl-5.10.1/miniperl -w
use ExtUtils::Constant qw (constant_types C_constant XS_constant);

my $types = {map {($_, 1)} qw(IV)};
my @names = (qw(GDBM_CACHESIZE GDBM_CENTFREE GDBM_COALESCEBLKS GDBM_FAST
	       GDBM_FASTMODE GDBM_INSERT GDBM_NEWDB GDBM_NOLOCK GDBM_OPENMASK
	       GDBM_READER GDBM_REPLACE GDBM_SYNC GDBM_SYNCMODE GDBM_WRCREAT
	       GDBM_WRITER));

print constant_types(), "\n"; # macro defs
foreach (C_constant ("GDBM_File", 'constant', 'IV', $types, undef, 8, @names) ) {
    print $_, "\n"; # C constant subs
}
print "\n#### XS Section:\n";
print XS_constant ("GDBM_File", $types);
__END__
   */

  switch (len) {
  case 9:
    /* Names all of length 9.  */
    /* GDBM_FAST GDBM_SYNC */
    /* Offset 7 gives the best switch position.  */
    switch (name[7]) {
    case 'N':
      if (memEQ(name, "GDBM_SYNC", 9)) {
      /*                      ^       */
#ifdef GDBM_SYNC
        *iv_return = GDBM_SYNC;
        return PERL_constant_ISIV;
#else
        return PERL_constant_NOTDEF;
#endif
      }
      break;
    case 'S':
      if (memEQ(name, "GDBM_FAST", 9)) {
      /*                      ^       */
#ifdef GDBM_FAST
        *iv_return = GDBM_FAST;
        return PERL_constant_ISIV;
#else
        return PERL_constant_NOTDEF;
#endif
      }
      break;
    }
    break;
  case 10:
    if (memEQ(name, "GDBM_NEWDB", 10)) {
#ifdef GDBM_NEWDB
      *iv_return = GDBM_NEWDB;
      return PERL_constant_ISIV;
#else
      return PERL_constant_NOTDEF;
#endif
    }
    break;
  case 11:
    /* Names all of length 11.  */
    /* GDBM_INSERT GDBM_NOLOCK GDBM_READER GDBM_WRITER */
    /* Offset 6 gives the best switch position.  */
    switch (name[6]) {
    case 'E':
      if (memEQ(name, "GDBM_READER", 11)) {
      /*                     ^           */
#ifdef GDBM_READER
        *iv_return = GDBM_READER;
        return PERL_constant_ISIV;
#else
        return PERL_constant_NOTDEF;
#endif
      }
      break;
    case 'N':
      if (memEQ(name, "GDBM_INSERT", 11)) {
      /*                     ^           */
#ifdef GDBM_INSERT
        *iv_return = GDBM_INSERT;
        return PERL_constant_ISIV;
#else
        return PERL_constant_NOTDEF;
#endif
      }
      break;
    case 'O':
      if (memEQ(name, "GDBM_NOLOCK", 11)) {
      /*                     ^           */
#ifdef GDBM_NOLOCK
        *iv_return = GDBM_NOLOCK;
        return PERL_constant_ISIV;
#else
        return PERL_constant_NOTDEF;
#endif
      }
      break;
    case 'R':
      if (memEQ(name, "GDBM_WRITER", 11)) {
      /*                     ^           */
#ifdef GDBM_WRITER
        *iv_return = GDBM_WRITER;
        return PERL_constant_ISIV;
#else
        return PERL_constant_NOTDEF;
#endif
      }
      break;
    }
    break;
  case 12:
    /* Names all of length 12.  */
    /* GDBM_REPLACE GDBM_WRCREAT */
    /* Offset 10 gives the best switch position.  */
    switch (name[10]) {
    case 'A':
      if (memEQ(name, "GDBM_WRCREAT", 12)) {
      /*                         ^        */
#ifdef GDBM_WRCREAT
        *iv_return = GDBM_WRCREAT;
        return PERL_constant_ISIV;
#else
        return PERL_constant_NOTDEF;
#endif
      }
      break;
    case 'C':
      if (memEQ(name, "GDBM_REPLACE", 12)) {
      /*                         ^        */
#ifdef GDBM_REPLACE
        *iv_return = GDBM_REPLACE;
        return PERL_constant_ISIV;
#else
        return PERL_constant_NOTDEF;
#endif
      }
      break;
    }
    break;
  case 13:
    /* Names all of length 13.  */
    /* GDBM_CENTFREE GDBM_FASTMODE GDBM_OPENMASK GDBM_SYNCMODE */
    /* Offset 5 gives the best switch position.  */
    switch (name[5]) {
    case 'C':
      if (memEQ(name, "GDBM_CENTFREE", 13)) {
      /*                    ^              */
#ifdef GDBM_CENTFREE
        *iv_return = GDBM_CENTFREE;
        return PERL_constant_ISIV;
#else
        return PERL_constant_NOTDEF;
#endif
      }
      break;
    case 'F':
      if (memEQ(name, "GDBM_FASTMODE", 13)) {
      /*                    ^              */
#ifdef GDBM_FASTMODE
        *iv_return = GDBM_FASTMODE;
        return PERL_constant_ISIV;
#else
        return PERL_constant_NOTDEF;
#endif
      }
      break;
    case 'O':
      if (memEQ(name, "GDBM_OPENMASK", 13)) {
      /*                    ^              */
#ifdef GDBM_OPENMASK
        *iv_return = GDBM_OPENMASK;
        return PERL_constant_ISIV;
#else
        return PERL_constant_NOTDEF;
#endif
      }
      break;
    case 'S':
      if (memEQ(name, "GDBM_SYNCMODE", 13)) {
      /*                    ^              */
#ifdef GDBM_SYNCMODE
        *iv_return = GDBM_SYNCMODE;
        return PERL_constant_ISIV;
#else
        return PERL_constant_NOTDEF;
#endif
      }
      break;
    }
    break;
  case 14:
    if (memEQ(name, "GDBM_CACHESIZE", 14)) {
#ifdef GDBM_CACHESIZE
      *iv_return = GDBM_CACHESIZE;
      return PERL_constant_ISIV;
#else
      return PERL_constant_NOTDEF;
#endif
    }
    break;
  case 17:
    if (memEQ(name, "GDBM_COALESCEBLKS", 17)) {
#ifdef GDBM_COALESCEBLKS
      *iv_return = GDBM_COALESCEBLKS;
      return PERL_constant_ISIV;
#else
      return PERL_constant_NOTDEF;
#endif
    }
    break;
  }
  return PERL_constant_NOTFOUND;
}

