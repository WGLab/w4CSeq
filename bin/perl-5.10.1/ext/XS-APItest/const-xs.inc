BOOT:
  {
#ifdef dTHX
    dTHX;
#endif
    HV *symbol_table = get_hv("XS::APItest::", GV_ADD);
#ifndef SYMBIAN
    HV *constant_missing;
#endif

    static const struct notfound_s values_for_notfound[] =
      {
#ifdef HV_DELETE
#else
        { "HV_DELETE", 9 },
#endif
#ifdef HV_DISABLE_UVAR_XKEY
#else
        { "HV_DISABLE_UVAR_XKEY", 20 },
#endif
#ifdef HV_FETCH_ISSTORE
#else
        { "HV_FETCH_ISSTORE", 16 },
#endif
#ifdef HV_FETCH_ISEXISTS
#else
        { "HV_FETCH_ISEXISTS", 17 },
#endif
#ifdef HV_FETCH_LVALUE
#else
        { "HV_FETCH_LVALUE", 15 },
#endif
#ifdef HV_FETCH_JUST_SV
#else
        { "HV_FETCH_JUST_SV", 16 },
#endif
#ifdef G_SCALAR
#else
        { "G_SCALAR", 8 },
#endif
#ifdef G_ARRAY
#else
        { "G_ARRAY", 7 },
#endif
#ifdef G_VOID
#else
        { "G_VOID", 6 },
#endif
#ifdef G_DISCARD
#else
        { "G_DISCARD", 9 },
#endif
#ifdef G_EVAL
#else
        { "G_EVAL", 6 },
#endif
#ifdef G_NOARGS
#else
        { "G_NOARGS", 8 },
#endif
#ifdef G_KEEPERR
#else
        { "G_KEEPERR", 9 },
#endif
#ifdef G_NODEBUG
#else
        { "G_NODEBUG", 9 },
#endif
#ifdef G_METHOD
#else
        { "G_METHOD", 8 },
#endif
#ifdef G_FAKINGEVAL
#else
        { "G_FAKINGEVAL", 12 },
#endif
        { NULL, 0 } };
	const struct notfound_s *value_for_notfound = values_for_notfound;

    static const struct iv_s values_for_iv[] =
      {
#ifdef HV_DELETE
        { "HV_DELETE", 9, HV_DELETE },
#endif
#ifdef HV_DISABLE_UVAR_XKEY
        { "HV_DISABLE_UVAR_XKEY", 20, HV_DISABLE_UVAR_XKEY },
#endif
#ifdef HV_FETCH_ISSTORE
        { "HV_FETCH_ISSTORE", 16, HV_FETCH_ISSTORE },
#endif
#ifdef HV_FETCH_ISEXISTS
        { "HV_FETCH_ISEXISTS", 17, HV_FETCH_ISEXISTS },
#endif
#ifdef HV_FETCH_LVALUE
        { "HV_FETCH_LVALUE", 15, HV_FETCH_LVALUE },
#endif
#ifdef HV_FETCH_JUST_SV
        { "HV_FETCH_JUST_SV", 16, HV_FETCH_JUST_SV },
#endif
#ifdef G_SCALAR
        { "G_SCALAR", 8, G_SCALAR },
#endif
#ifdef G_ARRAY
        { "G_ARRAY", 7, G_ARRAY },
#endif
#ifdef G_VOID
        { "G_VOID", 6, G_VOID },
#endif
#ifdef G_DISCARD
        { "G_DISCARD", 9, G_DISCARD },
#endif
#ifdef G_EVAL
        { "G_EVAL", 6, G_EVAL },
#endif
#ifdef G_NOARGS
        { "G_NOARGS", 8, G_NOARGS },
#endif
#ifdef G_KEEPERR
        { "G_KEEPERR", 9, G_KEEPERR },
#endif
#ifdef G_NODEBUG
        { "G_NODEBUG", 9, G_NODEBUG },
#endif
#ifdef G_METHOD
        { "G_METHOD", 8, G_METHOD },
#endif
#ifdef G_FAKINGEVAL
        { "G_FAKINGEVAL", 12, G_FAKINGEVAL },
#endif
#ifdef G_WANT
        { "G_WANT", 6, G_WANT },
#endif
#ifdef G_WANT
        /* This is the default value: */
#else
        { "G_WANT", 6, G_ARRAY|G_VOID },
#endif
        { NULL, 0, 0 } };
	const struct iv_s *value_for_iv = values_for_iv;
#ifndef SYMBIAN
	constant_missing = get_missing_hash(aTHX);
#endif
        while (value_for_iv->name) {
	    constant_add_symbol(aTHX_  symbol_table, value_for_iv->name,
				value_for_iv->namelen, newSViv(value_for_iv->value));
	    ++value_for_iv;
	}
	while (value_for_notfound->name) {

	    /* Need to add prototypes, else parsing will vary by platform.  */
	    SV **sv = hv_fetch(symbol_table, value_for_notfound->name,
			       value_for_notfound->namelen, TRUE);
	    if (!sv) {
		Perl_croak(aTHX_ 
			   "Couldn't add key '%s' to %%XS::APItest::",
			   value_for_notfound->name);
	    }
	    if (!SvOK(*sv) && SvTYPE(*sv) != SVt_PVGV) {
		/* Nothing was here before, so mark a prototype of ""  */
		sv_setpvn(*sv, "", 0);
	    } else if (SvPOK(*sv) && SvCUR(*sv) == 0) {
		/* There is already a prototype of "" - do nothing  */
	    } else {
		/* Someone has been here before us - have to make a real
		   typeglob.  */
		/* It turns out to be incredibly hard to deal with all the
		   corner cases of sub foo (); and reporting errors correctly,
		   so lets cheat a bit.  Start with a constant subroutine  */
		CV *cv = newCONSTSUB(symbol_table,
				     value_for_notfound->name,
				     &PL_sv_yes);
		/* and then turn it into a non constant declaration only.  */
		SvREFCNT_dec(CvXSUBANY(cv).any_ptr);
		CvCONST_off(cv);
		CvXSUB(cv) = NULL;
		CvXSUBANY(cv).any_ptr = NULL;
	    }
#ifndef SYMBIAN
	    if (!hv_store(constant_missing, value_for_notfound->name,
			  value_for_notfound->namelen, &PL_sv_yes, 0))
		Perl_croak(aTHX_  "Couldn't add key '%s' to missing_hash",
			   value_for_notfound->name);
#endif

	    ++value_for_notfound;
	}
    /* As we've been creating subroutines, we better invalidate any cached
       methods  */
    ++PL_sub_generation;
  }

void
constant(sv)
    PREINIT:
	STRLEN		len;
    INPUT:
	SV *		sv;
        const char *	s = SvPV(sv, len);
    PPCODE:
#ifdef SYMBIAN
	sv = newSVpvf("%"SVf" is not a valid XS::APItest macro", sv);
#else
	HV *constant_missing = get_missing_hash(aTHX);
	if (hv_exists(constant_missing, s, SvUTF8(sv) ? -(I32)len : (I32)len)) {
	    sv = newSVpvf("Your vendor has not defined XS::APItest macro %" SVf
			  ", used", sv);
	} else {
	    sv = newSVpvf("%"SVf" is not a valid XS::APItest macro",
			  sv);
	}
#endif
	PUSHs(sv_2mortal(sv));
