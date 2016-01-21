#
# ${R_HOME}/src/library/grDevices/src/cairo/Makefile


srcdir = .
top_srcdir = ../../../../..

top_builddir = ../../../../..
subdir = src/library/grDevices/src/cairo
R_HOME = $(top_builddir)

pkg = grDevices

include $(top_builddir)/Makeconf

## need to find bitmap.h
ALL_CPPFLAGS = $(R_XTRA_CPPFLAGS) -I/usr/include/cairo -I/usr/include/pixman-1 -I/usr/include/freetype2 -I/usr/include/libpng12   \
  -I$(top_srcdir)/src/modules/X11 $(CPPFLAGS) $(DEFS)
ALL_CFLAGS = $(ALL_CFLAGS_LO)

## the X11 module is not made at all if BUILD_X11 is false.
SOURCES = cairoBM.c
DEPENDS = $(SOURCES:.c=.d)
OBJECTS = $(SOURCES:.c=.o) $(top_builddir)/src/modules/X11/rbitmap.o

cairo_la = cairo$(SHLIB_EXT)
## AIX needs LIBM
cairo_la_LIBADD = -lpng12 -lfreetype -lfontconfig -lz -lcairo    -lpng -lz $(LIBR) $(LIBM) # $(LIBINTL)

all: Makedeps
	@$(MAKE) R

Makefile: $(srcdir)/Makefile.in  $(top_builddir)/config.status
	@cd $(top_builddir) && $(SHELL) ./config.status $(subdir)/$@

Makedeps: Makefile $(DEPENDS)
	@cat $(DEPENDS) >> Makefile
	@touch $@

R: Makefile
	@$(MAKE) $(cairo_la)
	@$(MKINSTALLDIRS) $(top_builddir)/library/$(pkg)/libs
	@cp $(cairo_la) $(top_builddir)/library/$(pkg)/libs

$(top_builddir)/src/modules/X11/rbitmap.o:
	(cd $(top_builddir)/src/modules/X11; $(MAKE) rbitmap.o)

$(cairo_la): $(OBJECTS)
	$(SHLIB_LINK) -o $@ $(OBJECTS) $(cairo_la_LIBADD)

mostlyclean: clean
clean:
	@-rm -rf .libs _libs Makedeps *.d *.o $(cairo_la)
distclean: clean
	@-rm -f Makefile
maintainer-clean: distclean

## Automagically generated dependencies:
cairoBM.o: cairoBM.c ../../../../../src/include/config.h \
 ../../../../../src/include/Defn.h \
 ../../../../../src/include/R_ext/Complex.h \
 ../../../../../src/include/Rinternals.h \
 ../../../../../src/include/R_ext/Arith.h \
 ../../../../../src/include/R_ext/libextern.h \
 ../../../../../src/include/R_ext/Boolean.h \
 ../../../../../src/include/R_ext/Error.h \
 ../../../../../src/include/R_ext/Memory.h \
 ../../../../../src/include/R_ext/Utils.h \
 ../../../../../src/include/R_ext/Print.h \
 ../../../../../src/include/Rinlinedfuns.h \
 ../../../../../src/include/Errormsg.h \
 ../../../../../src/include/R_ext/GraphicsEngine.h \
 ../../../../../src/include/R_ext/GraphicsDevice.h \
 ../../../../../src/include/Fileio.h cairoBM.h /usr/include/cairo/cairo.h \
 /usr/include/cairo/cairo-version.h /usr/include/cairo/cairo-features.h \
 /usr/include/cairo/cairo-deprecated.h /usr/include/cairo/cairo-svg.h \
 /usr/include/cairo/cairo.h /usr/include/cairo/cairo-pdf.h \
 /usr/include/cairo/cairo-ps.h cairoFns.c \
 ../../../../../src/modules/X11/bitmap.h
