# Makefile for Unfold.x

-include $(QE_ROOT)/make.inc

# location of needed modules and included files (if any)
MODFLAGS= $(MOD_FLAG)$(QE_ROOT)/UtilXlib \
	  $(MOD_FLAG)$(QE_ROOT)/FoX/finclude \
	  $(MOD_FLAG)$(QE_ROOT)/Modules \
	  $(MOD_FLAG)$(QE_ROOT)/upflib \
	  $(MOD_FLAG)$(QE_ROOT)/PW/src $(MOD_FLAG).
IFLAGS=

#location of needed libraries
LIBOBJS= $(QE_ROOT)/clib/clib.a

UNOBJS = \
unfold.o

UNKLST = \
unklist.o

UNLIBS = \
openfil_pp.o \
unfold_data.o \
pkm.o \
vecfbz.o \
specfun.o \
wbin.o \
wxml.o \
unfold_input_parameters.o \
unfold_read_cards.o


QEMODS=$(QE_ROOT)/Modules/libqemod.a \
       $(QE_ROOT)/KS_Solvers/libks_solvers.a \
       $(QE_ROOT)/FFTXlib/libqefft.a \
       $(QE_ROOT)/LAXlib/libqela.a \
       $(QE_ROOT)/UtilXlib/libutil.a \
       $(QE_ROOT)/upflib/libupf.a \
       $(QE_ROOT)/dft-d3/libdftd3qe.a

PWOBJS= $(QE_ROOT)/PW/src/libpw.a

TLDEPS=bindir mods libs liblapack pw

all : tldeps unfold.x unklist.x

unfold.x : $(UNOBJS) libunfold.a $(LIBOBJS) $(PWOBJS) $(QEMODS)
	$(LD) $(LDFLAGS) -g -o $@ \
	$(UNOBJS) libunfold.a $(PWOBJS) $(QEMODS) $(LIBOBJS) $(QELIBS)
	- ( mkdir -p ../bin ; cd ../bin; ln -fs ../src/$@ . )

unklist.x : $(UNKLST) libunfold.a $(PWOBJS) $(LIBOBJS) $(QEMODS)
	$(LD) $(LDFLAGS) -o $@ \
	unklist.o libunfold.a $(PWOBJS) $(QEMODS) $(LIBOBJS) $(QELIBS)
	- ( mkdir -p ../bin ; cd ../bin; ln -fs ../src/$@ . )

libunfold.a : $(UNLIBS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

tldeps :
	if test -n "$(TLDEPS)" ; then \
	( cd $(QE_ROOT) ; $(MAKE) $(TLDEPS) || exit 1 ) ; fi

clean :
	- /bin/rm -f *.x *.o *.a *~ *.F90 *.d *.mod *.i *.L

include make.depend
