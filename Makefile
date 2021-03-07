all: unfold

check:
ifndef QE_ROOT
	$(error QE_ROOT is not defined, compile with "make QE_ROOT=/path/to/qe")
endif

unfold: check
	if test -d src ; then \
        ( cd src ; $(MAKE) QE_ROOT=$(realpath $(QE_ROOT)) || exit 1) ; fi

clean:
	cd src && $(MAKE) clean

