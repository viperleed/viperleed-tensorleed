# Build the TensErLEED programs ref-calc (executable refcalc) and rfactor.
#
#   make [all|refcalc|rfactor|debug|clean] [PARAM=<file>] [FC=...]
#        [FFLAGS=...] [BLAS=...] [LDLIBS=...] [BUILDDIR=...] [PROGRAMS=...]
#
# BLAS are the BLAS/LAPACK link flags (default: system BLAS and LAPACK;
# e.g. BLAS=-lopenblas); LDLIBS, all libraries refcalc links (default:
# BLAS). rfactor uses neither BLAS nor LAPACK and links no libraries.
#
# PARAM is the dimension file the legacy ref-calc sources still INCLUDE
# (required for refcalc; rfactor sizes its arrays from its input and
# ignores it); it is copied to BUILDDIR, which also holds objects, module
# files and executables, so the source tree is never written to. 'debug'
# builds PROGRAMS (default: both) with DEBUG_FFLAGS into BUILDDIR/debug,
# e.g. make debug PROGRAMS=rfactor
# Free-form sources compile with -std=f2018; the legacy fixed-form sources
# are exempt.

FC = gfortran
FFLAGS = -O2
BLAS = -llapack -lpthread -lblas
LDLIBS = $(BLAS)
BUILDDIR = build
DEBUG_FFLAGS = -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -Wall
STD_FFLAGS = -std=f2018

B := $(BUILDDIR)

.PHONY: all refcalc rfactor debug clean

PROGRAMS = refcalc rfactor

all: $(PROGRAMS)

refcalc: $(B)/refcalc

rfactor: $(B)/rfactor

debug:
	$(MAKE) $(PROGRAMS) FFLAGS='$(DEBUG_FFLAGS)' BUILDDIR='$(BUILDDIR)/debug'

$(B):
	mkdir -p $@

$(B)/PARAM: $(PARAM) | $(B)
	@test -n "$(PARAM)" || { echo "PARAM=<file> is required" >&2; exit 1; }
	cp $(PARAM) $@

# ref-calc
$(B)/ref-calc.o: src/ref-calc.f src/GLOBAL $(B)/PARAM
	$(FC) $(FFLAGS) -I$(B) -o $@ -c $<

$(B)/lib.tleed.o: lib/lib.tleed.f src/GLOBAL | $(B)
	$(FC) $(FFLAGS) -Isrc -I$(B) -o $@ -c $<

$(B)/refcalc: $(B)/ref-calc.o $(B)/lib.tleed.o
	$(FC) $(FFLAGS) -o $@ $^ $(LDLIBS)

# rfactor
$(B)/rfacsb.f90.o: lib/rfacsb.f90 | $(B)
	$(FC) $(FFLAGS) $(STD_FFLAGS) -J$(B) -o $@ -c $<

$(B)/rfactor.f90.o: src/rfactor.f90 $(B)/rfacsb.f90.o
	$(FC) $(FFLAGS) $(STD_FFLAGS) -I$(B) -o $@ -c $<

$(B)/rfactor: $(B)/rfacsb.f90.o $(B)/rfactor.f90.o
	$(FC) $(FFLAGS) -o $@ $^

clean:
	rm -f $(B)/*.o $(B)/*.mod $(B)/PARAM $(B)/refcalc $(B)/rfactor
	rm -f $(B)/debug/*.o $(B)/debug/*.mod $(B)/debug/PARAM $(B)/debug/refcalc \
	      $(B)/debug/rfactor
