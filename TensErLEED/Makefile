# Build the TensErLEED programs ref-calc (executable refcalc) and rfactor.
#
#   make [all|refcalc|rfactor|debug|clean] PARAM=<file> [FC=...] [FFLAGS=...]
#        [LDLIBS=...] [BUILDDIR=...] [PROGRAMS=...]
#
# PARAM is the dimension file the legacy sources still INCLUDE; it is copied
# to BUILDDIR, which also holds objects, module files and executables, so the
# source tree is never written to. The programs need different PARAM files,
# so build one per call (make rfactor PARAM=...). 'debug' builds PROGRAMS
# (default: both) with DEBUG_FFLAGS into BUILDDIR/debug, e.g.
# make debug PROGRAMS=rfactor PARAM=...
# Free-form sources compile with -std=f2018; the legacy fixed-form sources
# are exempt.

FC = gfortran
FFLAGS = -O2
LDLIBS = -llapack -lpthread -lblas
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

$(B)/rfacsb.o: lib/rfacsb.f $(B)/rfacsb.f90.o
	$(FC) $(FFLAGS) -I$(B) -o $@ -c $<

$(B)/rfactor.o: src/rfactor.f $(B)/PARAM $(B)/rfacsb.f90.o
	$(FC) $(FFLAGS) -I$(B) -o $@ -c $<

$(B)/rfactor: $(B)/rfacsb.f90.o $(B)/rfacsb.o $(B)/rfactor.o
	$(FC) $(FFLAGS) -o $@ $^ $(LDLIBS)

clean:
	rm -f $(B)/*.o $(B)/*.mod $(B)/PARAM $(B)/refcalc $(B)/rfactor
	rm -f $(B)/debug/*.o $(B)/debug/*.mod $(B)/debug/PARAM $(B)/debug/refcalc \
	      $(B)/debug/rfactor
