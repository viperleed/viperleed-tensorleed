# Build the TensErLEED programs ref-calc (executable refcalc) and rfactor.
#
#   make [all|refcalc|rfactor|debug|clean] [PARAM=<file>] [FC=...]
#        [FFLAGS=...] [BLAS=...] [LDLIBS=...] [BUILDDIR=...] [PROGRAMS=...]
#   make print-<VARIABLE>     print the value of a variable, e.g. REFCALC_SOURCES
#
# BLAS are the BLAS/LAPACK link flags (default: system BLAS and LAPACK;
# e.g. BLAS=-lopenblas); LDLIBS, all libraries refcalc links (default:
# BLAS). rfactor uses neither BLAS nor LAPACK and links no libraries.
#
# PARAM is the dimension file the fixed-form ref-calc sources still INCLUDE
# (required for refcalc; rfactor sizes its arrays from its input and
# ignores it); it is copied to BUILDDIR, which also holds objects, module
# files and executables, so the source tree is never written to. 'debug'
# builds PROGRAMS (default: both) with DEBUG_FFLAGS into BUILDDIR/debug,
# e.g. make debug PROGRAMS=rfactor
#
# Each program's free-form sources (REFCALC_F90, RFACTOR_F90) are listed in
# dependency order, modules first, and compile in that order with
# -std=f2018; their objects go to BUILDDIR/obj/<program> and their module
# files to BUILDDIR/mod/<program>, so modules of the same name in the two
# programs do not clash. refcalc links them with the fixed-form
# src/ref-calc.f and lib/lib.tleed.f, which keep their legacy flags and are
# exempt from -std. <PROGRAM>_SOURCES lists all sources of a program.

FC = gfortran
FFLAGS = -O2
BLAS = -llapack -lpthread -lblas
LDLIBS = $(BLAS)
BUILDDIR = build
DEBUG_FFLAGS = -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -Wall
STD_FFLAGS = -std=f2018

B := $(BUILDDIR)

REFCALC_F90 = lib/kinds.f90
REFCALC_SOURCES = $(REFCALC_F90) src/ref-calc.f lib/lib.tleed.f
RFACTOR_F90 = lib/rfacsb.f90 src/rfactor.f90
RFACTOR_SOURCES = $(RFACTOR_F90)

vpath %.f90 lib src

# objects of the free-form sources $(2) of program $(1)
f90_objects = $(patsubst %.f90,$(B)/obj/$(1)/%.o,$(notdir $(2)))
REFCALC_F90_OBJ := $(call f90_objects,refcalc,$(REFCALC_F90))
RFACTOR_F90_OBJ := $(call f90_objects,rfactor,$(RFACTOR_F90))

# each free-form object depends on the one listed before it
chain = $(if $(word 2,$(1)),$(eval $(word 2,$(1)): $(firstword $(1)))$(call chain,$(wordlist 2,$(words $(1)),$(1))))
$(call chain,$(REFCALC_F90_OBJ))
$(call chain,$(RFACTOR_F90_OBJ))

.PHONY: all refcalc rfactor debug clean

PROGRAMS = refcalc rfactor

all: $(PROGRAMS)

refcalc: $(B)/refcalc

rfactor: $(B)/rfactor

debug:
	$(MAKE) $(PROGRAMS) FFLAGS='$(DEBUG_FFLAGS)' BUILDDIR='$(BUILDDIR)/debug'

print-%:
	@echo '$($*)'

$(B) $(B)/obj/refcalc $(B)/obj/rfactor $(B)/mod/refcalc $(B)/mod/rfactor:
	mkdir -p $@

$(B)/PARAM: $(PARAM) | $(B)
	@test -n "$(PARAM)" || { echo "PARAM=<file> is required" >&2; exit 1; }
	cp $(PARAM) $@

$(B)/obj/refcalc/%.o: %.f90 | $(B)/obj/refcalc $(B)/mod/refcalc
	$(FC) $(FFLAGS) $(STD_FFLAGS) -J$(B)/mod/refcalc -o $@ -c $<

$(B)/obj/rfactor/%.o: %.f90 | $(B)/obj/rfactor $(B)/mod/rfactor
	$(FC) $(FFLAGS) $(STD_FFLAGS) -J$(B)/mod/rfactor -o $@ -c $<

# ref-calc: the fixed-form objects use the modules, so they depend on them
# for real (an incremental make rebuilds them when a module changes)
$(B)/ref-calc.o: src/ref-calc.f src/GLOBAL $(B)/PARAM $(REFCALC_F90_OBJ)
	$(FC) $(FFLAGS) -I$(B) -I$(B)/mod/refcalc -o $@ -c $<

$(B)/lib.tleed.o: lib/lib.tleed.f src/GLOBAL $(REFCALC_F90_OBJ) | $(B)
	$(FC) $(FFLAGS) -Isrc -I$(B) -I$(B)/mod/refcalc -o $@ -c $<

$(B)/refcalc: $(REFCALC_F90_OBJ) $(B)/ref-calc.o $(B)/lib.tleed.o
	$(FC) $(FFLAGS) -o $@ $^ $(LDLIBS)

# rfactor
$(B)/rfactor: $(RFACTOR_F90_OBJ)
	$(FC) $(FFLAGS) -o $@ $^

clean:
	rm -f $(B)/*.o $(B)/PARAM $(B)/refcalc $(B)/rfactor \
	      $(B)/obj/*/*.o $(B)/mod/*/*.mod
	rm -f $(B)/debug/*.o $(B)/debug/PARAM $(B)/debug/refcalc $(B)/debug/rfactor \
	      $(B)/debug/obj/*/*.o $(B)/debug/mod/*/*.mod
