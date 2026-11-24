# ref-calc: design of the free-form port

The reference-calculation program `ref-calc` (executable `refcalc`) is the
main program `LEED` in `src/ref-calc.f` plus the routine library
`lib/lib.tleed.f`, both fixed-form, with array sizes from the `PARAM`
include file. This document fixes where each legacy routine and each COMMON
block goes in the free-form port, the storage-association tricks of the
legacy code that the port must keep or undo on purpose, and the conventions
of the new files.

The port preserves behaviour: `fd.out`, `amp.out`, the tensor files and the
stdout log keep their format and, within the test tolerances, their values;
floating-point operations stay in the same order where results depend on
it; every variable keeps its legacy precision (`sp` for default `REAL` and
`COMPLEX`, `dp` for `REAL*8`/`COMPLEX*16`). Routines are ported one at a
time or as a tightly coupled group, with the regression tests green after
each step.

## 1. Port stages

The port runs in stages; each fills the modules and routines assigned to it
below and switches the driver's calls over to them.

| Stage | Content |
|---|---|
| 1 | Build: Makefile rules for mixed fixed/free-form builds, module `kinds`, this document |
| 2 | Input reading into derived types (`READIN`, `READGEO`) |
| 3 | Energy-independent setup: Clebsch-Gordan and Gaunt coefficients, sublattice set-up |
| 4 | Atomic t-matrices, beam selection, Bravais layers |
| 5 | Composite layers (Beeby inversion) |
| 6 | Bulk (layer doubling) and layer stacking |
| 7 | Top layer, intensities, tensors, output files |
| 8 | Remove the remaining COMMON blocks and `PARAM`; arrays sized from the input, per energy where it saves memory (number of beams in use, LMAX); driver in free form |
| 9 | viperleed: for TensErLEED >= 3.0.0, no `PARAM` and no run-time compilation of `ref-calc`; writes the input-only sizes of section 5 into the input |

Stages 1 to 5 are done: `READIN` and `READGEO` live in `leed_input` as
`read_input` and `read_geometry`, their stdout echo in `leed_output`, and
the driver copies what they return into its COMMON blocks and
`PARAM`-sized arrays for the routines still in `lib.tleed.f`. The
energy-independent setup runs in `angular_coefficients` (`CAA`, `CLM`, the
`(L,M)` permutations), `atomic_scattering` (`PPP`) and `bravais_layer`
(`V`, `VL`, `JJS`, from the `lattice` type); the driver keeps its
`PARAM`-sized arrays for the legacy consumers and passes them to the
ported routines. The constants of `src/GLOBAL` and `pi` live in
`constants`; the fixed-form sources still include `GLOBAL`.

Stage 4 adds the atomic t-matrices to `atomic_scattering`
(`atomic_tmatrix`, with `thermal_phase_shifts` and `spherical_bessel` as
private helpers, and `accumulate_atomic_tmatrix` for the sum over
elements), the beam selection `select_beams` in `beam_selection`, and in
`bravais_layer` `lattice_sums`, `spherical_harmonics`, `beams_ylm` and
`bravais_layer_matrices` (with `intralayer_matrix` and `matrix_elements` as
private helpers). The new type `energy_point` (energy, k-parallel,
imaginary potential) in `leed_types` replaces the blank COMMON of the
routines of this stage; `TV` and `LMAX` stay plain arguments. The routines
keep the precision and packed buffers of the legacy ones (section 4.2) and
reproduce their results bit for bit (`tests/tier1/test_refcalc_routines.py`).

Stage 5 puts the composite layers into `composite_layer`:
`sort_subplanes` (`SRTLAY`), `interplanar_propagators` (`GHMAT`),
`subplane_scattering` (`TAUMAT`), `composite_layer_matrices` (`RTINV_SIMPLE`)
and `copy_subplane_positions` (`COPYPOS`), with the helpers of the legacy
routines as private procedures. `/MPT/` and `/MFB/` are gone: their members
are arguments. `beams_ylm` also serves `composite_layer_matrices` and the
driver's call before `FINAL1`, so the legacy `BEAMS_YLM` and `SPHRM_MOD`
left `lib/lib.tleed.f`. The same tests cover them.

## 2. Target files and modules

All new files are free-form `.f90`, listed in `REFCALC_F90` of the
Makefile in dependency order (modules first).

| File | Module | Contents | Stage |
|---|---|---|---|
| `lib/kinds.f90` | `kinds` | kind parameters `sp`, `dp` | 1 |
| `lib/constants.f90` | `constants` | physical and numerical constants: `bohr`, `hartree`, `memach` (the `PARAMETER`s of `src/GLOBAL`), `pi` (the legacy single-precision literal) | 3 |
| `lib/leed_types.f90` | `leed_types` | derived types replacing the COMMON blocks and the `PARAM` sizes (section 4) | created 2, extended by each later stage |
| `lib/leed_input.f90` | `leed_input` | reader: parses the input on stdin into `leed_types` | 2 |
| `lib/leed_output.f90` | `leed_output` | writer: stdout log, `fd.out`, `amp.out`, tensor files | created 2 (input echo), files 7 |
| `lib/angular_coefficients.f90` | `angular_coefficients` | Clebsch-Gordan/Gaunt coefficients, factorials, index tables | 3 |
| `lib/atomic_scattering.f90` | `atomic_scattering` | phase shifts to temperature-dependent atomic t-matrices | 3 (`CPPP`), 4, done |
| `lib/beam_selection.f90` | `beam_selection` | beams in use at an energy; exit-beam bookkeeping of the tensor loop | 4, done (`select_beams`), 7 (`DECIDE`) |
| `lib/bravais_layer.f90` | `bravais_layer` | lattice sums and reflection/transmission of Bravais layers | 3 (`SLIND`), 4, done |
| `lib/composite_layer.f90` | `composite_layer` | Beeby inversion for composite layers, interlayer propagators, tensor storage | 5, done |
| `lib/layer_stacking.f90` | `layer_stacking` | layer doubling for the bulk, stacking of the layers above it | 6 |
| `lib/surface_amplitudes.f90` | `surface_amplitudes` | top layer, outgoing amplitudes, intensities, tensor amplitudes | 7 |
| `src/ref-calc.f` | program `LEED` | driver, fixed form; gains `use` statements as stages land | every stage |
| `src/refcalc.f90` | program `refcalc` | free-form driver replacing `src/ref-calc.f` | 8 |

`lib/lib.tleed.f` stays the library the program `delta` links (section 3).

## 3. Routine map

`LEED` calls every routine of `lib/lib.tleed.f` directly or indirectly
except `DEBWAL` and `ZSU`; these two are not ported and leave the library in
stage 8. `delta` (`src/delta.f`, `lib/lib.delta.f`) uses `TSCATF`,
`PSTEMP`, `BESSEL` and `CPPP` and no other routine of the library. These
four stay textually unchanged in `lib/lib.tleed.f` for `delta` (guarded by
`tests/tier1/test_delta.py`); `ref-calc` gets its own module versions. Every
other routine leaves `lib/lib.tleed.f` when it is ported. `lib/lib.tleed.f`
never `use`s a module.

Column "I/O" marks routines that read or write files; their I/O moves to
`leed_input` or `leed_output`.

| Legacy routine | File | Target module | Stage | Notes |
|---|---|---|---|---|
| `LEED` (program) | `ref-calc.f` | driver | every stage; free form in 8 | holds all COMMON blocks; opens `fd.out` (unit 7), `amp.out` (unit 8) |
| `READIN` | `lib.tleed.f` | `leed_input` | 2 | I/O: stdin, echo to stdout; `SLIND` is called after the reader |
| `READGEO` | `lib.tleed.f` | `leed_input` | 2 | I/O: stdin, echo to stdout |
| `SLIND` | `lib.tleed.f` | `bravais_layer` | 3, done | `sublattice_equivalences`; takes the `lattice` type instead of `/SL/`; when a rotation through `360/IDEG` degrees maps a sublattice onto none (`NL > 1`, `IDEG` inconsistent with the lattices), it stops with an error (section 6) naming `IDEG`, the sublattice and the rotation; the legacy code set `JJS` to `NL + 1` there, which `FMAT` uses as an index past `FLMS` |
| `FACT` | `lib.tleed.f` | `angular_coefficients` | 3, done | function `scaled_factorial`; single-precision literals kept |
| `CA` | `lib.tleed.f` | `angular_coefficients` | 3, done | function `composite_layer_coefficient`; default-`REAL` intermediates and the literal `12.56637` kept |
| `CAAA` | `lib.tleed.f` | `angular_coefficients` | 3, done | `composite_layer_coefficients` |
| `BLM` | `lib.tleed.f` | `angular_coefficients` | 3, done | function `harmonics_product_integral`; invalid-argument message on stderr, non-zero exit |
| `CELMG` | `lib.tleed.f` | `angular_coefficients` | 3, done | `matrix_coefficients`; `YLMC`, `FAC1`, `FAC2` are local; `FAC` stays a workspace argument, filled by `matrix_coefficients` and again by `legendre_product_integrals` |
| `LXGENT` | `lib.tleed.f` | `angular_coefficients` | 3, done | `lm_permutations` |
| `CPPP` | `lib.tleed.f` | `atomic_scattering` | 3, done | `legendre_product_integrals`; used by `delta`: frozen copy stays |
| `TSCATF` | `lib.tleed.f` | `atomic_scattering` | 4, done | `atomic_tmatrix`; used by `delta`: frozen copy stays; the low-energy stop message goes to stderr; the sum over elements into the atomic t-matrix is `accumulate_atomic_tmatrix` |
| `PSTEMP` | `lib.tleed.f` | `atomic_scattering` | 4, done | `thermal_phase_shifts` (private); used by `delta`: frozen copy stays |
| `BESSEL` | `lib.tleed.f` | `atomic_scattering` | 4, done | `spherical_bessel` (private); used by `delta`: frozen copy stays |
| `DEBWAL` | `lib.tleed.f` | none | removed 8 | unused |
| `BEAMS` | `lib.tleed.f` | `beam_selection` | 4, done | `select_beams`; sets `NT`, the number of beams in use |
| `FMAT` | `lib.tleed.f` | `bravais_layer` | 4, done | `lattice_sums`; takes the `lattice` type and an `energy_point` instead of `/SL/` and the blank COMMON (`AK(2)`) |
| `RSMF_SIMPLE` | `lib.tleed.f` | `bravais_layer` | 4, done | `bravais_layer_matrices`; `TV` and `LMAX` are arguments; `ROP`/`TOP` and `cylm(nt,lmmax)` stay packed with leading dimension `NT` |
| `MFOLD_SIMPLE` | `lib.tleed.f` | `bravais_layer` | 4, done | `matrix_elements` (private) |
| `XM` | `lib.tleed.f` | `bravais_layer` | 4, done | `intralayer_matrix` (private) |
| `SPHRM_MOD` | `lib.tleed.f` | `bravais_layer` | 4, done | `spherical_harmonics`; the error message goes to stderr |
| `BEAMS_YLM` | `lib.tleed.f` | `bravais_layer` | 4, done | `beams_ylm`, on top of `spherical_harmonics`, with an `energy_point` instead of `E,AK2,AK3,VPI`; called by `bravais_layer_matrices` before it forms the prefactors `AMULT` from the returned `XA`, as `RSMF_SIMPLE` does. Also called by `composite_layer_matrices` and by the driver before `FINAL1` |
| `COPYPOS` | `lib.tleed.f` | `composite_layer` | 5, done | `copy_subplane_positions` |
| `RTINV_SIMPLE` | `lib.tleed.f` | `composite_layer` | 5, done | `composite_layer_matrices`; takes the `lattice` type and an `energy_point`; `/MPT/` and `/MFB/` are arguments (section 4.1); writes the tensor amplitudes to `tstore`, packed as `TSTORE(2,LMN,NT)` and indexed by the beam within its set, without the beam-set offset (with several beam sets of a composite bulk-type layer each call overwrites the columns of the previous one; kept: no output can show it, because the input reader rejects tensor output for a composite layer of the overlayer type (`LAY = 2`) with a non-trivial superlattice, `leed_input.f90`, so only the routine test sees the offset). The `SRTLAY`, `GHMAT` and `TAUMAT` calls stay inside it; the branches that `NEW.EQ.1` (`RTINV_SIMPLE`) and `NEW+1` (`SRTLAY`) never take are gone |
| `MFOLT_SIMPLE` | `lib.tleed.f` | `composite_layer` | 5, done | `reflection_transmission_element` (private); `GP` and `LM` (`lmax + 1`) are arguments |
| `TAUMAT` | `lib.tleed.f` | `composite_layer` | 5, done | `subplane_scattering`; takes the `lattice` type and an `energy_point` |
| `TAUY_SIMPLE` | `lib.tleed.f` | `composite_layer` | 5, done | `tau_ylm` (private) |
| `THMAT` | `lib.tleed.f` | `composite_layer` | 5, done | `beeby_matrix` (private) |
| `GHMAT` | `lib.tleed.f` | `composite_layer` | 5, done | `interplanar_propagators`, with the private helpers `add_reciprocal_terms` and `fill_symmetric_elements` (the symmetric fill-in shared with `propagator_from_sum`) and `decode_lm` (the (l,m) decoding); one work buffer with both shapes explicit (section 4.2); `NUGH` is modified in place, as in the legacy routine; `L2M = 2*MLMAX + 1` sizes the direct-sum loops, so extra terms are computed and unused (kept; no effect on any output) |
| `GHD` | `lib.tleed.f` | `composite_layer` | 5, done | `direct_propagator_sum` (private) |
| `GHSC` | `lib.tleed.f` | `composite_layer` | 5, done | `propagator_from_sum` (private) |
| `SB` | `lib.tleed.f` | `composite_layer` | 5, done | `spherical_hankel` (private; computes the Hankel functions of the first kind, `atomic_scattering` has its own `spherical_bessel`) |
| `SH` | `lib.tleed.f` | `composite_layer` | 5, done | `harmonics_table` (private); `REAL*8`/`COMPLEX*16` locals: `dp` |
| `SRTLAY` | `lib.tleed.f` | `composite_layer` | 5, done | `sort_subplanes`; called by `composite_layer_matrices` before `interplanar_propagators` |
| `CXMTXT` | `lib.tleed.f` | `composite_layer` | 5, done | `solve_pivoted` (private); the real-modulus pivot search and its elimination threshold `1.0E-10` are kept; only the solve path (`INOPT = -1`) that `TAUMAT` uses is ported: the inverse and determinant options and the singularity threshold, which only set the determinant mark, are dropped |
| `XMT` | `lib.tleed.f` | `composite_layer` | 5, done | `multiple_scattering_matrix` (private) |
| `DBLGAS_MOD` | `lib.tleed.f` | `layer_stacking` | 6 | |
| `SUBRAS_BLAS` | `lib.tleed.f` | `layer_stacking` | 6 | messages on stdout |
| `TLRTA_BLAS` | `lib.tleed.f` | `layer_stacking` | 6 | |
| `DBGT_MOD` | `lib.tleed.f` | `layer_stacking` | 6 | |
| `DECIDE` | `lib.tleed.f` | `beam_selection` | 7 | sets the k-parallel of blank COMMON; message on stdout |
| `GET_TOPLAY_PROPAGATORS` | `lib.tleed.f` | `surface_amplitudes` | 7 | |
| `GET_TOPLAY_MULTSCATT` | `lib.tleed.f` | `surface_amplitudes` | 7 | called with 7 arguments, declared with 6 |
| `ADREF2T` | `lib.tleed.f` | `surface_amplitudes` | 7 | |
| `MULTAMP_OPT` | `lib.tleed.f` | `surface_amplitudes` | 7 | |
| `RINT_SIMPLE` | `lib.tleed.f` | `surface_amplitudes`, output to `leed_output` | 7 | I/O: writes `fd.out`, `amp.out` |
| `FINAL1` | `lib.tleed.f` | `surface_amplitudes` | 7 | |
| `ZGE` | `lib.tleed.f` | `surface_amplitudes` | 7 | |
| `ZTU` | `lib.tleed.f` | `surface_amplitudes` | 7 | |
| `ZSU` | `lib.tleed.f` | none | removed 8 | unused |
| `FINALOV` | `lib.tleed.f` | `surface_amplitudes` | 7 | |
| `HEAD` | `lib.tleed.f` | `leed_output` | 7 | I/O: headers of `fd.out`, `amp.out` |
| `OPENOUT` | `lib.tleed.f` | driver (opening), `leed_output` | 7 | I/O: opens the tensor files |
| `OUTAMP` | `lib.tleed.f` | `leed_output` | 7 | I/O: tensor files |
| `OUTXIST` | `lib.tleed.f` | `leed_output` | 7 | I/O: tensor files |
| `DELIMIT` | `lib.tleed.f` | `leed_output` | 7 | I/O: tensor files |

## 4. COMMON blocks and storage association

No COMMON block survives the port: a ported routine receives what it took
from COMMON as arguments or derived-type components. As long as a legacy
caller remains, it passes its COMMON members as actual arguments; a block
disappears with its last legacy user, and stage 8 removes the driver's
declarations. No COMMON block is initialised by `DATA`, and none is
overlaid by `EQUIVALENCE`.

### 4.1 COMMON blocks

| Block | Declared in (members) | Becomes |
|---|---|---|
| blank COMMON | `LEED` (`E,AK21,AK31,VPI,TV,EMACH`); `GET_TOPLAY_PROPAGATORS` (`E,AK2,AK3`); `DBGT_MOD`, `DBLGAS_MOD`, `GHD`, `GHMAT`, `RTINV_SIMPLE`, `SUBRAS_BLAS`, `TAUMAT` (`E,AK2,AK3,VPI`); `RSMF_SIMPLE` (`E,AK2,AK3,VPI,TV`); `RINT_SIMPLE` (`E,CK2,CK3,VPI`); `FMAT` (`E,AK,VPI` with `AK(2)`) | `leed_types`: energy point (energy, k-parallel of the incident beam for which the layer matrices are computed, imaginary inner potential); `TV` a cell-area argument of `RSMF_SIMPLE`; `EMACH` the named constant `MEMACH` |
| `/SL/` | `LEED`, `READIN`, `GHD`, `GHMAT` (`ARA1,ARA2,ARB1,ARB2,RBR1,RBR2,NL1,NL2`); `FMAT`, `SLIND`, `TAUMAT` (`BR1,BR2,AR1,AR2,RAR1,RAR2,NL1,NL2`) | `leed_types`: lattice (substrate and superlattice vectors, superlattice reciprocal vectors, `NL1`, `NL2`) |
| `/MS/` | `LEED`, `READIN`, `RSMF_SIMPLE`, `SUBRAS_BLAS` (`LMAX,EPS,LITER`); `RTINV_SIMPLE` (`LMAX,EPS1,LITER`) | `leed_types`: control parameters |
| `/ADS/` | `LEED`, `GET_TOPLAY_PROPAGATORS`, `OUTXIST`, `READIN` (`ASE,VPIS,VPIO,VV`) | `leed_types`: `ASE`, `VPIS`, `VPIO` with the inner potential; `VV` (real inner potential at the current energy) with the energy point |
| `/INPOT/` | `LEED`, `READIN` (`EM,C1,...,C8,WORKFN`) | `leed_types`: inner-potential model, components default-initialised to zero |
| `/MPT/` | `LEED` (`NA,NS,ID,LAY,L1,NTAU,TSTS,TV1,DCUTS,NOPT,NEW`); `RTINV_SIMPLE` (`NA,NS,ID,LAY,LM,NTAU,TST,TV,DCUT,NOPT,NEW`) | gone: its members are arguments of `composite_layer_matrices`; `ID`, `NOPT` and `NEW` dropped |
| `/MFB/` | `RTINV_SIMPLE` (`GP,L1`); `MFOLT_SIMPLE` (`GP,LM`) | gone: `GP(2)` and `LM` (`lmax + 1`) are arguments of `reflection_transmission_element` |
| `/BT/` | `LEED` only (`IT,T,T0,DRPER,DRPAR,DR0`) | nothing |

Hazards of each block:

- **blank COMMON.** The block has different lengths in different routines
  (legal only for blank COMMON). `FMAT` declares `AK(2)` and so reads the
  two k-parallel scalars `AK2`, `AK3` of the other routines as one array:
  the port passes k-parallel as a two-component array to every routine.
  `LEED` names the COMMON k-parallel `AK21`, `AK31` and has *local* `AK2`,
  `AK3` for the true incident beam; `DECIDE` sets `AK21`, `AK31` to the
  k-parallel for which the layer matrices are computed next (the incident
  beam, or the reversed exit beam `NEXIT` of the tensor loop), so `AK2` in
  the library is not `AK2` in the driver. `RINT_SIMPLE` calls the same storage `CK2`,
  `CK3` and has a local `AK`. `TV` in the blank COMMON (set by `LEED` to
  the overlayer or substrate cell area before each `RSMF_SIMPLE` call) is
  another variable than `TV` of `/MPT/` in `RTINV_SIMPLE`. `EMACH` is set
  by `LEED` and read by no routine through COMMON.
- **`/SL/`.** `FMAT`, `SLIND` and `TAUMAT` rename the members positionally:
  their `BR1`, `BR2` are the substrate vectors `ARA1`, `ARA2`; their `AR1`,
  `AR2` the superlattice vectors `ARB1`, `ARB2`; their `RAR1`, `RAR2` the
  superlattice reciprocal vectors `RBR1`, `RBR2`. `READIN` also has local
  `RAR1`, `RAR2` that hold the *substrate* reciprocal vectors. The derived
  type uses one name per quantity everywhere.
- **`/MS/`.** `RTINV_SIMPLE` names `EPS` `EPS1`; `RTINV_SIMPLE` and
  `RSMF_SIMPLE` use only `LMAX`.
- **`/ADS/`.** `VPIS` and `VPIO` are both set to `VPI` by `READIN` and
  written to the tensor files by `OUTXIST`; `GET_TOPLAY_PROPAGATORS` reads
  `ASE` and `VPIO`.
- **`/INPOT/`.** `READIN` reads `EM`, `C1`...`C8` only for phase-shift
  format 2; for format 1 it reads `C1`...`C4` and leaves `EM`, `C5`...`C8`
  undefined. `LEED` branches on the sign of `EM` (`IF (EM) 140,120,110`) to
  choose the inner-potential formula, so format 1 relies on COMMON storage
  being zero at start-up, which the standard does not guarantee. The
  derived type initialises these components to zero explicitly, which keeps
  the observed behaviour.
- **`/MPT/`.** Used only to pass values from `LEED` to `RTINV_SIMPLE`, with
  renames: `L1` becomes `LM`, `TSTS` `TST`, `TV1` `TV`, `DCUTS` `DCUT`.
  `LEED` also passes `TSTS` as an actual argument to `READIN`, which sets
  it. `LEED` sets `ID`, `NOPT` and `NEW` to 1; `RTINV_SIMPLE` does not use
  `ID` or `NOPT`. `NEW` is read only by `RTINV_SIMPLE` (tests `NEW.EQ.1`)
  and `SRTLAY` (arithmetic IF on `NEW+1`), whose other branches are dead
  because `NEW` is always 1; `GHMAT` does not read it. The port keeps the
  live branches and no `NEW`.
- **`/MFB/`.** Hidden state between two library routines: `RTINV_SIMPLE`
  sets `GP` (the k-parallel offset of incident beam `JGP`) inside its beam
  loop just before calling `MFOLT_SIMPLE`, which reads it through COMMON;
  nothing passes it as an argument. `RTINV_SIMPLE` sets the second member
  (its `L1`, copied from `/MPT/`'s `LM`, equal to `LMAX + 1`);
  `MFOLT_SIMPLE` calls it `LM` and uses it as the number of angular
  momenta in its propagation phases. Note that `L1` of `/MFB/` in `RTINV_SIMPLE` and `L1` of
  `/MPT/` in `LEED` are different storage with the same value.
- **`/BT/`.** Dead: declared only in `LEED`. Its members are set in `LEED`
  and passed to `TSCATF` as ordinary arguments; no library routine declares
  the block (a comment in `RSMF_SIMPLE` still says it is passed to `MFOLD`).

### 4.2 Storage association the port must keep or undo deliberately

- **Matrices packed with leading dimension `NT` in `MKNT`-sized arrays.**
  `LEED` declares the beam-space matrices with the maximum beam count, e.g.
  `ROP(MKNT,MKNT)`, `RAMP(MKNT,MKNT)`, `S1(MKNT,MKNT)`, `CYLM(MKNT,MLMMAX)`,
  `RG(2,MNSUB,MKNT)`, and passes them to dummies declared with `NT`, the
  number of beams in use at the current energy (`BEAMS`), e.g.
  `RA(NT,NT)` in `RSMF_SIMPLE`, `R_TOP(NT,NT)` and `CYLM(NT,LMMAX)` in
  `RTINV_SIMPLE`, `RAMP(N,N)` in `DBLGAS_MOD` and `DBGT_MOD`. For
  `NT < MKNT` the callee's matrix is packed into the first `NT*NT`
  elements, and the driver's element `(I,J)` is not the callee's. The
  driver copies these matrices to and from the per-layer stores
  `REPSTO(MNLTYPE,MKNT,MKNT)` etc. over the full `MKNT` range, which moves
  the packed data intact. The port either keeps one consistent view (arrays
  allocated `NT` by `NT` per energy, the per-energy sizing of stage 8) or
  keeps packed buffers; mixing the two views breaks results only for
  energies with `NT < MKNT`. The same holds for the arrays `LEED` sizes
  with `PARAM` maxima and passes to `RTINV_SIMPLE`, which declares them with
  the dimensions of the current layer: `TAU(MLMT,MLEV)` as `TAU(LMT,LEV)`,
  `XEV(MLEV,MLEV2)` as `XEV(LEV,LEV2)`, `GH(MLMG,MLMMAX)` as
  `GH(LMG,LMMAX)`, `TH(MLMNI,MLMNI)` as `TH(LMNI,LMNI)` (`th(lmn,lmn)` in
  the port) and the vectors `TS`, `TAUG`, `TAUGM`, `XH`, `RG1`, `RG2` and
  `HGHD`. Likewise the arrays `LEED` sizes with the sublayer maxima `MNSUB`
  and `MNSUB2` and passes to `RTINV_SIMPLE` and on to `SRTLAY`, `GHMAT` and
  `COPYPOS`, which declare them with the sublayer count `NLAY` of the
  current layer and `NLAY2 = NLAY*(NLAY-1)/2`: `POS(MNSUB,3)` and
  `POSS(MNSUB,3)` as `POS(NLAY,3)`, `POSS(NLAY,3)`, `MGH(MNSUB,MNSUB)` as
  `MGH(NLAY,NLAY)` and `DRL(MNSUB2,3)` as `DRL(NLAY2,3)`, whose elements
  move for `NLAY < MNSUB`, and the vectors `LPS`, `LPSS`, `NORD`
  (`MNSUB`), `TEST`, `NUGH`, `NGEQ` (`MNSUB2`) and `IPL` (`MLMNI` as
  `LMN`), which keep their element order. The port keeps these
  as packed leading parts of the larger arrays (sized by stage 8 from the
  derived values, they are packed exactly).
- **Element-as-array passing.**
  - `TENSOV(2,MLMN,MKNT,MNCOMP)`: `LEED` passes `TENSOV(1,1,1,LCOMP(...))`
    to `RTINV_SIMPLE` and `FINALOV`, which declare `TSTORE(2,LMN,NT)`: one
    composite layer's slice, re-shaped and packed with `LMN` and `NT`. The
    port passes the explicit section `tensov(:,:,:,k)`, which is contiguous,
    to a dummy `tstore(2,lmn,nt)` that takes the packed leading part.
  - `EFSCAM`, `EFSCAP(MKNT,MKNT,MNSTACK)`: `LEED` passes
    `EFSCAM(1,1,ISTACK)` and `EFSCAP(1,1,ISTACK)` to `DBGT_MOD` (`S1(N,N)`,
    `S2(N,N)`) and to `MULTAMP_OPT` (`MATRIX(N,N)`): one stack slice,
    packed with leading dimension `NT`.
  - `PQ(2,NT)`: `RTINV_SIMPLE` passes `PQ(1,1+NA)` to `GHMAT`, which
    declares `PQ(2)`: the k-parallel offset of the first beam of the
    current beam set.
  The port passes explicit sections (`tensov(:,:,:,k)`,
  `efscam(:,:,istack)`; `composite_layer_matrices` takes `pq` whole and
  hands `pq(:,1+na)` on) to dummies of matching shape, or
  keeps contiguous packed buffers on purpose; never element-as-array.
- **One actual for two dummies.** `RTINV_SIMPLE` passes its work array
  `VT(LM2N)` to `GHMAT` twice, as `Y1(L2M,L2M)` (handed to `GHD`) and as
  `Y(LM,LM)` (handed to `SH`): one buffer seen with two shapes. The port
  gives `interplanar_propagators` one work buffer `work(l2m,l2m)`, which it
  hands whole to `direct_propagator_sum` (`y(l2m,l2m)`) and to the table
  routine and the reciprocal-sum routine (`y(lm,lm)`, the leading `lm*lm`
  elements). The rank-1 `VT(LM2N)` of `composite_layer_matrices` is passed
  to `work(l2m,l2m)` as one buffer, in the same way. `TG(2,LM2N)`, the
  work vector `S(LMS)` of `GHMAT`, is declared as a vector `tg(2*lm2n)`, so
  that only the length differs.
- **Argument-count mismatch.** `LEED` calls `GET_TOPLAY_MULTSCATT` with 7
  arguments (`..., NT, NEXIT`); the routine declares 6. The port's explicit
  interface removes the extra argument; behaviour is unchanged since the
  routine never sees it.

## 5. Sizes and the input-only values

Every `PARAM` value is classified in `tests/harness/param_size.py`: strict
(`READIN`/`READGEO` stop unless it equals the count in the input), bound (at
least the count in the input) or explicit (reaches `ref-calc` only through
`PARAM`). The tier-1 tests check that each strict value equals and each
bound value covers what the input implies, and that a minimal `PARAM`
derived from the input gives identical outputs. The port sizes every array
from these derived values (stage 8), never from a maximum; `NCA(LMAX)` and
`NLMS(LMAX)` (the `DATA` tables in `LEED`) give the sizes of `CAA` and
`CLM`.

The four explicit values cannot be derived from the input and change
results or parsing. `LEED` copies them from `PARAM` at start-up
(`IDEG = MIDEG`, `NL1 = MNL1`, `NL2 = MNL2`, `NPUN = MNPUN`,
`NT0 = MNT0`, `NPSI = MNPSI`):

| Value | Enters |
|---|---|
| `IDEG` (`MIDEG`) | `READIN` to `SLIND` (`JJS(NL,IDEG)`) and `FMAT` (sector sum over `IDEG` sublattice orientations, `SCC`, `SA`) |
| `NL1`, `NL2` (`MNL1`, `MNL2`) | `/SL/`; `SLIND`'s coset grid; `READIN` reads `NL1*NL2` beam sets |
| `NPUN` (`MNPUN`) | `READIN` reads `NPU(1:NPUN)`; `HEAD` (headers of `fd.out`, `amp.out`) and `BEAMS`; `NT0 = MNT0`, equal to it, is the number of exit beams of the tensor loop in `LEED` and the length of `XIST` written by `OUTXIST` |
| `NPSI` (`MNPSI`) | `READIN`'s phase-shift energy loops, and so where `READGEO` starts reading |

For TensErLEED >= 3.0.0 the input carries these values in three records,
each next to the data it sizes; the reader requires all three, and text
after the numbers is a comment, as elsewhere in the input. The reader
reads them list-directed, so any field width works (blank- or
comma-separated values, and no `/` before the last number); the format
below is what the converter and viperleed write:

| Record | Written as | Place |
|---|---|---|
| `IDEG`, `NL1`, `NL2` | `3I5` | directly after the `VPI`/`ASE`/`WORKFN` record, before the beam sets |
| `NPUN` | `I5` | directly before the `NPU` record(s) |
| `NPSI` | `I5` | directly after the `IFORM` record, before the phase-shift block |

Line 2 does not move, and no comment of these records contains the text
`   LMAX`: viperleed edits the energy range and `LMAX` of an existing input
by line position and by that text. Example (Ir(100)-(2x1)-O, comments
starting in column 46 as for `IFORM`):

```
     5.00     1.00     2.65   V0i,ASE,WORKFN
    2    2    1                              IDEG NL1 NL2 (TL >= 3.0.0)
...
   12                                        NPUN (TL >= 3.0.0)
    2    4    6   10   12   14   18   46   48   52   54   58
...
  1                                          IFORM - formatted input and output
   44                                        NPSI (TL >= 3.0.0)
```

The ported reader reads only this format. As long as `PARAM` exists, it
stops when a value differs from its `PARAM` counterpart. viperleed writes
these records from stage 9 on; until then the test suite converts each
frozen input for the tree under test (`tests/harness/fin_convert.py`),
while the legacy program reads the original. The output formats do not
change.

## 6. Conventions for the new files

- **Standard and headers.** New files are free-form Fortran 2018 and
  compile with `-std=f2018` (`STD_FFLAGS` in the Makefile;
  `tests/unit/test_makefile.py` checks every free-form file the Makefile lists); the
  fixed-form `src/ref-calc.f` and `lib/lib.tleed.f` are exempt. Each file
  starts with a Doxygen `\file`/`\brief` header as in `lib/kinds.f90`;
  history and attribution comments of the ported legacy routines move with
  them. Headers do not claim a language version ("Fortran 90 version").
- **Kinds.** Every real and complex uses `sp` or `dp` from `kinds`; legacy
  default `REAL`/`COMPLEX` become `sp`, `REAL*8`/`COMPLEX*16` become `dp`.
  `lib/kinds.f90` is identical to `module kinds` in `lib/rfacsb.f90`
  (`tests/unit/test_kinds.py`) until `rfactor` takes it from there too.
- **Interfaces.** Every procedure is a module procedure with
  `implicit none` and an `intent` on every dummy; arrays are assumed-shape
  or explicit-shape matching the actual argument (section 4.2). No COMMON,
  no `save`-carried state, no module variables as hidden globals.
- **Names.** Descriptive names for routines and variables, the legacy name
  in a comment where it helps to compare with the legacy code and the
  TensErLEED paper (Blum and Heinz, Comput. Phys. Commun. 134, 392 (2001)).
- **I/O.** Only `leed_input` reads and only `leed_output` writes. The
  driver opens and closes every file; reader and writer routines take the
  unit number as an argument (stdin and stdout included), so no unit number
  is shared implicitly. A computational routine may call a writer routine
  at the point where the legacy code wrote, to keep the legacy order of the
  output; the writer takes values only, and the computational routine
  opens, reads and positions no file.
- **Memory.** Arrays that scale with the number of beams squared or larger
  are allocatable, never automatic; workspace buffers are allocated once
  and passed explicitly; no array temporaries for large arrays (checked
  with `-Warray-temporaries`).
- **Errors.** Invalid input or inconsistent sizes stop with a message on
  stderr and a non-zero exit code.
- **File names.** viperleed (before TensErLEED 3.0.0) builds `ref-calc`
  from the first matches of `lib/lib.tleed*` and `src/ref-calc*`, and
  `delta` from `lib/lib.tleed*`, `lib/lib.delta*` and `src/delta*`; it builds
  `rfactor` from `lib/rfacsb*` and `src/rfactor*`. New files therefore never
  start with `lib.tleed`, `ref-calc`, `delta`, `rfacsb` or `rfactor`.

## 7. Open points

Choices made here:

1. `SLIND` in stage 3 (energy-independent set-up) although `READIN`
   (stage 2) called it; the driver calls the port after the ported reader
   with the reader's `lattice` and still copies the lattice into `/SL/`
   for `FMAT`, `GHD`, `GHMAT` and `TAUMAT` (`SLIND` writes nothing on stdout,
   so the stdout order is unchanged).
2. `BEAMS` in stage 4 (first stage running the per-energy computation),
   as `select_beams` in `beam_selection`; `DECIDE` in stage 7 with the
   exit-beam loop of the tensor output. Settled for stage 4: `BEAMS` does
   not call `DECIDE`, and `select_beams` is self-contained. Open for stage
   7: `DECIDE` joins `beam_selection` and then sets the k-parallel of an
   `energy_point` instead of blank COMMON.
3. `CPPP` in stage 3 (called once before the energy loop), `TSCATF`,
   `PSTEMP`, `BESSEL` in stage 4.
4. `SPHRM_MOD` in `bravais_layer` (stage 4), as `spherical_harmonics`. The
   legacy copy stayed in `lib.tleed.f` for `RTINV_SIMPLE` until stage 5,
   which replaced its use by `beams_ylm` in `bravais_layer`, built on
   `spherical_harmonics`; the legacy `SPHRM_MOD` and `BEAMS_YLM` are gone.
5. `leed_output` created in stage 2 for the input echo on stdout of
   `READIN`/`READGEO`; the output files in stage 7.
6. The driver stays fixed-form `src/ref-calc.f` until stage 8, which
   replaces it with `src/refcalc.f90`.
7. The input-only values of section 5 change the input format (records
   and places in section 5). Tier-1 inputs stay frozen: the test suite
   gives a tree with the ported reader a converted copy of each input
   (values from the run's `PARAM`) and the legacy program the original;
   which one a tree gets is decided from the tree, not from the input.
