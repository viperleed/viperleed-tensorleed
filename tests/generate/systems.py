"""Definitions of the test systems used by the generators (tiers 1 and 2).

Each system is a viperleed example input set plus PARAMETERS overrides.
Run inside an environment with viperleed installed (generators only; pytest never
imports this module).
"""
from dataclasses import dataclass, field
from pathlib import Path

from harness.paths import location

EXAMPLES = location("viperleed_examples")


@dataclass
class SystemDef:
    """A test system: its inputs, the PARAMETERS changes, and what its runs must cover."""
    name: str
    source_zip: Path | None  # None: inputs from base_files
    # PARAMETERS overrides: key (left of '=', incl. flags) -> value.
    # Existing lines with the same key are commented out. A value of None
    # only comments out the key.
    parameters: dict = field(default_factory=dict)
    # Extra input files written into the calc directory: name -> content
    extra_files: dict = field(default_factory=dict)
    # Code paths that the frozen FIN must exhibit (checked by
    # make_refcalc_cases.check_coverage): see COVERAGE_CHECKS there.
    must_cover: list = field(default_factory=list)
    covers: list = field(default_factory=list)  # informative, -> case.toml
    description: str = ""
    # Runs; None = select automatically from viperleed's energy -> LMAX map:
    # lmax_runs: [(energy, lmax)] single-energy runs,
    # multi_energy_run: ([energies], lmax)
    lmax_runs: list | None = None
    multi_energy_run: tuple | None = None
    base_files: dict = field(default_factory=dict)  # name -> content, if no source_zip


_AG = EXAMPLES / "Ag(100)" / "input_files.zip"
_IR = EXAMPLES / "Ir(100)-O" / "input_files.zip"

_AG_BASE = {"THEO_ENERGIES": "45 210 3", "LMAX": "6-9"}
_AG_COVERS = ["bravais", "energy-dependent-lmax", "tensors", "layer-doubling"]

# Cu(111)-(1x1), defined here (no viperleed example): a = 2.556 A, ideal
# interlayer distance a*sqrt(2/3), ABC stacking in a 120-degree cell. The bulk
# repeat vector has a lateral shift whose mirror image is the other (hcp)
# hollow, so in-plane sign errors are not equivalent as on fcc(100).
_CU111_POSCAR = """Cu(111) 1x1 clean
  1.0
    2.5560000000000000    0.0000000000000000    0.0000000000000000
   -1.2780000000000000    2.2135609320730250    0.0000000000000000
    0.0000000000000000    0.0000000000000000   20.0000000000000000
   Cu
    6
Direct
  0.0000000000000000  0.0000000000000000  0.8000000000000000
  0.3333333333333333  0.6666666666666666  0.6956517369574367
  0.6666666666666666  0.3333333333333333  0.5913034739148733
  0.0000000000000000  0.0000000000000000  0.4869552108723099
  0.3333333333333333  0.6666666666666666  0.3826069478297465
  0.6666666666666666  0.3333333333333333  0.2782586847871831
"""
_CU111_PARAMETERS = """RUN = 0
BULK_LIKE_BELOW = 0.43
SITE_DEF Cu = surf 1
T_DEBYE = 343
T_EXPERIMENT = 100
VIBR_AMP_SCALE = *surf 1.3
"""
_CU111_IVBEAMS = """  h         k
  1.000000   0.000000
  0.000000   1.000000
 -1.000000   0.000000
  1.000000   1.000000
  2.000000   0.000000
"""

SYSTEMS = {
    "ag100": SystemDef(
        "ag100", _AG, dict(_AG_BASE),
        must_cover=["nsub_all_1", "tensor_layers=5"],
        covers=_AG_COVERS,
        description="Ag(100), normal incidence, Bravais layers, one bulk layer"),
    "ag100_oblique": SystemDef(
        "ag100_oblique", _AG,
        {**_AG_BASE, "BEAM_INCIDENCE": "THETA 12, PHI 33"},
        must_cover=["theta_nonzero"],
        covers=_AG_COVERS + ["oblique-incidence", "decide-kpar"],
        description="Ag(100), oblique incidence (theta 12, phi 33)"),
    "ag100_bulk2": SystemDef(
        "ag100_bulk2", _AG,
        {**_AG_BASE,
         # bulk detection would reset N_BULK_LAYERS to 1
         "BULK_LIKE_BELOW": None,
         "BULK_REPEAT": "[0 0 4.07294]",   # 2 x the single-layer repeat z
         "LAYER_CUTS": "0.35 0.45 0.55 0.65 0.75",
         "N_BULK_LAYERS": "2",
         # the lower bulk layer (atom 6) gets its own site with another
         # vibration amplitude, so the two bulk layers differ numerically
         # (with identical layers, an A<->B mix-up in DBLGAS passes)
         "SITE_DEF Ag": "surf 1, bulkb 6"},
        extra_files={"VIBROCC": (
            "= Vibration Amplitudes\n"
            "Ag_surf = 0.09\n"
            "Ag_def = 0.069\n"
            "Ag_bulkb = 0.10\n")},
        must_cover=["bulk_two_layers", "bulk_layers_distinct_sites"],
        covers=_AG_COVERS + ["two-layer-bulk"],
        description="Ag(100), bulk repeat unit of two distinct layers (different sites)"),
    "ag100_mix": SystemDef(
        "ag100_mix", _AG,
        {**_AG_BASE, "ELEMENT_MIX Ag": "Ag Pd"},
        extra_files={"VIBROCC": (
            "= Vibration Amplitudes\n"
            "Ag_surf = Ag 0.09, Pd 0.09\n"
            "Ag_def = Ag 0.069, Pd 0.069\n"
            "= Occupations\n"
            "Ag_surf = Ag 0.7, Pd 0.3\n"
            "Ag_def = Ag 1.0, Pd 0.0\n")},
        must_cover=["two_elements", "mixed_site"],
        covers=_AG_COVERS + ["element-mix", "ata-tmatrix"],
        description="Ag(100), surface site occupied by 0.7 Ag / 0.3 Pd"),
    "ir100_2x1_o": SystemDef(
        "ir100_2x1_o", _IR,
        {"THEO_ENERGIES": "50 170 4", "LMAX": "6-8"},
        must_cover=["nsub_gt_1", "beam_sets_gt_1", "two_elements"],
        covers=["composite-layer", "beeby", "superstructure", "several-beam-sets",
                "two-elements", "mideg2", "energy-dependent-lmax", "tensors"],
        description="Ir(100)-(2x1)-O, composite layers, several beam sets"),
    "ir100_2x1_o_thick": SystemDef(
        "ir100_2x1_o_thick", _IR,
        # O alone, then a layer of two Ir planes (4 sublayers): the vectors
        # between the Ir atoms of the two planes are pairwise equal (NGEQ in
        # SRTLAY/GHMAT). (O + two Ir planes in one layer amplified FMA noise
        # in the ALM ~7x beyond the T_* tolerance.)
        # Bulk: one layer of the two lowest Ir planes (repeat 2 x 1.91739 A
        # without lateral shift), i.e. a composite layer with bulk periodicity
        # (LAY = 2, NSUB = 2) under a (2x1) overlayer (a path no other case reaches).
        {"THEO_ENERGIES": "50 170 4", "LMAX": "6-8",
         "BULK_LIKE_BELOW": None,
         "BULK_REPEAT": "[0 0 3.83478]",
         "LAYER_CUTS": "0.25 0.35 0.45 0.55 0.73",
         "N_BULK_LAYERS": "1"},
        must_cover=["nsub_ge_4", "bulk_composite", "beam_sets_gt_1", "two_elements"],
        covers=["composite-layer", "beeby", "equal-interplanar-vectors",
                "composite-bulk-layer", "rtinv-substrate", "superstructure",
                "several-beam-sets", "two-elements", "energy-dependent-lmax", "tensors"],
        description="Ir(100)-(2x1)-O, layer of two Ir planes, composite bulk layer"),
    "cu111": SystemDef(
        "cu111", None, {"THEO_ENERGIES": "50 140 3", "LMAX": "5-7"},
        base_files={"POSCAR": _CU111_POSCAR, "PARAMETERS": _CU111_PARAMETERS,
                    "IVBEAMS": _CU111_IVBEAMS},
        must_cover=["nsub_all_1", "hexagonal"],
        covers=["bravais", "fcc111", "hexagonal", "lateral-bulk-shift", "tensors",
                "layer-doubling"],
        description="Cu(111)-(1x1), fcc stacking (non-equivalent mirrored bulk shift)"),
}

# Large systems for the opt-in benchmark (tests/benchmark): inputs only, two
# single-energy runs (lowest and highest energy). Choice: see
# tests/benchmark/README.md.
_FE2O3 = EXAMPLES / "Fe2O3(012)-(1x1)" / "input_files.zip"
_CUTE = EXAMPLES / "Cu(111)-Te" / "input_files.zip"
_BENCH_PARAMS = {"THEO_ENERGIES": "100 400 10",  # step sets the phase-shift grid
                 "LMAX": "8-14", "LOG_DEBUG": None}

# Candidate not used: peak RSS 83 MB at 400 eV (LMAX 10) vs 249 MB for
# Cu(111)-Te (see tests/benchmark/README.md).
BENCHMARK_CANDIDATES = {
    "bench_fe2o3": SystemDef(
        "bench_fe2o3", _FE2O3, dict(_BENCH_PARAMS),
        covers=["benchmark"], description="Fe2O3(012)-(1x1), benchmark"),
}

BENCHMARKS = {
    "bench_cu111_te": SystemDef(
        "bench_cu111_te", _CUTE, dict(_BENCH_PARAMS),
        covers=["benchmark"], description="Cu(111)-Te (5 0, 1 2), benchmark"),
}
