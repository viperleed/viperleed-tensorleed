"""The constants of lib/constants.f90 equal those of the legacy code, bit for bit.

BOHR, HARTREE and MEMACH are PARAMETERs of src/GLOBAL (default REAL); pi is
the literal 3.14159265 that the legacy routines write; composite_layer's
estimate of reciprocal-lattice points uses the literal 3.1415926. Both sides are
compiled by gfortran, so that the comparison is of the single-precision bits
the compiler produces, not of decimal texts.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2026 ViPErLEED developers'
__created__ = '2026-11-11'
__license__ = 'GPLv3+'

import re
import subprocess

from harness.build import program_sources

_NAMES = ("bohr", "hartree", "memach")
_GLOBAL_PARAMETER = re.compile(r"^\s*PARAMETER\s*\(\s*(\w+)\s*=\s*([^)\s]+)\s*\)",
                               re.I | re.M)
_LEGACY_PI = "3.14159265"

# the literal of composite_layer's estimate of the number of reciprocal-lattice
# points, which differs from pi in the last digit (bits as gfortran gives them)
_PI_ESTIMATE_BITS = "40490FDA"

_PI_ESTIMATE_PROGRAM = """\
program print_pi_estimate
  use kinds, only: sp
  use composite_layer, only: pi_estimate
  implicit none
  write (*, '(z8.8)') transfer(pi_estimate, 1_4)
end program print_pi_estimate
"""

_PROGRAM = """\
program print_bits
  use kinds, only: sp
  use constants, only: bohr, hartree, memach, pi
  implicit none
  real(sp), parameter :: legacy(4) = [{legacy}]
  real(sp), parameter :: module_values(4) = [bohr, hartree, memach, pi]
  integer :: i
  do i = 1, 4
    write (*, '(2z9.8)') transfer(module_values(i), 1_4), transfer(legacy(i), 1_4)
  end do
end program print_bits
"""


def _legacy_values(tl_root):
    """The GLOBAL PARAMETERs bohr, hartree, memach and pi, as Fortran literals."""
    text = (tl_root / "TensErLEED" / "src" / "GLOBAL").read_text()
    found = {name.lower(): value for name, value in _GLOBAL_PARAMETER.findall(text)}
    assert set(_NAMES) <= set(found), f"GLOBAL defines {sorted(found)}"
    return [found[n] + "_sp" for n in _NAMES]


def test_constants_equal_legacy_values(tl_root, tmp_path):
    """bohr, hartree, memach (GLOBAL) and pi (legacy literal) have the same bits
    in module constants."""
    lib = tl_root / "TensErLEED" / "lib"
    legacy = [*_legacy_values(tl_root), _LEGACY_PI + "_sp"]
    source = tmp_path / "print_bits.f90"
    source.write_text(_PROGRAM.format(legacy=", ".join(legacy)))
    exe = tmp_path / "print_bits"
    proc = subprocess.run(["gfortran", "-std=f2018", "-J", str(tmp_path), "-o", str(exe),
                           str(lib / "kinds.f90"), str(lib / "constants.f90"),
                           str(source)], capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr[-3000:]
    out = subprocess.run([str(exe)], capture_output=True, text=True, check=True)
    rows = [line.split() for line in out.stdout.splitlines()]
    diffs = [f"{name}: module {got} != legacy {want}"
             for name, (got, want) in zip((*_NAMES, "pi"), rows) if got != want]
    assert not diffs, "\n".join(diffs)


def test_composite_layer_pi_estimate_bits(tl_root, tmp_path):
    """The pi of the estimate of reciprocal-lattice points in composite_layer
    (a private parameter, exposed in a copy of the module) is 3.1415926 in
    single precision, bits 0x40490FDA; 3.14159265 (constants::pi) would be
    0x40490FDB and change which propagators are summed in direct space."""
    modules = [src for src in program_sources(tl_root, "refcalc") if src.suffix == ".f90"]
    stop = next(i for i, src in enumerate(modules) if src.name == "composite_layer.f90")
    source = modules[stop].read_text()
    assert "  private\n" in source
    exposed = tmp_path / "composite_layer.f90"
    exposed.write_text(source.replace("  private\n", "  private\n  public :: pi_estimate\n", 1))
    probe = tmp_path / "print_pi_estimate.f90"
    probe.write_text(_PI_ESTIMATE_PROGRAM)
    exe = tmp_path / "print_pi_estimate"
    proc = subprocess.run(["gfortran", "-std=f2018", "-w", "-J", str(tmp_path), "-o", str(exe),
                           *map(str, modules[:stop]), str(exposed), str(probe), "-llapack", "-lblas"],
                          capture_output=True, text=True, cwd=tmp_path)
    assert proc.returncode == 0, proc.stderr[-3000:]
    out = subprocess.run([str(exe)], capture_output=True, text=True, check=True)
    assert out.stdout.strip() == _PI_ESTIMATE_BITS
