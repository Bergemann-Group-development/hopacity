"""The standalone HOP library (libhop.so) through the hopacity wrapper.

Needs the library of the imported hopacity: an installed wheel, or a source tree after
``make -C standalone``. Chunking is only
neutral within one process: a spawned worker reads only the lines within dlam of its own
chunk, so a far Lorentz wing reaching in from outside that window is missing (as in a
Dispatch run). Hence nprocs=1 here.

References ``golden/hop_eos.npz`` and ``hop_opacity.npz``
------------------------------------------------------------
Purpose: catch unintended changes to HOP's equation of state, continuum and line-list
opacity on a 4 x 6 (T, rho) grid (3500-12000 K, 1e-10 to 1e-7 g/cm^3) over 5160-5170 AA
(vacuum, R = 20000), vmic = 1 km/s. Densities are in 1/cm^3, opacities in 1/cm. There is
no analytic answer for these states, so the reference records agreement with a reviewed
result, not correctness.

``hop_eos`` was generated at df4c20bf, ``hop_opacity`` at 6a71ff283 of the Dispatch
repository, both with gfortran 15.1.0, ``OPTS=hoppe`` (-O3 -mcpu=apple-m1), from (md5
prefix): absdat/TS_absdat.txt 5b921bc86ca0, abund/abund_magg22.txt abcd0926adda,
isotopes/isotopes_blp.txt d63a83de4230, OpenMULTI3D's master_line.txt 2821ea404b92 (of which
``fixtures/master_line_5000_5170.txt`` is an exact cut, see ``hop_reference.py``), and
molec/molec_data_BC.txt a4e2d3274cda (``hop_eos``) or cfb2d3e18c35 (``hop_opacity``).

Tolerance rtol = 1e-5. The references are exact on the generating build. A change of the
optimisation flags alone (-O1 to -O3) moved line-opacity cells by float32 roundoff of at
most 5e-6 relative, so a tighter bound would fail on a different build without any change
in the physics. Changes to the physics or input data below 1e-5 go undetected.

Regenerate with ``python tests/hop_reference.py hop_eos hop_opacity`` (review the
summary), then ``--write``, and give the reason in the commit message.
"""
import numpy as np
import pytest

from hop_reference import (HOP_DLAM, HOP_LINELIST, HOP_VMIC, TESTS_DIR, assert_matches_golden,
                           hop_eos_reference, hop_grid, hop_library, hop_opacity_reference,
                           hop_wavelengths, load_golden)

pytestmark = pytest.mark.skipif(
    hop_library() is None,
    reason="libhop.so not built - run 'make hop_standalone'")

RTOL = 1e-5


def test_eos_matches_reference():
    assert_matches_golden(hop_eos_reference(), load_golden("hop_eos"), rtol=RTOL)


def test_opacity_matches_reference():
    assert_matches_golden(hop_opacity_reference(), load_golden("hop_opacity"), rtol=RTOL)


def test_chunking_is_exact_without_a_linelist():
    """Continuum and sampled-metal opacity do not depend on the wavelength grid.

    Each wavelength is computed independently there, so slicing the axis is exact.  This
    is the half of the chunking behaviour that a threaded hop_opacity may rely on.
    """
    from hopacity import HopEOS

    temp, rho = hop_grid()
    lam = hop_wavelengths()
    # the solar table of the OpenMULTI3D checkout HOP sits in
    absmet = TESTS_DIR.parents[1] / 'input' / 'absmet_files' / 'absmet'
    if not absmet.is_file():
        pytest.skip('solar absmet table not present')
    with HopEOS(feh=0.0, absmet_file=absmet) as hop:
        hop.load_state(temp, rho=rho)
        for kw in ({}, dict(include_absmet=True)):
            chi1, scat1 = hop.opacity(lam, nchunks=1, **kw)
            chi3, scat3 = hop.opacity(lam, nchunks=3, **kw)
            np.testing.assert_array_equal(chi3, chi1)
            np.testing.assert_array_equal(scat3, scat1)


def test_chunking_perturbs_line_opacity():
    """Line opacity DOES depend on the wavelength grid it is handed - pinned here.

    ``f0 = find_closest(sp%lam, lam0)`` (bbtrans_mod.f90:614) is the grid index nearest
    the line centre, and both the skip test and the analytic profile extent are keyed on
    the continuum *at that index* (linelist_mod.f90:373,389,419)::

        if (arg <= logref(f0,i,j,k)) cycle
        ref  = tol * cont(f0,i,j,k) / bi
        vm2  = max(real(arg - logref(f0,i,j,k)), a*Rij*over_pi/ref)

    A narrower chunk puts f0 on a different wavelength, so ``ref`` shifts and every line
    is truncated somewhere slightly different.  The effect sits at the ``tolerance``
    level, and widening dlam does not remove it (it is not a line-selection effect).

    Consequence for the threading work: the wavelength axis must NOT be split across
    threads, or results would depend on OMP_NUM_THREADS.  Thread the cell loop instead -
    cont/logref are indexed (f0,i,j,k), so per-cell work is untouched by that split.
    """
    from hopacity import HopEOS

    temp, rho = hop_grid()
    lam = hop_wavelengths()
    with HopEOS(feh=0.0) as hop:
        hop.load_state(temp, rho=rho)
        kw = dict(linelist=str(HOP_LINELIST), dlam=HOP_DLAM, vmic=HOP_VMIC,
                  line_scattering=True)
        chi1, _ = hop.opacity(lam, nchunks=1, **kw)
        chi3, _ = hop.opacity(lam, nchunks=3, **kw)
    assert np.any(chi3 != chi1), "chunking artifact vanished - re-read the docstring"
    rel = np.abs(chi3 - chi1) / np.abs(chi1)
    assert rel.max() < 1e-3, f"chunking artifact grew to {rel.max():.2e}"


def test_linelist_rebuild_is_transparent():
    """Widening the line-list window silently rebuilds and re-solves - same answer.

    A line list can only ever be *added* to a HOP instance, so asking for a wider window
    than was read makes the wrapper free the instance, re-init it from its recipe, replay
    the kept species and re-solve the state (_rebuild/_replay_state).  Nothing about that
    may be visible in the result: reading a wide window directly must give the same
    numbers as narrowing to it after a narrow read.

    This is the contract to preserve if the rebuild ever moves into the C shim.
    """
    from hopacity import HopEOS

    temp, rho = hop_grid()
    lam = hop_wavelengths()
    kw = dict(linelist=str(HOP_LINELIST), vmic=HOP_VMIC, line_scattering=True)

    with HopEOS(feh=0.0) as hop:                      # wide window from the start
        hop.load_state(temp, rho=rho, species=["Fe I"])
        direct, _ = hop.opacity(lam, dlam=30.0, **kw)
        n_fe_direct = hop.species_density("Fe I")

    with HopEOS(feh=0.0) as hop:                      # narrow first, then widen
        hop.load_state(temp, rho=rho, species=["Fe I"])
        hop.opacity(lam, dlam=5.0, **kw)
        rebuilt, _ = hop.opacity(lam, dlam=30.0, **kw)
        n_fe_rebuilt = hop.species_density("Fe I")

    np.testing.assert_array_equal(rebuilt, direct)
    np.testing.assert_array_equal(n_fe_rebuilt, n_fe_direct)


def test_reopen_after_close():
    """HOP keeps process-global state; closing must leave it reusable.

    The reference generator and the two reference tests above all open an instance in
    the same process, so this is a precondition for the rest of the module.
    """
    from hopacity import HopEOS

    seen = []
    for _ in range(2):
        with HopEOS(feh=0.0) as hop:
            seen.append(hop.load_state([5777.0], rho=[2.5e-7]).ne[0])
    assert seen[0] == seen[1]


def test_legacy_absdat_with_occupation_weighted_eos():
    """The legacy absdat gives H I four partition-function levels, fewer than the
    hydrogenic levels the occupation-weighted EOS (mhd_eos, on by default) fills.

    Filling them without resizing the array corrupts the heap and aborts the process.
    """
    from hopacity import HopEOS
    from hopacity.core import _DATA_DIR

    with HopEOS(absdat_file=_DATA_DIR / "absdat" / "MULTI_absdat.txt",
                abund_file=_DATA_DIR / "abund" / "legacy_abund.txt",
                molec_file=_DATA_DIR / "molec" / "molec_data.txt") as hop:
        hop.load_state(np.array([5000.0, 6000.0]), ne=np.array([1e13, 1e14]))
        chi, _ = hop.opacity(np.array([2000.0, 6563.0]))
    assert np.all(np.isfinite(chi)) and np.all(chi > 0)


def test_option_errors_carry_the_librarys_message():
    """hop_set_option distinguishes an unknown name from an unparseable value, and
    hop_last_error explains which - neither was possible with the old status-1-for-both
    three-setter API."""
    from hopacity import HopEOS, HopError

    with HopEOS(feh=0.0) as hop:
        with pytest.raises(HopError, match="unknown option 'nonesuch'"):
            hop.set_option("nonesuch", 1)
        with pytest.raises(HopError, match="cannot read 'maybe' as a logical"):
            hop.set_option("mhd_eos", "maybe")
        # every listed name is settable, and the value crosses as a string either way
        names = hop.option_names()
        assert {"mhd_eos", "ref_lam", "abundances", "voigt_kernel"} <= set(names)
        hop.set_option("ref_lam", 5000)      # int for a real-valued option
        hop.set_option("ref_lam", 5000.0)
        hop.set_option("mhd_eos", False)
        hop.set_option("voigt_kernel", 1)  # inert since the QLG-bot switch; the ABI still carries it


def test_species_error_names_the_species():
    from hopacity import HopEOS, HopError

    with HopEOS(feh=0.0) as hop:
        with pytest.raises(HopError, match="unknown species 'Nonium XII'"):
            hop.keep_species("Nonium XII")


def test_second_instance_is_refused():
    """The single-instance guard is what the wrapper's rebuild logic relies on."""
    from hopacity import HopEOS, HopError

    with HopEOS(feh=0.0):
        with pytest.raises(HopError, match="already open"):
            HopEOS(feh=0.0)


def test_line_strengths_are_the_ranking_hop_itself_uses():
    """hop_get_linelist hands back the per-line strengths hop_opacity worked out.

    The strength is ``max over cells of (line-centre extinction / continuum extinction)`` --
    the one line ranking that is comparable *between* species, because it carries the
    abundance, the ionization fraction, the dissociation equilibrium and the partition
    function that a gf-and-Boltzmann estimate does not. It is the same quantity a full run
    writes as the trailing column of ``linelist_out.txt``; the point of the accessor is to
    get it from a solved state in memory instead of from a stored file.

    Pinned here by what it must *order*, not by its values: a solar-type photosphere over
    5000-5002 AA is headed by Fe I 5001.87, then Ni I 5000.34, then Ti I 5001.01, and the
    molecular bands in the window sit orders of magnitude below them. Ranked by the naive
    proxy those same bands come out on top, so this is the assertion that would catch the
    accessor silently returning nothing useful.
    """
    from hopacity import HopEOS, air2vac

    lam = air2vac(np.arange(5000.0, 5002.0, 0.005))
    with HopEOS(feh=0.0) as hop:
        hop.load_state([5876.0], rho=[2.609e-07])   # tau_Ross ~ 2/3 of a solar MARCS model
        lines = hop.line_strengths(lam, str(HOP_LINELIST), dlam=2.0, vmic=HOP_VMIC)

    assert len(lines) > 0
    assert set(lines.columns) == {"ion", "lam0", "evi", "loggf", "strength"}
    assert (lines.strength >= 0).all()

    inside = lines[(lines.lam0 >= 5000.0) & (lines.lam0 <= 5002.0)]
    top = inside.sort_values("strength", ascending=False)
    assert list(top.ion.head(3)) == ["Fe I", "Ni I", "Ti I"]
    assert abs(top.lam0.iloc[0] - 5001.87) < 0.02

    # the molecules present in this window are far below the atomic lines, which is the
    # whole distinction the proxy cannot make
    molecules = inside[inside.ion.isin(["C2", "CN", "CH", "MgH"])]
    if len(molecules):
        assert molecules.strength.max() < 0.05 * top.strength.iloc[0]


def test_line_strengths_need_a_state():
    """The strengths are a by-product of an opacity call, so there has to be a state."""
    from hopacity import HopEOS, HopError

    lam = np.array([5000.0, 5001.0])
    with HopEOS(feh=0.0) as hop:
        with pytest.raises(HopError):
            hop.line_strengths(lam, str(HOP_LINELIST))


def test_linelist_numeric_reader_keeps_supported_and_fallback_syntax(tmp_path):
    """The fast reader covers normal decimals/exponents and preserves list-I/O fallback."""
    from hopacity import HopEOS, air2vac

    lines_path = tmp_path / "numeric-forms.txt"
    lines_path.write_text(
        "'Fe I'\n"
        "'0' 3\n"
        "5000.000 0.0 2.0\n"
        "5000.500\t0.0\t2.0D+00 0 1.0E+08 0 1 2\n"
        "5001.000, 0.0, 2.0\n"
    )
    lam = air2vac(np.arange(4999.9, 5001.1, 0.005))
    with HopEOS(feh=0.0) as hop:
        hop.load_state([5876.0], rho=[2.609e-07])
        lines = hop.line_strengths(lam, lines_path, dlam=0.2, vmic=HOP_VMIC)

    inside = lines[(lines.lam0 >= 4999.9) & (lines.lam0 <= 5001.1)]
    np.testing.assert_allclose(inside.lam0, [5000.0, 5000.5, 5001.0], atol=2e-8)
    np.testing.assert_allclose(inside.loggf, 2.0)


MG_B = """\
'  12.000            '    1         2
'Mg I '
  5172.684  2.712 -0.393  8.00   3.0  1.0E+08 'x' 'x'  0.0  1.0 'Mg I'
  5183.604  2.717 -0.167  8.00   5.0  1.0E+08 'x' 'x'  0.0  1.0 'Mg I'
"""


def test_long_data_paths(tmp_path):
    """The bundled data files sit wherever the package is installed, and a virtual
    environment in a deep directory gives paths beyond 256 characters.

    The tables are read once per process, so both runs go to fresh worker processes.
    """
    from hopacity import HopEOS
    from hopacity.core import _DATA_DIR

    deep = tmp_path.joinpath(*[f"padding_directory_{i:02d}" for i in range(20)])
    deep.mkdir(parents=True)
    (deep / "Anstee.txt").symlink_to(_DATA_DIR / "Anstee-1802-2014.txt")
    (deep / "vcs").symlink_to(_DATA_DIR / "vcs")
    assert len(str(deep / "Anstee.txt")) > 400
    linelist = tmp_path / "mgb.txt"
    linelist.write_text(MG_B)
    lam = np.linspace(5170.0, 5186.0, 40)
    chi = []
    for paths in ({}, {"anstee_file": str(deep / "Anstee.txt"), "vcs_data_dir": str(deep / "vcs")}):
        with HopEOS(feh=0.0, **paths) as hop:
            hop.load_state([5000.0], rho=[1e-8])
            chi.append(hop.opacity(lam, linelist=linelist, dlam=10.0, vmic=1.0, nprocs=2)[0])
    np.testing.assert_array_equal(chi[1], chi[0])


def test_a_dying_worker_raises(tmp_path):
    """HOP's error stop ends a worker process; the parent must report it, not wait."""
    from hopacity import HopEOS, HopError

    linelist = tmp_path / "mgb.txt"
    linelist.write_text(MG_B)
    with HopEOS(feh=0.0, anstee_file=str(tmp_path / "missing.dat")) as hop:
        hop.load_state([5000.0], rho=[1e-8])
        with pytest.raises(HopError, match="worker process died"):
            hop.opacity(np.linspace(5170.0, 5186.0, 40), linelist=linelist, dlam=10.0,
                        nprocs=2)


LYMAN_BETA = """'H I  lambda(AA) nblo nbup elo(eV) eup(eV) loggf ident raddamp'
'1'   1
  1025.722     1    3    0.000   12.087   -0.801   'Lybeta   ' 5.70E+08
"""


def test_lyman_wing_gives_way_to_rayleigh_scattering(tmp_path):
    """Redward of 1400 AA a Lyman line's natural damping wing is H I Rayleigh scattering.

    The continuum's Rayleigh cross-section (Dalgarno) is the sum of these far red wings
    over the series, so the line profile leaves natural damping out there. Without that
    cut a constant-width Lorentz wing levels off far from the line instead of falling
    as lambda^-4: in a cool, thin cell Lyman beta alone added about 1e-3 cm^2/g of true
    absorption from 1400 AA to 20 um, several times the continuum in the infrared, and
    only within the chunk that loaded the line. Pressure broadening is negligible in
    this cell, so redward of 1400 AA the line may add only float32 roundoff, in every
    chunk.
    """
    from hopacity import HopEOS, wavelength_grid

    linelist = tmp_path / "Hlinedata"     # the name selects the Hlinedata column layout
    linelist.write_text(LYMAN_BETA)
    lam = wavelength_grid(1000.0, 30000.0, R_sampling=500, in_air=False)
    with HopEOS(feh=0.0) as hop:
        hop.load_state([3800.0], rho=[1e-11])
        cont, _ = hop.opacity(lam)
        chi1, _ = hop.opacity(lam, linelist=linelist, dlam=100.0)
        chi4, _ = hop.opacity(lam, linelist=linelist, dlam=100.0, nchunks=4)
    red = lam > 1400.0
    assert np.all(chi1[0, ~red & (lam > 1030.0)] > cont[0, ~red & (lam > 1030.0)])
    np.testing.assert_allclose(chi1[0, red], cont[0, red], rtol=1e-6)
    np.testing.assert_allclose(chi4, chi1, rtol=1e-6)
