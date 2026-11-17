"""Sampled line-opacity tables (hopacity.absmet) built with the standalone HOP library."""
import numpy as np
import pytest

from hop_reference import TESTS_DIR, hop_library

pytestmark = pytest.mark.skipif(
    hop_library() is None,
    reason="libhop.so not built - run 'make hop_standalone'")

LYMAN_ALPHA = """'H I  lambda(AA) nblo nbup elo(eV) eup(eV) loggf ident raddamp'
'1'   1
  1215.671     1    2    0.000   10.199   -0.080   'Lyalpha  ' 4.70E+08
"""


def test_every_marcs_node_has_a_density():
    """Every (T, P6) node of the MARCS grids is a gas state of HOP's equation of state.

    At 22000 K the EOS turns over at log rho = 0.5, far above the densities the nodes need
    (log P6 = 9 is reached below log rho = 0), and fails at 15000 K there; the density
    scan stops past the highest node so that this does not count against the nodes.
    """
    from hopacity import HopEOS
    from hopacity.absmet import MARCS_LOG_P6, MARCS_TEMP, p6_density

    with HopEOS(feh=-3.0, alpha_enhancement=0.4) as hop:
        log_rho = p6_density(hop, MARCS_TEMP, MARCS_LOG_P6)
        np.testing.assert_allclose(hop.log_p6, np.broadcast_to(MARCS_LOG_P6, log_rho.shape),
                                   atol=0.01)
    assert np.all(np.diff(log_rho, axis=1) > 0)


def test_lyman_alpha_without_electrons(tmp_path):
    """In cool metal-poor gas the Stark profile of Lyman alpha underflows at the -4000 cm^-1
    detuning where its red wing is matched to Allard et al.; the match is then skipped
    instead of turning the line into NaN."""
    from hopacity import HopEOS

    linelist = tmp_path / "Hlinedata"     # the name selects the Hlinedata column layout
    linelist.write_text(LYMAN_ALPHA)
    lam = np.linspace(1215.6, 1215.8, 201)
    with HopEOS(feh=-3.0) as hop:
        hop.load_state([500.0, 500.0], rho=[1e-10, 1e-9])    # n_e = 1e-10 to 1e-9 cm^-3
        chi, _ = hop.opacity(lam, linelist=linelist, dlam=100.0)
    assert np.all(np.isfinite(chi))


def test_table_round_trip(tmp_path):
    """A table written and read back keeps its grids, abundances and float32 opacities."""
    from hopacity.absmet import build_absmet, read_absmet, write_absmet

    lam = np.linspace(5160.0, 5170.0, 101)
    table = build_absmet(TESTS_DIR / "fixtures" / "master_line_5000_5170.txt",
                         temp=[4000.0, 6300.0], log_p6=[2.0, 4.0], lam=lam, dlam=5.0,
                         feh=-1.0)
    assert table.log_kappa.shape == (101, 2, 2)
    assert table.abundances[25] == pytest.approx(7.51 - 1.0, abs=1e-6)   # Fe, Magg et al. 2022
    assert table.log_kappa.max() > -10       # the Mg b and Fe lines
    back = read_absmet(write_absmet(tmp_path / "absmet.bin", table))
    np.testing.assert_allclose(back.lam, lam, rtol=1e-15)
    np.testing.assert_allclose(back.log_p6, table.log_p6, atol=1e-6)
    np.testing.assert_array_equal(back.temp, table.temp)
    np.testing.assert_array_equal(back.log_kappa, table.log_kappa.astype(np.float32))
    np.testing.assert_array_equal(back.abundances, table.abundances.astype(np.float32))
