"""Turbospectrum-format line lists read by HOP (linelist_reader_mod, linelist_mod).

A Turbospectrum block opens with its species code (``'26.000' 1 4``, ``'0106.000012' 1 2``).
HOP reads fdamp, raddmp and gamst as the tumult gvdw, grad and gstk, and resolves the
Anstee-Barklem-O'Mara cross sections Turbospectrum would (anstee_mod), so the
linelist_out it writes in the tumult layout must give the opacity the list itself gives.
The Fe I records are from the Gaia-ESO NLTE list of Turbospectrum_NLTE (4 Sep 2026).

Needs the library of the imported hopacity and its default data, the ABO table included.
"""
import subprocess
import sys

import numpy as np
import pytest

from hop_reference import HOP_VMIC, hop_library

pytestmark = pytest.mark.skipif(
    hop_library() is None,
    reason="libhop.so not built - run 'make -C standalone'")

FE_I = """'  26.000            '    1    4
'Fe I    NLTE'
  5161.738 4.294 -5.246 -7.510 7.0 2.09E+08 0.000 's' 'p' 0.0 1.0 'Fe I b1D t5P*'  0 0
  5162.272 4.178 0.020 798.279 11.0 2.40E+08 0.000 'p' 'd' 0.0 1.0 'Fe I y5F* g5F'  0 0
  5162.376 2.609 -4.592 -7.710 3.0 4.47E+08 0.000 's' 'p' 0.0 1.0 'Fe I b3F x5D*'  0 0
  5163.383 3.929 -4.613 -7.550 3.0 2.14E+08 0.000 'p' 'd' 0.0 1.0 'Fe I z3D* f7D'  0 0
"""


def _ch(code):
    return f"""'{code} '    1    2
'CH test'
  5163.000 0.580 -1.261 .00 40. 0.177E+07  'X' 'X'   0.0    1.0 'x'
  5164.000 0.600 -1.000 .00 40. 0.177E+07  'X' 'X'   0.0    1.0 'x'
"""


def _wavelengths():
    from hopacity import air2vac
    return air2vac(np.arange(5160.0, 5166.0, 0.01))


def _opacities(lists, write=None):
    from hopacity import HopEOS
    lam = _wavelengths()
    with HopEOS() as hop:
        hop.load_state([4500.0, 6000.0], rho=[1e-7, 3e-7])
        if write is not None:
            hop.write_linelist(write, lam, lists[0], dlam=2.0, vmic=HOP_VMIC)
        return [hop.opacity(lam, linelist=ll, dlam=2.0, vmic=HOP_VMIC)[0] for ll in lists]


def test_tumult_rewrite_of_a_turbospectrum_list_gives_its_opacity(tmp_path):
    """linelist_out states the values HOP resolved, so it reproduces the list's opacity."""
    ts = tmp_path / 'ts.list'
    ts.write_text(FE_I + _ch('0106.000012'))
    out = tmp_path / 'linelist_out.txt'
    chi_ts, = _opacities([ts], write=out)
    chi_out, = _opacities([out])
    np.testing.assert_allclose(chi_out, chi_ts, rtol=1e-4)

    text = out.read_text().splitlines()
    assert "'000012' 2" in text                    # 1H 12C, three digits per constituent
    gvdw = {round(float(r.split()[0]), 1): float(r.split()[3])
            for r in text if r.strip() and not r.lstrip().startswith("'")}
    # an s-p line with a log gamma6 gets Anstee-Barklem-O'Mara values, as in Turbospectrum
    assert gvdw[5162.4] > 20
    # a packed sigma.alpha from the list is kept
    assert gvdw[5162.3] == pytest.approx(798.279)


def test_isotopologue_labels_of_both_layouts_agree(tmp_path):
    """Hydrogen has no isotope table: an explicit 001 reads as the natural mixture."""
    lists = []
    for name, text in [('natural', _ch('0106.000012')), ('explicit', _ch('0106.001012')),
                       ('tumult', "'CH'\n'000012' 2\n  5163.000 0.580 -1.261 0 1.77e6 0\n"
                                  "  5164.000 0.600 -1.000 0 1.77e6 0\n")]:
        lists.append(tmp_path / name)
        lists[-1].write_text(text)
    chi = _opacities(lists)
    assert np.max(chi[0]) > np.min(chi[0]) * 1.01   # the CH lines are there
    np.testing.assert_array_equal(chi[1], chi[0])
    np.testing.assert_array_equal(chi[2], chi[0])


def test_ragged_block_is_not_taken_for_fixed_width(tmp_path):
    """Records whose lengths average the first one's put a fixed stride on the block's end.

    Here the stride also meets the last record, and a windowed read taken as fixed-width
    then starts on the closing quote of the record before, which no Turbospectrum record
    parses (HOP stops). Run in a child process, since that stop ends the process.
    """
    pad = [0, 2, -2, 2, -2, 0]                      # lengths 72+p, summing to 6*72
    records = []
    for k, p in enumerate(pad):
        rec = f"  {5160.5 + 0.5 * k:.3f} 2.609 -2.000 -7.710 3.0 4.47E+08 'x' 'y' 0.0 1.0 '"
        records.append(rec + ' ' * (70 + p - len(rec)) + "'")
    lengths = [len(r) + 1 for r in records]
    assert sum(lengths) == len(pad) * lengths[0] and len(set(lengths)) > 1
    path = tmp_path / 'ragged.list'
    path.write_text("'  26.000 '  1  6\n'Fe I'\n" + '\n'.join(records) + '\n' + FE_I)

    script = f"""
import numpy as np
from hopacity import HopEOS, air2vac
with HopEOS() as hop:
    hop.load_state([5000.0], rho=[1e-7])
    lines = hop.line_strengths(air2vac(np.arange(5161.45, 5161.6, 0.01)), {str(path)!r}, dlam=0.05)
print(sorted(round(x, 3) for x in lines.lam0[lines.loggf == -2]))
"""
    result = subprocess.run([sys.executable, '-c', script], capture_output=True, text=True,
                            cwd=tmp_path, timeout=120)
    assert result.returncode == 0, result.stderr[-2000:] + result.stdout[-2000:]
    assert result.stdout.split('\n')[-2] == '[5161.5]'
