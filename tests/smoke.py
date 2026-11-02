"""Checks that an installed hopacity loads its library and bundled data files: an EOS
solve, then continuum and line opacity on a short window of a two-line list."""
import tempfile
from pathlib import Path

import numpy as np

from hopacity import HopEOS, wavelength_grid

# Two Mg I b lines in the Turbospectrum format
LINELIST = """\
'  12.000            '    1         2
'Mg I '
  5172.684  2.712 -0.393  8.00   3.0  1.0E+08 'x' 'x'  0.0  1.0 'Mg I'
  5183.604  2.717 -0.167  8.00   5.0  1.0E+08 'x' 'x'  0.0  1.0 'Mg I'
"""


def main():
    lam = wavelength_grid(5165.0, 5190.0, R_sampling=20000)
    with tempfile.TemporaryDirectory() as tmp:
        linelist = Path(tmp) / 'mgb.txt'
        linelist.write_text(LINELIST)
        with HopEOS(feh=0.0) as hop:
            hop.load_state([5000.0, 6500.0], rho=[1e-8, 3e-7])
            cont, _ = hop.opacity(lam)
            chi, scat = hop.opacity(lam, linelist=linelist, dlam=10.0, vmic=1.0)
            ne = hop.ne
    assert np.all(ne > 0), ne
    assert np.all(np.isfinite(chi)) and np.all(chi >= cont), 'line opacity below continuum'
    peak = (chi / cont).max()
    assert peak > 10, f'the Mg b lines are missing: peak line/continuum {peak:.3g}'
    print(f'hopacity smoke test passed: n_e = {ne}, peak line/continuum = {peak:.3g}')


if __name__ == '__main__':
    main()
