"""The HOP EOS and opacity states the references in ``golden/`` are computed on, and the
one path that writes those references.

Inputs are the default data files of ``hopacity.HopEOS`` and the line list
``fixtures/master_line_5000_5170.txt``, through the ``libhop.so`` of the imported
hopacity. The fixture is a verbatim cut of OpenMULTI3D's ``input/linelists/master_line.txt``
(md5 prefix 2821ea404b92): every species block with its records between 4960 and 5040 AA
and between 5120 and 5210 AA (air), counts adjusted. It covers every window the tests read
with their largest ``dlam`` (30 AA), and gives the same opacity, scattering and line
strengths there as the full list, bit for bit.

A reference records agreement with a reviewed earlier result, not correctness.
Regenerate one only after a deliberate change to the physics or its input data::

    python tests/hop_reference.py hop_opacity          # summary of what would change
    python tests/hop_reference.py hop_opacity --write  # rewrite the .npz

The summary is the review: per quantity, the largest relative change and the fraction of
values that moved. Record the reason in ``test_hop.py``'s docstring.
"""
import argparse
import sys
from pathlib import Path

import numpy as np

TESTS_DIR = Path(__file__).resolve().parent
GOLDEN = TESTS_DIR / "golden"
HOP_LINELIST = TESTS_DIR / "fixtures" / "master_line_5000_5170.txt"

# a coarse T x rho node grid (4 x 6 = 24 cells) and a narrow 10 AA window: small enough to
# solve in seconds, wide enough to cross the Mg b lines and a good deal of Fe I
HOP_DLAM = 10.0        # line-wing search radius [AA]
HOP_VMIC = 1.0         # [km/s]; must be pinned, it moves the line opacity a lot


def hop_library():
    """The libhop.so of the imported hopacity (a source tree's, or an installed wheel's),
    or None when it is not built."""
    from hopacity import HopError, find_library

    try:
        return find_library()
    except HopError:
        return None


def hop_grid():
    """The (temp [K], rho [g/cm^3]) node grid, both shaped (ntemp, nrho)."""
    from hopacity import logstep

    tab_temp = logstep(3500.0, 12000.0, 0.35)
    tab_rho = logstep(1e-10, 1e-7, 1.2)
    return np.meshgrid(tab_temp, tab_rho, indexing="ij")


def hop_wavelengths():
    """The vacuum wavelength grid [AA] of the opacity reference."""
    from hopacity import wavelength_grid

    return wavelength_grid(5160.0, 5170.0, R_sampling=20000)


def hop_eos_reference():
    """Thermodynamic state, two species densities and the (T, ne) -> rho inversion."""
    from hopacity import HopEOS

    temp, rho = hop_grid()
    with HopEOS(feh=0.0) as hop:
        hop.load_state(temp, rho=rho, species=["Fe I", "CH"])
        out = {
            "temp": hop.temp, "rho": hop.rho, "ne": hop.ne, "pe": hop.pe, "pg": hop.pg,
            "energy": hop.energy, "eh": hop.eh, "chi500": hop.chi500,
            "grph": np.array(hop.grph),
            "n_FeI": hop.species_density("Fe I"),
            "n_CH": hop.species_density("CH"),
        }
        # the (T, ne) branch of load_state must invert the (T, rho) branch
        hop.load_state(temp, ne=out["ne"])
        out["rho_from_ne"] = hop.rho
    return out


def hop_opacity_reference():
    """Continuum and line-list opacity [1/cm] and scattering on the node grid."""
    from hopacity import HopEOS

    temp, rho = hop_grid()
    lam = hop_wavelengths()
    with HopEOS(feh=0.0) as hop:
        hop.load_state(temp, rho=rho)
        chi_c, scat_c = hop.opacity(lam)
        chi_l, scat_l = hop.opacity(lam, linelist=str(HOP_LINELIST), dlam=HOP_DLAM,
                                    vmic=HOP_VMIC, line_scattering=True)
        return {
            "lam": lam,
            "chi_continuum": chi_c, "scat_continuum": scat_c,
            "chi_linelist": chi_l, "scat_linelist": scat_l,
        }


REFERENCES = {"hop_eos": hop_eos_reference, "hop_opacity": hop_opacity_reference}


def load_golden(name):
    """``golden/<name>.npz`` as a plain dict."""
    path = GOLDEN / f"{name}.npz"
    if not path.exists():
        raise FileNotFoundError(f"{path} is missing; write it with "
                                f"'python tests/hop_reference.py {name} --write'")
    with np.load(path, allow_pickle=False) as data:
        return {k: data[k] for k in data.files}


def assert_matches_golden(actual, golden, rtol):
    """Compare a computed reference dict with its golden, key by key."""
    missing = set(golden) - set(actual)
    extra = set(actual) - set(golden)
    assert not missing and not extra, f"golden keys drifted: missing={missing} extra={extra}"
    for key in sorted(golden):
        np.testing.assert_allclose(np.asarray(actual[key]), golden[key], rtol=rtol,
                                   err_msg=f"golden mismatch for {key!r}")


def summarise(old, new):
    """One line per quantity: largest relative change and the fraction of values moved."""
    lines = []
    for key in sorted(set(old) | set(new)):
        if key not in old or key not in new:
            lines.append(f"{key:24s} {'added' if key in new else 'removed'}")
            continue
        a, b = old[key], new[key]
        if a.shape != b.shape:
            lines.append(f"{key:24s} shape {a.shape} -> {b.shape}")
            continue
        moved = a != b
        if not moved.any():
            lines.append(f"{key:24s} identical")
            continue
        rel = np.abs(b - a) / np.maximum(np.abs(a), np.finfo(float).tiny)
        lines.append(f"{key:24s} max rel {rel.max():.2e}, median moved {np.median(rel[moved]):.2e}, "
                     f"{moved.mean():.1%} of {a.size} moved")
    return lines


def main(argv=None):
    parser = argparse.ArgumentParser(description="Review or rewrite HOP's golden references.")
    parser.add_argument("names", nargs="+", choices=sorted(REFERENCES))
    parser.add_argument("--write", action="store_true", help="rewrite the .npz files")
    args = parser.parse_args(argv)
    for name in args.names:
        new = {k: np.asarray(v) for k, v in REFERENCES[name]().items()}
        path = GOLDEN / f"{name}.npz"
        print(f"== {name}")
        if path.exists():
            print("\n".join("   " + line for line in summarise(load_golden(name), new)))
        else:
            print("   (no existing reference)")
        if args.write:
            np.savez_compressed(path, **new)
            print(f"   wrote {path}")


if __name__ == "__main__":
    sys.exit(main())
