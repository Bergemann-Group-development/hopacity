"""Sampled metal-line opacity tables (absmet) computed with HOP, and their file format.

An absmet table holds log10 of the line opacity per gram of stellar matter on a grid of
wavelength, temperature and van der Waals perturber pressure
P6 = P(H I) + 0.41 P(He I) + 0.85 P(H2), in the layout of the MARCS tables written by
``add_metal_opac.f`` (Gustafsson et al. 2008, Sect. 5.3), which HOP's ``osmet_mod``,
MULTI and Turbospectrum read.  The line profile of a metal line depends on the gas
through T and, by van der Waals broadening, P6; a node (T, P6) also fixes the gas state,
the density at which the equation of state gives that P6, and with it the ionisation and
dissociation equilibria that set the number of absorbers per gram.

:func:`build_absmet` computes such a table: at every node it solves HOP's equation of
state for that density and takes the opacity of a line list minus the continuum, at the
nodes' own wavelengths (opacity sampling).  Unless given otherwise, the grids are those
of the MARCS tables: 17 temperatures from 500 to 30000 K, log10 P6 from -8 to 9 and
153910 wavenumbers spaced by a factor 1 + 1/20000 from 910 AA to 200 um.
"""
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .core import HopEOS, HopError

__all__ = ['AbsmetTable', 'MARCS_TEMP', 'MARCS_LOG_P6', 'marcs_wavelengths',
           'read_absmet', 'write_absmet', 'p6_density', 'build_absmet']

MARCS_TEMP = np.array([500., 1000., 1500., 2000., 2500., 3300., 4000., 4800., 5600., 6300.,
                       7200., 8400., 10000., 12500., 15000., 22000., 30000.])
MARCS_LOG_P6 = np.arange(-8.0, 10.0)
LOG_KAPPA_FLOOR = -30.0   # the tables' value for no opacity

_N_SPECIES = 184          # length of the species record
_N_ELEMENTS = 92          # length of the abundance record


def marcs_wavelengths():
    """Vacuum wavelengths [AA] of the MARCS tables, ascending: wavenumbers 1e8/910 cm^-1
    divided by 1.00005 per step, 153910 of them."""
    wnum = 1e8 / 910.0 * 1.00005 ** -np.arange(153910)
    return 1e8 / wnum


@dataclass
class AbsmetTable:
    """An absmet table.  ``log_kappa`` is (wavelength, temperature, P6), log10 cm^2/g.

    A table built by :func:`build_absmet` leaves the species record empty; its first
    header record names the line list instead.
    """
    lam: np.ndarray           #: vacuum wavelengths [AA], ascending
    temp: np.ndarray          #: temperatures [K]
    log_p6: np.ndarray        #: log10 P6 [dyn/cm2]
    log_kappa: np.ndarray     #: log10 line opacity per gram [cm^2/g]
    abundances: np.ndarray    #: log eps (H = 12) by Z = 1..92
    vmic: float               #: microturbulence [km/s]
    species: np.ndarray       #: species codes 100 Z + stage - 1 (0 pads the record)
    text: tuple = ('', '')    #: the two 80-character header records


def _records(raw, endian):
    """Split the bytes of an unformatted sequential file into records."""
    out, pos = [], 0
    dt = np.dtype(endian + 'u4')
    while pos + 4 <= len(raw):
        n = int(np.frombuffer(raw, dt, 1, pos)[0])
        if pos + 8 + n > len(raw) or np.frombuffer(raw, dt, 1, pos + 4 + n)[0] != n:
            raise ValueError
        out.append(raw[pos + 4:pos + 4 + n])
        pos += 8 + n
    return out


def read_absmet(path):
    """Read an absmet table in either byte order.  The ±30 sentinels of a table are kept
    as they are."""
    raw = Path(path).read_bytes()
    for endian in '<>':
        try:
            recs = _records(raw, endian)
            break
        except ValueError:
            continue
    else:
        raise HopError(f'{path}: not an unformatted sequential file')
    i4, f4, f8 = (np.dtype(endian + t) for t in ('i4', 'f4', 'f8'))
    np6 = int(np.frombuffer(recs[4], i4)[0])
    ntemp = int(np.frombuffer(recs[6], i4)[0])
    rec = np.dtype([('w', f8), ('k', f4, (ntemp, np6))])
    body = np.frombuffer(b''.join(recs[12:]), rec)
    body = body[body['w'] != 0]
    # the file runs in decreasing wavenumber, so in increasing wavelength
    return AbsmetTable(lam=1e8 / body['w'].astype(float),
                       temp=np.frombuffer(recs[7], f4).astype(float),
                       log_p6=np.log10(np.frombuffer(recs[5], f4).astype(float)),
                       log_kappa=body['k'].astype(float),
                       abundances=np.frombuffer(recs[11], f4).astype(float),
                       vmic=float(np.frombuffer(recs[3], f4)[0]),
                       species=np.frombuffer(recs[2], i4).copy(),
                       text=tuple(r.decode('ascii', 'replace').rstrip() for r in recs[:2]))


def write_absmet(path, table):
    """Write ``table`` little-endian in the add_metal_opac.f layout, log_kappa clipped
    to ±30."""
    def record(arr):
        payload = np.ascontiguousarray(arr).tobytes()
        n = np.array([len(payload)], '<u4').tobytes()
        return n + payload + n

    nlam, ntemp, np6 = table.log_kappa.shape
    wnum = 1e8 / np.asarray(table.lam, float)
    species = np.zeros(_N_SPECIES, '<i4')
    species[:len(table.species)] = table.species
    abund = np.full(_N_ELEMENTS, -99.0)
    abund[:len(table.abundances)] = table.abundances
    rec = np.empty(nlam, [('w', '<f8'), ('k', '<f4', (ntemp, np6))])
    rec['w'] = wnum
    rec['k'] = np.clip(table.log_kappa, -30.0, 30.0)
    end = np.zeros(1, rec.dtype)          # wavenumber 0 ends the table
    with open(path, 'wb') as fo:
        for txt in table.text:
            fo.write(record(np.frombuffer(f'{txt:80.80s}'.encode('ascii'), np.uint8)))
        fo.write(record(species))
        fo.write(record(np.array([table.vmic], '<f4')))
        fo.write(record(np.array([np6], '<i4')))
        fo.write(record((10.0 ** np.asarray(table.log_p6)).astype('<f4')))
        fo.write(record(np.array([ntemp], '<i4')))
        fo.write(record(np.asarray(table.temp, '<f4')))
        fo.write(record(np.array([nlam], '<i4')))
        fo.write(record(np.array([wnum.max(), wnum.min()], '<f8')))
        fo.write(record(np.zeros(12, np.uint8)))
        fo.write(record(abund.astype('<f4')))
        for r in rec:
            fo.write(record(r))
        fo.write(record(end))
    return path


def p6_density(hop, temp, log_p6, log_rho=None):
    """log10 density [g/cm^3] at which HOP's equation of state gives log10 P6 ``log_p6``,
    for every pair of ``temp`` and ``log_p6`` (shape (ntemp, np6)).

    A scan over ``log_rho`` (default -22 to 2 in steps of 0.1) brackets each node, P6 rising with density at fixed T, and
    a linear interpolation in log P6 places it; a final solve at those densities has to
    reproduce the targets to 0.01 dex, and leaves that state loaded.  Nodes outside the
    scan raise.
    """
    log_rho = np.linspace(-22.0, 2.0, 241) if log_rho is None else np.asarray(log_rho, float)
    temp = np.asarray(temp, float)
    log_p6 = np.asarray(log_p6, float)
    tt, rr = np.meshgrid(temp, log_rho, indexing='ij')
    hop.load_state(tt, rho=10.0 ** rr)
    scan = hop.log_p6
    out = np.empty((temp.size, log_p6.size))
    for i in range(temp.size):
        # the equation of state fails at densities far above the highest node (log rho
        # > 0), so the scan ends at the first density past it
        above = np.nonzero(scan[i] >= log_p6.max())[0]
        n = above[0] + 1 if above.size else log_rho.size
        if not np.all(np.diff(scan[i, :n]) > 0):
            raise HopError(f'P6 does not rise with density at T = {temp[i]} K')
        out[i] = np.interp(log_p6, scan[i, :n], log_rho[:n], left=np.nan, right=np.nan)
    if not np.all(np.isfinite(out)):
        raise HopError('a (T, P6) node lies outside the density scan')
    tt, _ = np.meshgrid(temp, log_p6, indexing='ij')
    hop.load_state(tt, rho=10.0 ** out)
    err = np.abs(hop.log_p6 - log_p6[None, :]).max()
    if err > 0.01:
        raise HopError(f'P6 reproduced only to {err:.3f} dex')
    return out


def build_absmet(linelist, *, vmic=1.0, temp=MARCS_TEMP, log_p6=MARCS_LOG_P6, lam=None,
                 line_margin=100.0, nprocs=1, text=None, **composition):
    """Compute an absmet table with HOP.

    linelist : line list(s) as :meth:`HopEOS.opacity` takes them.  Include the hydrogen
        lines (``Hlinedata``) for a table that replaces HOP's line list in the NLTE
        iterations; HOP's continuum carries none of them, only the bound-free opacity.
    vmic : microturbulence [km/s] of the line profiles.
    temp, log_p6, lam : the grids, by default those of the MARCS tables.
    line_margin, nprocs : passed to :meth:`HopEOS.opacity`.
    text : the two header records; by default the line list and the composition.
    composition : keyword arguments of :class:`HopEOS` (``feh``, ``alpha_enhancement``,
        ``abund_file``, ``abundances``, ...).

    The scripts calling this with ``nprocs > 1`` need an ``if __name__ == '__main__':``
    guard (see :meth:`HopEOS.opacity`).  The table takes about half a minute on 8 processes
    for the MARCS grids and 870000 lines.
    """
    lam = marcs_wavelengths() if lam is None else np.asarray(lam, float)
    temp = np.asarray(temp, float)
    log_p6 = np.asarray(log_p6, float)
    with HopEOS(**composition) as hop:
        p6_density(hop, temp, log_p6)
        rho = hop.rho[..., None]
        cnt, _ = hop.opacity(lam, nprocs=nprocs)
        tot, _ = hop.opacity(lam, linelist=linelist, line_margin=line_margin, vmic=vmic, nprocs=nprocs)
        abund = hop.abundances
    kappa = np.maximum(tot.astype(float) - cnt, 0.0) / rho
    if not np.all(np.isfinite(kappa)):
        raise HopError('HOP returned non-finite line opacity')
    log_kappa = np.log10(np.maximum(kappa, 10.0 ** LOG_KAPPA_FLOOR))
    if text is None:
        llists = [linelist] if isinstance(linelist, (str, Path)) else list(linelist)
        text = ('hopacity absmet: ' + ' '.join(Path(f).name for f in llists),
                ' '.join(f'{k}={v}' for k, v in composition.items()) or 'default composition')
    return AbsmetTable(lam=lam, temp=temp, log_p6=log_p6,
                       log_kappa=np.moveaxis(log_kappa, -1, 0), abundances=abund,
                       vmic=float(vmic), species=np.zeros(0, np.int32), text=text)
