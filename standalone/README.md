# Standalone HOP: EOS and opacity from Python

`libhop.so` is the HOP opacity/EOS engine built without the DISPATCH framework, driven
from Python through `tumult.hop`. It computes the same EOS state and the same opacities a
`make_eos=.true.` run of `dispatch.x` tabulates — without launching dispatch at all.

## Build

```bash
make -C HOP/standalone -j8            # -> HOP/standalone/build/libhop.so
make -C HOP/standalone clean          # rebuild from scratch after editing HOP/*.f90
```

gfortran only (the framework singletons are replaced by the `-DHOP_STANDALONE` branch of
`hop_env_mod.f90`). `OPT=-O2` overrides the default `-O1`. The Python side finds the
library at `HOP/standalone/build/libhop.so`, or wherever `$HOP_LIB` points.

Every source the library compiles lives under `HOP/` — nothing is pulled in from the
experiment's `SRC/`, so the directory can be copied or vendored on its own.

Check it loads:

```bash
python -c "import tumult; print(tumult.HopEOS(feh=0.0).load_state([5777.], rho=[2.5e-7]).ne)"
```

## The EOS

`load_state` solves the EOS on a set of `(T, rho)` points — any shape, as long as they
share it, since the cells are independent — and the state is then readable as the
properties `temp`, `rho`, `ne`, `pe`, `pg`, `energy`, `eh` and `chi500`, each in the shape
you passed in.

```python
import numpy as np
import matplotlib.pyplot as plt
from tumult.hop import HopEOS, logstep

fig, axes = plt.subplots(1, 3, figsize=(13, 3.6), constrained_layout=True)

tab_temp = logstep(2500.0, 40000.0, 0.05)  # ln-space node spacing
tab_rho = logstep(1e-12, 1e-5, 0.25)  # ln-space node spacing
temp, rho = np.meshgrid(tab_temp, tab_rho, indexing='ij')         # both (ntemp, nrho)

with HopEOS(feh=0.0) as hop:                     # one instance per process
    hop.load_state(temp, rho)

    for ax, plane, label in zip(axes, (hop.ne, hop.pg, hop.chi500),
                                (r'$\log n_e$ [cm$^{-3}$]', r'$\log P_g$ [dyn cm$^{-2}$]',
                                r'$\log \chi_{500}$ [cm$^{-1}$]')):
        m = ax.pcolormesh(np.log10(tab_temp), np.log10(tab_rho), np.log10(plane).T,
                          shading='nearest', cmap='magma', edgecolors='face')
        fig.colorbar(m, ax=ax, label=label)
        ax.set(xlabel=r'$\log T$ [K]', ylabel=r'$\log \rho$ [g cm$^{-3}$]')
```

`HopEOS(...)` takes the `&composition_params` settings: `feh`, `abund_file`,
`absdat_file`, `molec_file`, `tmolim`, `alpha_enhancement`, and any module option as a
keyword (`mhd_eos=False`, `ref_lam=5000.0`, `abundances='Fe 7.5'`, ...).

## Opacity

`opacity` is the one entry point for all of it — continuum, sampled metal-line background
and line list — and returns two float32 arrays of shape `<state shape> + (nlam,)`, in
1/cm:

- `chi`, the total extinction: continuum with Thomson + Rayleigh scattering folded in (as
  HOP's `scattering_absorption` does), plus the absmet background with `include_absmet`,
  plus the LTE line list with `linelist`;
- `scat`, the scattering part of it: Thomson + Rayleigh, and with `line_scattering=True`
  also the (1-eps)-weighted absmet and line-list opacity (the van Regemorter split
  `composition_type%opacity` applies).

Divide by `hop.rho[..., None]` for the per-gram opacity a dispatch `opa`/`scat` table
holds.

```python
from tumult.hop import wavelength_grid

def main():
    lam = wavelength_grid(3000.0, 3100.0, R=1000)   # &spectrum_params; air -> vacuum

    with HopEOS(feh=0.0) as hop:
        hop.load_state(temp, rho)
        chi, scat = hop.opacity(lam, linelist='input/linelists/master_line.txt',
                                nchunks=8, nprocs=8, line_scattering=True)
        over_rho = 1 / hop.rho[..., None]
        opa = chi * over_rho # (nlam, ntemp, nrho), cm^2/g
        sca = scat * over_rho

    fig, ax = plt.subplots(figsize=(7, 4.2), constrained_layout=True)
    i = 20
    j = 0
    ax.loglog(lam, opa[i, j], label=f'absorption')
    ax.loglog(lam, sca[i, j], label=f'scattering')
    ax.set(xlabel=r'$\lambda_{\rm vac}$ [AA]', ylabel=r'$\chi/\rho$ [cm$^2$ g$^{-1}$]',
            title=rf'T = {tab_temp[i]:.0f} K, $\rho$ = {tab_rho[j]:.1e} g cm$^{{-3}}$')
    ax.legend()

if __name__ == '__main__':       # required: the workers are spawned, see below
    main()
```

The plot shows the Balmer jump at 3646 AA, the Paschen jump at 8204 AA, and the H- bound-free
hump turning over into free-free beyond ~1.6 um. The example above (572 wavelengths,
56 x 65 nodes) takes about two minutes on 8 cores.

### Parallelism

The wavelength axis is split into `nchunks` pieces exactly as dispatch splits it over its
`n_nu` rt_patches, and the pieces run in `nprocs` processes. HOP keeps process-wide global
state, so each worker holds its own instance and — with a line list — reads only the lines
its own chunk needs. Cost scales with `nlam x ntemp x nrho`, so raise `nprocs` before you
widen the window.

Two things to get right:

- **Guard your script's top level with `if __name__ == '__main__':`.** The workers are
  spawned, not forked (a fork would inherit the parent's Fortran globals), so without the
  guard every worker re-runs the script on import and the pool spawns forever. Notebooks
  and pytest are unaffected.
- **Pass `vmic` if your atmosphere has microturbulence.** `&opac_table_params` defaults it
  to -1 (thermal broadening only), but `read_external_mod` overwrites `opac_table_vmic`
  with the mean of the atmosphere's own vmic whenever any of it is positive. Getting this
  wrong leaves the line opacity up to an order of magnitude off in the cool, dense corner.
  It does not affect the EOS or the continuum, and is applied without re-solving the state.

### Line-list opacity

Give `linelist` a file, a list of files, or a folder (the `&linelist_params`
`linelist_file` / `linelist_folder`), and `dlam` [AA] as the line-wing search radius:

```python
chi, scat = hop.opacity(lam, linelist='input/linelists/master_line.txt',
                        dlam=100.0, vmic=1.0, line_scattering=True,
                        nchunks=8, nprocs=8)
```

Only the lines within `dlam` of the requested wavelengths are read. A line list can only
ever be *added* to a HOP instance, so calling `opacity` with a different one (or over a
wider window than before) rebuilds the instance and re-solves the state — repeated calls
with the same list and window cost nothing extra. Note that chunking is therefore not
quite neutral for line opacity: a chunk sees only the lines within `dlam` of its own
wavelengths, so a far Lorentz wing reaching in from outside that window is missing —
exactly as it is in a dispatch run, where each rt_patch reads its own chunk.
