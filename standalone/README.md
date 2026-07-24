# Standalone HOP: EOS and opacity tables from Python

`libhop.so` is the HOP opacity/EOS engine built without the DISPATCH framework, driven
from Python through `tumult.hop`. It produces the same tables a `make_eos=.true.` run of
`dispatch.x` writes — `eos`, `opa`, `scat` — without launching dispatch at all.

## Build

```bash
make -C HOP/standalone -j8            # -> HOP/standalone/build/libhop.so
make -C HOP/standalone clean          # rebuild from scratch after editing HOP/*.f90
```

gfortran only (the framework singletons are replaced by the `-DHOP_STANDALONE` branch of
`hop_env_mod.f90`). `OPT=-O2` overrides the default `-O1`. The Python side finds the
library at `HOP/standalone/build/libhop.so`, or wherever `$HOP_LIB` points.

Check it loads:

```bash
python -c "import tumult; print(tumult.HopEOS(feh=0.0).load_state([5777.], rho=[2.5e-7]).ne)"
```

## EOS table

`make_eos_table` solves the EOS on a set of `(T, rho)` points and returns the six planes
`(temp, rho, ne, pg, E, chi500)` as float32 `(6,) + temp.shape`.

```python
import numpy as np
import matplotlib.pyplot as plt
from tumult.hop import HopEOS, logstep, make_eos_table

tab_temp = logstep(2500.0, 40000.0, 0.05)  # ln-space node spacing
tab_rho = logstep(1e-12, 1e-5, 0.25)  # ln-space node spacing
temp, rho = np.meshgrid(tab_temp, tab_rho, indexing='ij')         # both (ntemp, nrho)

with HopEOS(feh=0.0) as hop:                     # one instance per process
    eos = make_eos_table(hop, temp, rho)
_, _, ne, pg, energy, chi500 = eos

fig, axes = plt.subplots(1, 3, figsize=(13, 3.6), constrained_layout=True)
for ax, plane, label in zip(axes, (ne, pg, chi500),
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

## Opacity table

`make_opac_table` reproduces `generate_opacity_table` + `save_opac_table` on the same
`temp`/`rho` points: `opa` is the total absorption chi/rho [cm^2/g], `scat` is
(Thomson + Rayleigh)/rho, both float32 `(nlam,) + temp.shape`. It also returns the EOS
table, so this one call replaces the whole dispatch run.

```python
from tumult.hop import make_opac_table, wavelength_grid

def main():
    lam = wavelength_grid(3000.0, 20000.0, R=300)   # &spectrum_params; air -> vacuum
    out = make_opac_table(temp, rho, lam, nchunks=8, nprocs=8, feh=0.0)
    opa, scat, eos = out['opa'], out['scat'], out['eos']

    fig, ax = plt.subplots(figsize=(7, 4.2), constrained_layout=True)
    j = np.argmin(np.abs(tab_rho - 1e-7))
    for t in (4000, 6000, 10000):
        i = np.argmin(np.abs(tab_temp - t))
        ax.loglog(lam, opa[:, i, j], label=f'absorption, T = {tab_temp[i]:.0f} K')
    ax.loglog(lam, scat[:, i, j], 'k--', lw=1,
              label=f'scattering, T = {tab_temp[i]:.0f} K')
    ax.set(xlabel=r'$\lambda_{\rm vac}$ [AA]', ylabel=r'$\chi/\rho$ [cm$^2$ g$^{-1}$]',
           title=rf'$\rho$ = {tab_rho[j]:.1e} g cm$^{{-3}}$')
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
  It does not affect the EOS or the continuum.

### Line-list opacity

Give `linelist` a file, a list of files, or a folder (the `&linelist_params`
`linelist_file` / `linelist_folder`), and `dlam` [AA] as the line-wing search radius:

```python
out = make_opac_table(temp, rho, lam, linelist='input/linelists/master_line.txt',
                      dlam=100.0, vmic=1.0, nchunks=8, nprocs=8)
```

For a single state rather than a table, `HopEOS.init_linelist(...)` then
`hop.absorption(lam, include_linelist=True, dlam=...)` gives chi [1/cm] per cell.

## Writing dispatch-layout files

With `folder=`, the chunks go straight to disk in the dispatch save layout
(`opa_<iz><chunk>.bin` + `opa_meta.txt`, likewise `scat` and `eos`), so the parent never
holds the whole table and the usual readers work unchanged:

```python
make_opac_table(temp, rho, lam, folder='data/my_table/save', nprocs=8)

from tumult.io_binary import load_patch_files, concat_patches
dct, _ = load_patch_files('data/my_table/save', 'opa')
opa = concat_patches(dct, n_nu=8, fdim=0, zdim=2, lazy=False)[0]   # (nlam, ntemp, nrho)
```