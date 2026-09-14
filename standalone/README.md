# Standalone HOP: EOS and opacity from Python, and line-list pruning

`libhop.so` is the HOP opacity/EOS engine built without the DISPATCH framework, driven
from Python through `tumult.hop`. It computes the same EOS state and the same opacities a
`make_eos=.true.` run of `dispatch.x` tabulates — without launching dispatch at all.
`hop_prune.x` is a program on the same objects that selects the transitions a synthesis
needs (see [Line-list pruning](#line-list-pruning)).

## Build

```bash
make -C HOP/standalone -j8            # -> build/libhop.so, libhop.h, build/hop_prune.x
make -C HOP/standalone clean          # rebuild from scratch after editing HOP/*.f90
```

gfortran only. `OPT=-O2` overrides the default `-O1`. The wrapper finds the library at
`HOP/standalone/build/libhop.so`, or wherever `$HOP_LIB` points. The top-level `make`
builds `hop_prune.x` alongside `dispatch.x` when the compiler is gfortran.

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

- `chi`, the total extinction: continuum with Thomson + Rayleigh scattering folded in,
  plus the absmet background with `include_absmet`, plus the LTE line list with
  `linelist`;
- `scat`, the scattering part of it: Thomson + Rayleigh, and with `line_scattering=True`
  also the (1-eps)-weighted absmet and line-list opacity (the van Regemorter split).

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
`n_nu` rt_patches, and the pieces run in `nprocs` processes. Each worker holds its own
instance and — with a line list — reads only the lines its own chunk needs. Cost scales
with `nlam x ntemp x nrho`, so raise `nprocs` before you widen the window.

Two things to get right:

- **Guard your script's top level with `if __name__ == '__main__':`.** The workers are
  spawned, not forked, so without the guard every worker re-runs the script on import and
  the pool spawns forever. Notebooks and pytest are unaffected.
- **Pass `vmic` if your atmosphere has microturbulence.** `&opac_table_params` defaults it
  to -1 (thermal broadening only), but a dispatch run overrides that with the mean of the
  atmosphere's own vmic whenever any of it is positive — so pass that mean here. Getting
  this wrong leaves the line opacity up to an order of magnitude off in the cool, dense
  corner. It affects neither the EOS nor the continuum, and is applied without re-solving
  the state.

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
with the same list and window cost nothing extra.

## Line strengths

`line_strengths` returns one row per line — `ion`, `lam0` [AA, air], `evi` [eV], `loggf`,
`strength` — for every line HOP used over the given wavelengths:

```python
hop.load_state(T, rho=rho)                       # a whole atmosphere, not one point
lines = hop.line_strengths(lam, 'input/linelists/master_line.txt', dlam=2.0, vmic=1.0)
lines.sort_values('strength', ascending=False).head()
```

`strength` is HOP's own diagnostic: the line's centre extinction over the continuum
extinction, maximised over the cells of the loaded state
(`linelist_mod`: `ll%strength(l) = max(..., voigt_peak*Rij / cont(f0,i,j,k))`). It is the
same quantity a full `dispatch.x` run writes as the trailing column of `linelist_out.txt`,
and it is the criterion HOP itself uses to decide which lines are worth computing.

**Why it is the number to rank by.** It carries the element abundance, the ionization
fraction, the molecular dissociation equilibrium and the partition function. A
`lam0 * 10**loggf * exp(-evi*EK/T)` estimate carries none of them, so it cannot compare two
species: over 5000–5002 AA in the Sun it puts TiO and C2 bands at the top of a window whose
real order is Fe I 5001.87, Ni I 5000.34, Ti I 5001.01 — a rank correlation with the truth of
0.2, no better than chance. Against `data/sed_3000-9000`'s own dump, with that run's line
list, `line_strengths` matches all 155 lines in the window at Spearman 0.978, with a
near-constant 1.21 offset from using a 1D MARCS structure in place of the run's own.

**Cost**, for a single-file list on a 56-layer model: 0.08 s for 20 AA, 0.32 s for 100 AA,
1.5 s for 500 AA. A folder of per-species lists is roughly 5x that, dominated by reading the
files. Sample finely enough to resolve a line core (0.005 AA is safe in the optical); a
coarser grid under-reads the peaks.

The strengths are a by-product of computing opacity, so `line_strengths` runs one `opacity`
call and then reads them back. They are zero before that call, and they *accumulate* across
calls that share a loaded line list — a line's value is its own peak either way, but the set
of lines returned grows to everything used since the list was read. Load a band once and
compute over it in one call if you want exactly that band.

The three C entry points behind it are `hop_linelist_count` (collects and sizes),
`hop_get_linelist` (the columns, with a species index) and `hop_linelist_species` (the
newline-separated species ids that index refers to). Call them in that order.

## Line-list pruning

`hop_prune.x input.nml` writes the transitions a synthesis on one atmosphere would keep
from its configured line list, as a line list in the ordinary input format. A ladder of
runs on that atmosphere (a curve of growth, an abundance fit) then reads a few hundred
lines each as its `linelist_file` instead of re-reading and re-selecting the whole master
list. The program needs neither Python nor the DISPATCH runtime.

The input reuses the run's own namelist groups, so one set of settings describes the
pruning and the runs it feeds, plus a small `&prune_params` group:

```fortran
&atmos_params        atmos_format='marcs' atmos_file='model.mod' vmic=1.0 use_abnd=.true. /
&composition_params  abundances='Fe=+1.0' relative_abundances=.true. /
&linelist_params     linelist_file='master.txt' dlam=5.0 /
&spectrum_params     window_file='windows.txt' R=300000.0 /
&prune_params        output_file='pruned.txt' sample_count=8 /
```

- The atmosphere is parsed by `HOP/atmos_reader_mod.f90`, shared with `dispatch.x`
  (MARCS, text 1D and MULTI 1D so far; a 3D cube is refused until a thermodynamic-space
  sampling policy exists for it). `sample_count` depths spread over the model's whole
  index range are evaluated (0 = every depth); spanning the range is what matters, since a
  transition's strength is a maximum over the states.
- The composition is applied in the run's order: the standard abundance file scaled by
  `feh`, replaced by the model's own abundances with `use_abnd`, then the explicit
  `abundances` override. Prune at the strongest rung of a ladder: an element's own lines
  grow with its abundance, so that selection covers the weaker rungs for them.
- Transitions are read only inside the windows widened by `dlam`, and the opacity over the
  run's constant-`R` grid marks the ones that clear HOP's negligibility test on the
  sampled states. `dlam` and `R` therefore have to match the runs.
- The output goes to `<output_file>.<pid>.partial` and is renamed into place, so an
  interrupted run never leaves a truncated list; `<output_file>.provenance.nml` echoes every
  input group the selection depended on, plus a `&prune_result` group with the sampled
  depths, the effective vmic and [Fe/H], the vacuum windows and the transition count. Invalid input exits with status 1 and one line on
  stderr.

Against the `linelist_out.txt` a `dispatch.x` run writes with the same inputs, the pruned
set is a superset: the program evaluates the file's own depth points, whose deepest layer
is hotter than any cell centre of dispatch's interpolated mesh, so a few high-excitation
and near-threshold transitions survive here that the run drops. The run re-selects from
the pruned list, so nothing it needs is lost. With the run's `dlam` equal to the pruning
`dlam`, spectra from the pruned and the master list agree to float32 rounding.

`tumult.prune` drives the program from the synthesis workflows: one child per uncached
atmosphere node on the workflow's CPU budget, killed on cancellation.
