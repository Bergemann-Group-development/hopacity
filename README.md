# hopacity: HOP's EOS and opacities from Python

HOP is the equation-of-state and opacity code of OpenMULTI3D. `hopacity` is its Python package. For any set of temperatures and densities it solves the equation of state (electron density, pressures, internal energy, number densities of atoms, ions and molecules) and computes the extinction coefficient at the wavelengths you choose: the continuum with Thomson and Rayleigh scattering, sampled metal-line opacity, and line opacity from a line list in Turbospectrum format. It needs neither DISPATCH nor OpenMULTI3D, and returns numpy arrays, so you can write opacity tables in whatever format your own code reads.

## Install

hopacity needs Python 3.11 or later. Install it into a virtual environment of its own, so it does not interfere with other Python projects:

```bash
python3 -m venv ~/hopacity-env
source ~/hopacity-env/bin/activate        # once in every new terminal
python -m pip install --upgrade pip
```

### From a release wheel

The [releases](https://github.com/Bergemann-Group-development/hopacity/releases) of this repository carry ready-built wheels for Linux x86_64 and macOS on Apple silicon (macOS 14 and later). They contain the compiled library, its Fortran and OpenMP runtimes and HOP's default input data, so no compiler is needed. Download the one for your platform and install it:

```bash
pip install hopacity-0.1.0-py3-none-macosx_14_0_arm64.whl     # the file you downloaded
```

### From the source code

On other platforms, or to build from a clone, you need a Fortran compiler and `make`: gfortran (`brew install gcc` on macOS, `sudo apt install gfortran make` on Debian or Ubuntu), or the Intel compilers with `COMPILER=ifort` or `COMPILER=ifx` set in the environment. In the activated environment:

```bash
git clone https://github.com/Bergemann-Group-development/hopacity.git
pip install ./hopacity
```

`pip` compiles the library (under a minute) and installs the package. A package built this way uses the compiler's runtime libraries, so it runs on machines where that compiler is installed.

pip installs numpy, cffi and pandas with hopacity. The plotting examples below also use matplotlib: `pip install matplotlib`.

### Check that it works

Run this from any directory except the clone. In the clone, Python imports its `hopacity/` source folder, which has no compiled library, instead of the installed package.

```bash
python -c "from hopacity import HopEOS; print(HopEOS(feh=0.0).load_state([5777.], rho=[2.5e-7]).ne)"
```

It prints the electron density of a solar-photosphere state, `[2.02576419e+13]` (cm⁻³).

## Equation of state

`load_state` solves the EOS for temperatures `temp` [K] and densities `rho` [g/cm³]. The arrays can have any shape, as long as both have the same one: each point (a "cell") is solved independently. The results are properties of the same shape:

| Property | Quantity | Unit |
|---|---|---|
| `temp`, `rho` | temperature and mass density | K, g cm⁻³ |
| `ne` | electron number density | cm⁻³ |
| `pe`, `pg` | electron and gas pressure | dyn cm⁻² |
| `energy` | internal energy | erg g⁻¹ |
| `chi500` | extinction at 5000 Å (air) | cm⁻¹ |

`species_density('Fe I')` returns the number density [cm⁻³] of one atom, ion or molecule; pass `species=['Fe I', 'CO']` to `load_state` to keep the ones you will ask for. Instead of `rho` you can give `ne`, and `rho` is then solved for.

This example solves the EOS on a table of 56 temperatures from 2500 to 40000 K and 65 densities from 10⁻¹² to 10⁻⁵ g cm⁻³. `logstep(start, stop, step)` makes nodes evenly spaced in the natural logarithm, `step` apart.

```python
import matplotlib.pyplot as plt
import numpy as np
from hopacity import HopEOS, logstep

tab_temp = logstep(2500.0, 40000.0, 0.05)      # K
tab_rho = logstep(1e-12, 1e-5, 0.25)           # g/cm^3
temp, rho = np.meshgrid(tab_temp, tab_rho, indexing='ij')   # both (56, 65)

with HopEOS(feh=0.0) as hop:
    hop.load_state(temp, rho=rho)
    planes = {r'$\log_{10} n_e$ [cm$^{-3}$]': hop.ne,
              r'$\log_{10} P_{\rm g}$ [dyn cm$^{-2}$]': hop.pg,
              r'$\log_{10} \chi_{500}$ [cm$^{-1}$]': hop.chi500}

fig, axes = plt.subplots(1, 3, figsize=(13, 3.6), layout='constrained')
for ax, (label, plane) in zip(axes, planes.items()):
    mesh = ax.pcolormesh(np.log10(tab_temp), np.log10(tab_rho), np.log10(plane).T,
                         shading='nearest', cmap='magma')
    fig.colorbar(mesh, ax=ax, label=label)
    ax.set(xlabel=r'$\log_{10} T$ [K]', ylabel=r'$\log_{10} \rho$ [g cm$^{-3}$]')
fig.savefig('eos.pdf')
```

It takes about a second.

### Composition and options

`HopEOS` sets the chemical composition. `feh` scales the solar abundances (Magg et al. 2022 by default) of all elements heavier than Li to a metallicity [Fe/H]; `alpha_enhancement` is added to the α elements O, Ne, Mg, Si, S, Ar, Ca and Ti; `abundances='Mg=7.70 Fe=7.40'` sets single elements to absolute abundances log ε (on the scale log ε(H) = 12), or adds the values to the scaled ones if you also pass `relative_abundances=True`. `abund_file`, `absdat_file`, `isotope_file` and `molec_file` replace the bundled abundance, continuum, isotope and molecular data. Any other HOP setting is a keyword argument too, and `hop.option_names()` lists them.

HOP keeps its state in the process, so only one `HopEOS` can be open at a time. Use it in a `with` block, as here, or call `close()` before opening the next.

## Opacity

`opacity(lam)` computes the extinction at the vacuum wavelengths `lam` [Å] for every cell of the loaded state. It returns two float32 arrays of shape `<state shape> + (len(lam),)`, in cm⁻¹:

- `chi`, the total extinction: the continuum including Thomson and Rayleigh scattering, plus the lines of `linelist` if you give one, plus the sampled metal-line opacity with `include_absmet=True`;
- `scat`, its scattering part: Thomson and Rayleigh scattering, and with `line_scattering=True` also the scattering part of the line opacity (from the van Regemorter approximation of the photon destruction probability).

Divide by `hop.rho[..., None]` for the opacity per gram, in cm² g⁻¹.

`wavelength_grid(blue, red, R_sampling=R)` gives wavelengths at constant resolving power `R`. It takes the end points as air wavelengths, as line lists and spectra usually quote them, and returns vacuum wavelengths; with `in_air=False` the end points are vacuum wavelengths. `air2vac` and `vac2air` convert single wavelengths.

### Continuum

```python
import matplotlib.pyplot as plt
import numpy as np
from hopacity import HopEOS, logstep, wavelength_grid


def main():
    tab_temp = logstep(2500.0, 40000.0, 0.05)
    tab_rho = logstep(1e-12, 1e-5, 0.25)
    temp, rho = np.meshgrid(tab_temp, tab_rho, indexing='ij')
    lam = wavelength_grid(1500.0, 25000.0, R_sampling=200)     # 566 wavelengths

    with HopEOS(feh=0.0) as hop:
        hop.load_state(temp, rho=rho)
        chi, scat = hop.opacity(lam, nprocs=4)       # cm^-1, shape (56, 65, 566)
        kappa = chi / hop.rho[..., None]             # cm^2/g

    j = np.searchsorted(tab_rho, 1e-8)
    fig, ax = plt.subplots(layout='constrained')
    for t in (4500.0, 6000.0, 9000.0, 12000.0):
        i = np.searchsorted(tab_temp, t)
        ax.loglog(lam, kappa[i, j], label=f'{tab_temp[i]:.0f} K')
    ax.set(xlabel=r'$\lambda_{\rm vac}$ [$\rm\AA$]', ylabel=r'$\chi/\rho$ [cm$^2$ g$^{-1}$]',
           title=rf'$\rho$ = {tab_rho[j]:.1e} g cm$^{{-3}}$')
    ax.legend(handlelength=1.5)
    fig.savefig('continuum.pdf')


if __name__ == '__main__':      # needed with nprocs > 1, see "Parallel runs"
    main()
```

The whole table takes about two seconds on 4 cores. In the figure, the Balmer and Paschen edges near 3650 and 8200 Å rise with temperature, and at 4500 K the ultraviolet opacity of metal photoionisation dominates. The edges are rounded because high hydrogen levels merge into the continuum at these densities.

### Line lists

The master line lists that OpenMULTI3D uses are published on [MPCDF Keeper](https://keeper.mpdl.mpg.de/d/3191475aebe340638422/?p=%2Flinelists%2Fmaster_linelists). `vald_atoms.txt` (19 MB) holds the atomic lines from VALD, `Hlinedata` (10 kB) the hydrogen lines, and the other files one molecule each (about 1 GB in all). hopacity downloads them:

```bash
python -m hopacity.linelists list                                   # files and sizes
python -m hopacity.linelists fetch linelists --file vald_atoms.txt --file Hlinedata
python -m hopacity.linelists fetch linelists                        # all of them
```

`fetch` puts the files in the folder you name (here `linelists`), skips files that are already complete there, and resumes an interrupted download when run again. From Python, `hopacity.linelists.fetch('linelists', ['vald_atoms.txt'])` does the same.

Give `opacity` a line-list file, a list of files, or a folder (meaning all files in it), and with `dlam` the distance in Å from the wavelength range within which lines are read in, so that the wings of lines just outside it are included:

```python
import matplotlib.pyplot as plt
import numpy as np
from hopacity import HopEOS, logstep, wavelength_grid


def main():
    tab_temp = logstep(2500.0, 40000.0, 0.05)
    tab_rho = logstep(1e-12, 1e-5, 0.25)
    temp, rho = np.meshgrid(tab_temp, tab_rho, indexing='ij')
    lam = wavelength_grid(4840.0, 4880.0, R_sampling=100000)   # around H beta
    llist = ['linelists/vald_atoms.txt', 'linelists/Hlinedata']

    with HopEOS(feh=0.0) as hop:
        hop.load_state(temp, rho=rho)
        cnt, _ = hop.opacity(lam, nprocs=4)
        chi, scat = hop.opacity(lam, linelist=llist, dlam=10.0, vmic=1.0, nprocs=4)

    i, j = np.searchsorted(tab_temp, 6000.0), np.searchsorted(tab_rho, 1e-7)
    fig, ax = plt.subplots(layout='constrained')
    ax.semilogy(lam, chi[i, j] / cnt[i, j], color='black')
    ax.set(xlabel=r'$\lambda_{\rm vac}$ [$\rm\AA$]', ylabel=r'$\chi / \chi_{\rm cont}$',
           title=rf'$T$ = {tab_temp[i]:.0f} K, $\rho$ = {tab_rho[j]:.1e} g cm$^{{-3}}$')
    fig.savefig('lines.pdf')


if __name__ == '__main__':
    main()
```

This computes 825 wavelengths for all 3640 cells in a few seconds. The figure shows Hβ at 4862.7 Å (vacuum) and the metal lines around it, as line over continuum extinction.

`vmic` [km s⁻¹] is the microturbulence that broadens the lines; without it they have thermal widths only. It does not affect the EOS or the continuum. The first call that reads a line list reads in the lines near the wavelengths it is given; a later call with another list or over a wider range reads them again, and calls with the same list and range do not.

A line-list file is in [Turbospectrum](https://github.com/bertrandplez/Turbospectrum_NLTE) format, which [`tests/smoke.py`](tests/smoke.py) shows for two Mg I lines.

`include_absmet=True` adds the opacity of a sampled metal-line table, which the option `absmet_file` names (`HopEOS(absmet_file='path/to/absmet')`). Tables are on the same [Keeper share](https://keeper.mpdl.mpg.de/d/3191475aebe340638422/?p=%2Fabsmet_files).

### Parallel runs

`nprocs` runs the calculation in that many processes, each computing part of the wavelength range. The cost grows with the number of wavelengths times the number of cells, so for large tables set `nprocs` to the number of cores. Two things follow from the separate processes:

- A script that passes `nprocs` > 1 must put its work in a function called under `if __name__ == '__main__':`, as in the examples above. Each process imports the script, and without that guard each would start the calculation again. Jupyter notebooks do not need it.
- Each process reads the lines within `dlam` of its own part of the wavelength range. Far wings of strong lines more than `dlam` away from a part are missing there, so choose `dlam` larger than the wings you need. With `nprocs=1` all lines within `dlam` of the whole range contribute everywhere.

If HOP stops with an error inside one of these processes (for example on a missing data file), `opacity` raises `HopError`, and HOP's message is printed above the Python traceback. With `nprocs=1` the same error ends the Python process itself, after printing the message.

## Line strengths

`line_strengths(lam, linelist)` lists the lines that contribute over the wavelengths `lam` for the loaded state, as a pandas DataFrame with the columns `ion`, `lam0` (air wavelength, Å), `evi` (lower-level excitation energy, eV), `loggf` and `strength`. `strength` is the line-centre extinction over the continuum extinction, at the cell of the state where that ratio is largest. Unlike an estimate from gf value and excitation energy alone, it includes the abundance, ionisation and molecule formation, so it ranks lines of different elements and molecules against each other.

```python
from hopacity import HopEOS, wavelength_grid

lam = wavelength_grid(5165.0, 5190.0, R_sampling=100000)
with HopEOS(feh=0.0) as hop:
    hop.load_state([4500.0, 5000.0, 5800.0, 6500.0], rho=[1e-8, 5e-8, 2e-7, 3e-7])
    lines = hop.line_strengths(lam, 'linelists/vald_atoms.txt', dlam=2.0, vmic=1.0)
print(lines.sort_values('strength', ascending=False).head())
```

```text
      ion      lam0    evi  loggf       strength
16   Mg I  5183.604  2.717 -0.180  136965.578125
15   Mg I  5172.684  2.712 -0.402   83206.320312
143  Fe I  5167.489  1.485 -1.118   37440.898438
14   Mg I  5167.321  2.709 -1.030   19746.525391
163  Fe I  5171.596  1.485 -1.793    7913.413086
```

The four states roughly follow a solar photosphere. `strength` is evaluated at the line centre, so the wavelength grid need not resolve the lines; it only sets the range, widened by `dlam`, from which lines are taken. The strengths accumulate over calls with the same line list, so compute a wavelength band in one call to get exactly its lines.

## Line-list pruning

`hop_prune.x` is a Fortran program built with the library (see [Development](#development)) for OpenMULTI3D syntheses. From a model atmosphere and an OpenMULTI3D namelist, it writes the lines that a synthesis on that atmosphere uses, as a line list in the input format, so that a series of runs on the same atmosphere (a curve of growth, an abundance fit) reads a few hundred lines instead of the whole master list. It needs neither Python nor DISPATCH.

Its input is the run's own namelist groups plus `&prune_params`:

```fortran
&atmos_params        atmos_format='marcs' atmos_file='model.mod' vmic=1.0 use_abnd=.true. /
&composition_params  abundances='Fe=+1.0' relative_abundances=.true. /
&linelist_params     linelist_file='master.txt' dlam=5.0 /
&spectrum_params     window_file='windows.txt' R_sampling=300000.0 /
&prune_params        output_file='pruned.txt' sample_count=8 /
```

- The atmosphere is read with the reader `dispatch.x` uses (MARCS, text 1D and MULTI 1D; 3D models are refused). `sample_count` depths spread over the model are evaluated (0 means every depth).
- The composition is applied in the run's order: the abundance file scaled by `feh`, replaced by the model's own abundances with `use_abnd`, then `abundances`. Prune at the highest abundance of a series, since an element's lines grow with its abundance.
- Lines are read inside the windows widened by `dlam` and kept if their opacity on the run's wavelength grid passes HOP's negligibility test on the sampled depths, so `dlam` and `R_sampling` have to match the runs.
- Each record ends with the line's strength and its reach, the distance in Å out to which its extinction exceeds the tolerance times the continuum. A run reading a list whose lines carry a reach evaluates each line only where it reaches.
- The output is written to `<output_file>.<pid>.partial` and renamed when complete. `<output_file>.provenance.nml` records every input group the selection used and a `&prune_result` group with the sampled depths, the effective vmic and [Fe/H], the windows and the line count. Invalid input exits with status 1 and one line on stderr.

The pruned list contains every line the corresponding `dispatch.x` run keeps, and a few more: the program evaluates the model's own depth points, whose deepest layer is hotter than any cell of the run's interpolated mesh. In OpenMULTI3D, `tumult.prune` runs the program for the synthesis workflows.

## Development

HOP is developed in the OpenMULTI3D repository, of which this repository is a mirror of the `HOP/` folder.

```bash
make -C standalone -j8      # standalone/build/libhop.so (copied to hopacity/), hop_prune.x, atom_convert.x
make -C standalone clean    # rebuild from scratch after editing *.f90
```

Invoked here, `COMPILER=ifort` or `COMPILER=ifx` selects the Intel compilers and `OPT=-O2` overrides the optimisation (default `-O3`). In an OpenMULTI3D checkout the top-level `make` builds the library alongside `dispatch.x`, with the same compiler and option group, and `pip install -e .` there provides `hopacity` together with `tumult`. The package loads `hopacity/libhop.so`, or the library `$HOP_LIB` points to.

`standalone/libhop.h` declares the C entry points and is committed; `make -C standalone bindings` rewrites it from the `BIND(C)` declarations in `hop_c_api.f90` with `gen_bindings.py` after an interface change. Every `BIND(C)` procedure has to be public, or gfortran 16 does not export it. OpenMULTI3D's `tests/native/test_hop_abi.py` checks both.

`pip wheel --no-deps -w dist .` builds a wheel; `delocate-wheel` (macOS) or `auditwheel repair` (Linux) copies the compiler's runtimes into it. The release workflow (`.github/workflows/wheels.yml`) does this with cibuildwheel for every pushed `v*` tag, runs `tests/smoke.py` against each wheel and attaches the wheels to the release.

## License

hopacity and HOP are free software: you can redistribute them and/or modify them under the terms of the GNU General Public License version 3 as published by the Free Software Foundation (`LICENSE`). They are distributed in the hope that they will be useful, but without any warranty; without even the implied warranty of merchantability or fitness for a particular purpose. The copyright holders are listed in `AUTHORS`.
