module HOP

using Libdl
using Distributed
using Mmap
using FileWatching.Pidfile: mkpidlock

export HopEOS, load_state!, load_linelist!, opacity, opacity_threaded, opacity_threaded!, wavelength_grid, logstep, freq_split, close!, meshgrid, compute_grid

const _HOP_DIR = dirname(dirname(@__DIR__))   # …/HOP, up from HOP.jl/src/
const _DATA_DIR = joinpath(_HOP_DIR, "hopacity", "data")
const _LIB  = Ref{Ptr{Cvoid}}(C_NULL)
const _PATH = Ref{String}("")
const _SYMS = Dict{Symbol,Ptr{Cvoid}}()
const _PROC_LOCK = ReentrantLock()   # one live HOP instance per process

const ENERGY_OFFSET = (13.595 + 5.0) / 2.380491e-24 * 1.60218e-12
const HOP_TILE_ELEMS_AUTO = 30 * 25 * 50_000
const HOP_TILE_ELEMS_CEIL = 1_000_000_000

# ============================================================================
# Recipe: everything needed to rebuild an identical instance
# ============================================================================

struct HopRecipe
    libpath::String
    absdat_file::String
    abund_file::String
    isotope_file::String
    molec_file::String
    feh::Float64
    tmolim::Float64
    alpha_enhancement::Float64
    r_process_enhancement::Float64
    verbose::Int
    options::Dict{String,Any}   
end

# ============================================================================
# HopEOS
# ============================================================================

mutable struct HopEOS
    recipe::HopRecipe
    owns::Bool
    shape::Dims                                   # shape of the loaded state
    state::Union{Nothing,NamedTuple}              # (temp, rho, ne, vmic) replay args
    species::Vector{String}
    linelist::Union{Nothing,Tuple}                # (files, dlam, lam_min, lam_max)
end

"""
    HopEOS(; feh=0.0, absdat_file=…, abund_file=…, isotope_file=…, molec_file=…,
           tmolim=-99.0, alpha_enhancement=0.0, r_process_enhancement=0.0,
           verbose=false, libpath=find_library(), options...)

Allocate and initialise the HOP composition. Data-file paths default to the
standard files under `HOP/hopacity/data/`. `options` are module-level HOP
settings applied *before* `hop_init` (the `case` labels of
`hop_set_option_{str,real,int}`), e.g. `hydrogen_profiles="lemke"`,
`vcs_data_dir="…"`, `abundances="…"`, `disable_stk=true`.

One live instance per process; use `close!` (or let the finalizer run).
"""
function HopEOS(; 
                feh=0.0,
                absdat_file = joinpath(_DATA_DIR, "absdat", "TS_absdat.txt"),
                abund_file  = joinpath(_DATA_DIR, "abund", "abund_magg22.txt"),
                isotope_file = joinpath(_DATA_DIR, "isotopes", "isotopes_blp.txt"),
                molec_file   = joinpath(_DATA_DIR, "molec", "molec_data_BC.txt"),
                tmolim = -99.0, alpha_enhancement = 0.0, r_process_enhancement = 0.0,
                verbose = false, libpath = find_library(), options...)
    opts = Dict{String,Any}(String(k) => v for (k, v) in options)

    get!(opts, "vcs_data_dir", joinpath(_DATA_DIR, "vcs"))
    get!(opts, "anstee_file", joinpath(_DATA_DIR, "Anstee-1802-2014.txt"))
    for p in (absdat_file, abund_file, isotope_file, molec_file)
        isfile(p) || error("data file not found: $p")
    end

    recipe = HopRecipe(
        libpath, absdat_file, abund_file, isotope_file, molec_file,
        Float64(feh), Float64(tmolim), Float64(alpha_enhancement),
        Float64(r_process_enhancement), _verbosity(verbose), opts
    )

    eos = HopEOS(recipe, false, (), nothing, String[], nothing)

    _init_fortran!(eos)
    finalizer(close!, eos)
    eos
end

# ============================================================================
# Library handle
# ============================================================================

standalone_dir() = joinpath(_HOP_DIR, "standalone")

# Own build dir, so the Makefile's -O1 default (what other users get) and these
# flags never share objects.
const BUILD_DIR = "build-opt"
const BUILD_OPT = "-O3 -march=native -mprefer-vector-width=256"

function __init__()
    dir = standalone_dir()
    isdir(dir) || return nothing
    mkpidlock(joinpath(dir, "build.lock"); stale_age=600) do
        run(`make -s -C $dir -j$(min(Sys.CPU_THREADS, 8)) BUILD=$BUILD_DIR OPT=$BUILD_OPT`)
    end
    nothing
end

function find_library()
    lib = get(ENV, "HOP_LIB", "")
    path = isempty(lib) ? joinpath(standalone_dir(), BUILD_DIR, "libhop.so") : lib
    isfile(path) || error("libhop.so not found at $path; build it with " *
                          "`make -C $(standalone_dir()) BUILD=$BUILD_DIR OPT=\"$BUILD_OPT\"` " *
                          "or set \$HOP_LIB")
    path
end

_ensure_lib!(path::String=find_library()) = begin
    if _LIB[] == C_NULL || _PATH[] != path
        _LIB[] = Libdl.dlopen(path, Libdl.RTLD_GLOBAL | Libdl.RTLD_NOW)
        _PATH[] = path
        empty!(_SYMS)
    end
    _LIB[]
end

_fp(name::Symbol) = get!(() -> Libdl.dlsym(_LIB[], name), _SYMS, name)

# ============================================================================
# ccall wrappers
# ============================================================================

_hop_init(absdat, abund, iso, molec, feh, tmolim, alpha, rproc) = ccall(
    _fp(:hop_init), Cint,
    (Cstring, Cstring, Cstring, Cstring, Cdouble, Cfloat, Cfloat, Cfloat),
    absdat, abund, iso, molec, feh, tmolim, alpha, rproc
)

_hop_free() = ccall(_fp(:hop_free), Cint, ())

_hop_set_option(name, val) = ccall(_fp(:hop_set_option), Cint, (Cstring, Cstring), name, val)
_hop_last_error(buf, buflen) = ccall(_fp(:hop_last_error), Cint, (Ptr{UInt8}, Cint), buf, buflen)

_hop_load_state(n, temp, ne, rho, vmic, use_rho, use_ne) = ccall(
    _fp(:hop_load_state), Cint,
    (Cint, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Cfloat, Cint, Cint),
    n, temp, ne, rho, vmic, use_rho, use_ne
)

_hop_get_n() = ccall(_fp(:hop_get_n), Cint, ())
_hop_get_field(name, n, out) = ccall(_fp(:hop_get_field), Cint, (Cstring, Cint, Ptr{Float64}), name, n, out)
_hop_get_grph() = ccall(_fp(:hop_get_grph), Cdouble, ())
_hop_keep_species(id) = ccall(_fp(:hop_keep_species), Cint, (Cstring,), id)
_hop_set_vmic(n, vmic) = ccall(_fp(:hop_set_vmic), Cint, (Cint, Cfloat), n, vmic)
_hop_init_linelist(file, dlam, nlam, lam) = ccall(
    _fp(:hop_init_linelist), Cint, (Cstring, Cdouble, Cint, Ptr{Float64}), file, dlam, nlam, lam
)

_hop_opacity(nlam, lam, chi, scat, absmet, llist, dlam, lscat) = ccall(
    _fp(:hop_opacity), Cint,
    (Cint, Ptr{Float64}, Ptr{Float32}, Ptr{Float32}, Cint, Cint, Cdouble, Cint),
    nlam, lam, chi, scat, absmet, llist, dlam, lscat
)

# --- threaded opacity: OpenMP over cells, inside Fortran ---
_hop_opacity_omp(nlam, lam, chi, scat, absmet, llist, dlam, lscat, nthreads) = ccall(
    _fp(:hop_opacity_omp), Cint,
    (Cint, Ptr{Float64}, Ptr{Float32}, Ptr{Float32}, Cint, Cint, Cdouble, Cint, Cint),
    nlam, lam, chi, scat, absmet, llist, dlam, lscat, nthreads)

_verbosity(v::Bool) = v ? 1 : -1
_verbosity(v) = Int(v)

# ============================================================================
# HopEOS utility functions
# ============================================================================

_init_fortran!(eos::HopEOS) = begin
    _ensure_lib!(eos.recipe.libpath)
    lock(_PROC_LOCK)
    try
        set_option("verbose", eos.recipe.verbose)          # before hop_init
        for (name, value) in eos.recipe.options
            set_option(name, value)
        end
        status = _hop_init(eos.recipe.absdat_file, eos.recipe.abund_file,
                           eos.recipe.isotope_file, eos.recipe.molec_file,
                           eos.recipe.feh, eos.recipe.tmolim,
                           eos.recipe.alpha_enhancement, eos.recipe.r_process_enhancement)
        status == 0 || error("hop_init failed with status $status " *
                             "(1 = another instance is already open in this process)")
        eos.owns = true
    finally
        unlock(_PROC_LOCK)
    end
    eos
end

_last_error() = begin
    buf = Vector{UInt8}(undef, 512)
    status = _hop_last_error(buf, Cint(length(buf)))
    status == 0 || return ""
    i = something(findfirst(==(0x00), buf), length(buf) + 1) - 1
    String(buf[1:i])
end

function set_option(name::AbstractString, value)
    key = lowercase(String(name))
    valstr = value isa AbstractString ? String(value) :
             value isa Bool           ? (value ? "true" : "false") : string(value)
    status = _hop_set_option(key, valstr)
    status == 0 && return
    msg = _last_error()
    error("hop_set_option($name = $valstr) failed with status $status" *
          (isempty(msg) ? "" : ": $msg"))
end

function load_state!(eos::HopEOS; temp, rho=nothing, ne=nothing, vmic=-1.0, species=nothing)
    species !== nothing && keep_species!(eos, species)
    ta = Float64.(collect(temp))
    shape = size(ta)
    use_rho = rho !== nothing
    use_ne  = ne  !== nothing
    (use_rho || use_ne) || error("load_state! needs rho and/or ne")
    ra = use_rho ? Float64.(collect(rho)) : zeros(Float64, shape)
    na = use_ne  ? Float64.(collect(ne))  : zeros(Float64, shape)
    (size(ra) == shape && size(na) == shape) || error("temp, rho, ne must share shape")
    (all(isfinite, ta) && all(>(0), ta)) || error("temp must be finite and positive")

    lock(_PROC_LOCK)
    try
        status = _hop_load_state(prod(shape), vec(ta), vec(na), vec(ra), Cfloat(vmic), Cint(use_rho), Cint(use_ne))
        status == 0 || error("hop_load_state failed with status $status")
    finally
        unlock(_PROC_LOCK)
    end
    eos.shape = shape
    eos.state = (temp = ta, rho = use_rho ? ra : nothing, ne = use_ne ? na : nothing, vmic = Float64(vmic))
    eos
end

function keep_species!(eos::HopEOS, ids)
    for sid in (ids isa AbstractString ? (ids,) : ids)
        status = _hop_keep_species(String(sid))
        status == 0 || error("unknown species id $sid (use e.g. \"Fe I\", \"CH\")")
        sid in eos.species || push!(eos.species, String(sid))
    end
    eos
end

function close!(eos::HopEOS)
    # only the instance that owns the Fortran state may free it
    if eos.owns
        eos.owns = false
        lock(_PROC_LOCK)
        try
            _hop_free()
        finally
            unlock(_PROC_LOCK)
        end
        eos.state = nothing
        eos.shape = ()
    end
    nothing
end

# ============================================================================
# Extract outputs
# ============================================================================

_field(eos::HopEOS, name::AbstractString) = begin
    nn = n(eos)
    nn == 0 && error("no state loaded; call load_state! first")
    out = Vector{Float64}(undef, nn)
    status = _hop_get_field(name, nn, out)
    status == 0 || error("hop_get_field($name) failed with status $status")
    reshape(out, eos.shape)
end

temp(eos::HopEOS)   = _field(eos, "temp")
ne(eos::HopEOS)     = _field(eos, "ne")
rho(eos::HopEOS)    = _field(eos, "rho")
pe(eos::HopEOS)     = _field(eos, "pe")
pg(eos::HopEOS)     = _field(eos, "pg")
energy(eos::HopEOS) = _field(eos, "e")
eh(eos::HopEOS)     = _field(eos, "eh")
chi500(eos::HopEOS) = _field(eos, "chi500")
grph(eos::HopEOS)   = _hop_get_grph()

n(eos::HopEOS) = Int(_hop_get_n())

# ============================================================================
# Line lists
# ============================================================================

_linelist_files(files) = begin
    if files isa AbstractString
        p = String(files)
        fs = isdir(p) ? sort!(filter(isfile, joinpath.(p, readdir(p)))) : [p]
    else
        fs = String.(collect(files))
    end
    for f in fs
        isfile(f) || error("line-list file not found: $f")
    end
    Tuple(fs)
end

_load_linelist!(eos::HopEOS, files::Tuple, dlam, lam) = begin
    lam_min, lam_max = Float64(minimum(lam)), Float64(maximum(lam))
    la = Float64[lam_min, lam_max]
    for f in files
        status = _hop_init_linelist(f, Float64(dlam), Cint(2), la)
        status == 2 && error("hop_init_linelist needs an initialised instance")
        status == 0 || error("hop_init_linelist($f) failed with status $status")
    end
    eos.linelist = (files, Float64(dlam), lam_min, lam_max)
    eos
end

"""
    load_linelist!(eos, linelist; dlam=100.0, lam)

Load a line list (file, folder, or iterable) into `eos`.
Only lines within `[minimum(lam)-dlam, maximum(lam)+dlam]` are loaded.
Re-solves the state if one is already loaded. Calling with a different list
or wavelength range requires a fresh `HopEOS`.
"""
function load_linelist!(eos::HopEOS, linelist; dlam=100.0, lam)
    eos.owns || error("this HopEOS does not own an initialised instance")
    files = _linelist_files(linelist)
    lamv  = _check_lam(lam)
    spec  = (files, Float64(dlam), Float64(minimum(lamv)), Float64(maximum(lamv)))
    if eos.linelist !== nothing
        eos.linelist == spec && return eos
        error("a line list is already loaded; HOP appends and cannot clear it — " *
              "build a fresh HopEOS for a different list or wavelength range")
    end
    lock(_PROC_LOCK)
    try
        _load_linelist!(eos, files, dlam, lamv)
    finally
        unlock(_PROC_LOCK)
    end
    # mirror hop.py: if a state was already solved, re-solve so comp sees the new lines
    if eos.state !== nothing
        s = eos.state
        load_state!(eos; temp=s.temp, rho=s.rho, ne=s.ne, vmic=s.vmic)
    end
    eos
end

_set_vmic!(eos::HopEOS, vmic) = (_hop_set_vmic(n(eos), Cfloat(vmic)); nothing)

_opacity_call(lam::Vector{Float64}, include_absmet, include_llist, dlam, line_scattering) = begin
    nn = Int(_hop_get_n())
    chi  = Matrix{Float32}(undef, length(lam), nn)
    scat = Matrix{Float32}(undef, length(lam), nn)
    status = _hop_opacity(length(lam), lam, chi, scat, Cint(include_absmet),
                          Cint(include_llist), Float64(dlam), Cint(line_scattering))
    status == 0 || error("hop_opacity failed with status $status")
    chi, scat
end

# ============================================================================
# Compute opacities
# ============================================================================

_statelen(s::NamedTuple)   = length(s.temp)
_stateshape(s::NamedTuple) = size(s.temp)

_reshape_out(chi::AbstractMatrix, shape::Dims) = begin
    nlam = size(chi, 1)
    nd = length(shape)
    a = reshape(chi, (nlam, shape...))
    permutedims(a, (ntuple(i -> i + 1, nd)..., 1))
end

_cellchunk(n, w, nc) = begin
    nv, md = divrem(n, nc)
    ((w - 1) * nv + min(md, w - 1), w * nv + min(md, w))
end

"""
    opacity(eos, lam; linelist=nothing, dlam=100.0, include_absmet=false,
            line_scattering=false, vmic=nothing)

Compute extinction `chi` and scattering `scat` [1/cm] for the given vacuum wavelengths.
Returns two `Float32` arrays shaped `<state shape> × nλ`. Computes serially.
Passing `linelist` automatically loads it and re-solves the state.
"""
function opacity(eos::HopEOS, lam;
                 linelist=nothing, dlam=100.0, include_absmet=false,
                 line_scattering=false, vmic=nothing)
    lamv = _check_lam(lam)
    eos.state === nothing && error("no state loaded; call load_state! first")
    eos.owns || error("this HopEOS does not own an initialised instance")

    # mirror hop.py: passing a list here loads it and re-solves so comp sees the lines
    linelist !== nothing && load_linelist!(eos, linelist; dlam=dlam, lam=lamv)

    lock(_PROC_LOCK)
    try
        vmic !== nothing && _set_vmic!(eos, vmic)
        include_llist = eos.linelist !== nothing

        chi, scat = _opacity_call(lamv, include_absmet, include_llist, dlam, line_scattering)
        shape = _stateshape(eos.state)
        return _reshape_out(chi, shape), _reshape_out(scat, shape)
    finally
        unlock(_PROC_LOCK)
    end
end

"""
    opacity_threaded(eos, lam; dlam=100.0, line_scattering=false, vmic=nothing,
                     nthreads=nothing, include_llist=nothing)

Compute extinction and scattering using Julia threads over the loaded state's T-ρ cells.
Returns newly allocated arrays shaped `<state shape> × nλ`. Does not support `include_absmet`.
"""
function opacity_threaded(eos::HopEOS, lam;
                          dlam=100.0, line_scattering=false, vmic=nothing,
                          nthreads=nothing, include_llist=nothing)
    lamv  = _check_lam(lam)
    eos.state === nothing && error("no state loaded; call load_state! first")
    nlam  = length(lamv)
    ncell = Int(_hop_get_n())
    chi  = Matrix{Float32}(undef, nlam, ncell)
    scat = Matrix{Float32}(undef, nlam, ncell)
    opacity_threaded!(chi, scat, eos, lamv; dlam=dlam, line_scattering=line_scattering,
                      vmic=vmic, nthreads=nthreads, include_llist=include_llist)
    shape = _stateshape(eos.state)
    return _reshape_out(chi, shape), _reshape_out(scat, shape)
end

"""
    opacity_threaded!(chi, scat, eos, lam; dlam=100.0, line_scattering=false,
                      vmic=nothing, nthreads=nothing, include_llist=nothing)

In-place variant of `opacity_threaded`. Fills the provided `(nλ, ncell)` `Float32` 
matrices in column-major order with zero allocation.
"""
function opacity_threaded!(chi::AbstractMatrix{Float32}, scat::AbstractMatrix{Float32},
                           eos::HopEOS, lam;
                           dlam=100.0, line_scattering=false, vmic=nothing,
                           nthreads=nothing, include_llist=nothing)
    lamv = _check_lam(lam)
    eos.state === nothing && error("no state loaded; call load_state! first")
    eos.owns || error("this HopEOS does not own an initialised instance")
    nlam  = length(lamv)
    ncell = Int(_hop_get_n())
    size(chi)  == (nlam, ncell) || error("chi buffer must be ($nlam, $ncell), got $(size(chi))")
    size(scat) == (nlam, ncell) || error("scat buffer must be ($nlam, $ncell), got $(size(scat))")
    ilist = something(include_llist, eos.linelist !== nothing)
    lscat = line_scattering ? 1 : 0

    lock(_PROC_LOCK)
    try
        vmic !== nothing && _set_vmic!(eos, vmic)
        nthreads = clamp(something(nthreads, 1), 1, ncell)

        st = _hop_opacity_omp(nlam, lamv, chi, scat, Cint(0), Cint(ilist), Float64(dlam), Cint(lscat), Cint(nthreads))
        st == 0 || error("hop_opacity_omp failed with status $st: $(_last_error())")
    finally
        unlock(_PROC_LOCK)
    end
    return nothing
end

# ============================================================================
# Grid axis utility functions
# ============================================================================

function logstep(start, stop, step)
    lmin, lmax = log(start), log(stop)
    m = ceil(Int, (lmax - lmin) / step)
    exp.(lmin .+ step .* (0:m-1))
end

function meshgrid(ax...)
    grids = []
    space = Iterators.product(ax...)
    for i in eachindex(ax)
        append!(grids, [getindex.(space, i)])
    end

    grids
end

function freq_split(nnu::Integer, i::Integer, n_windows::Integer)
    nv, md = divrem(nnu, n_windows)
    j = mod(abs(i) - 1, n_windows) + 1
    # Return 1-based closed interval (start:end) instead of 0-based half-open
    ((j - 1) * nv + min(md, j - 1) + 1, j * nv + min(md, j))
end

function wavelength_grid(aa_blue, aa_red; R=0, daa=0, nnu=0, in_air=true)
    if R != 0
        blue, red = in_air ? (air2vac(aa_blue), air2vac(aa_red)) : (Float64(aa_blue), Float64(aa_red))
        fac = 1.0 + 1.0 / R
        m = ceil(Int, log(red / blue) / log(fac)) + 1
        lam = fill(fac, m); lam[1] = blue
        return accumulate(*, lam)
    end
    if daa != 0
        nnu = round(Int, (aa_red - aa_blue) / daa) + 1
    elseif nnu != 0
        daa = (aa_red - aa_blue) / (nnu - 1)
    else
        error("give R, daa or nnu together with aa_blue/aa_red")
    end
    lam = aa_blue .+ daa .* (0:nnu-1)
    in_air ? air2vac.(lam) : Float64.(lam)
end

# ============================================================================
# Wavelength checks
# ============================================================================

_n_air(air::Float64) = begin
    s2 = (1e4 / air)^2
    1.00008336624212083 + 0.02408926869968 / (130.1065924522 - s2) +
        0.0001599740894897 / (38.92568793293 - s2)
end

_check_lam(lam) = begin
    v = Float64.(collect(vec(lam)))
    (ndims(lam) <= 1 && length(v) >= 1) || error("lam must be a 1D array")
    (length(v) == 1 || all(>(0), diff(v))) || error("lam must be strictly ascending (vacuum Å)")
    v
end

function air2vac(air::Real)
    a = Float64(air)
    w = a * _n_air(a)
    w = a * _n_air(w)
    a * _n_air(w)
end

# ============================================================================
# Multiprocessing setup
# ============================================================================

function setup_workers!(needed)
    hop_proj = joinpath(@__DIR__, "..")
    spawned_workers = Distributed.addprocs(needed; exeflags="--project=$(hop_proj)")
    
    @eval Distributed.@everywhere $spawned_workers begin
        _old_out = stdout
        _old_err = stderr
        redirect_stdout(devnull)
        redirect_stderr(devnull)
    end
    
    @eval Distributed.@everywhere $spawned_workers using HOP

    @eval Distributed.@everywhere $spawned_workers begin
        redirect_stdout(_old_out)
        redirect_stderr(_old_err)
    end

    spawned_workers
end

function remove_workers!(spawned_workers)
	if !isempty(spawned_workers)
        Distributed.rmprocs(spawned_workers)
    end
end

# ============================================================================
# Grid computation using tiling + wavelength splits + OMP
# ============================================================================

function _hop_fill_block!(eos, t, rho, i0, i1, j0, j1, lam,
                          chi, scat, E2d, pg2d, ne2d, chi500_2d,
                          chiM, scatM; opacity_kwargs...)
    ti, rj   = t[i0:i1], rho[j0:j1]
    ntb, nrb = length(ti), length(rj)
    ncell    = ntb * nrb
    nlam     = length(lam)
    tt = repeat(collect(ti),  1, nrb)          # (ntb, nρ_b), tt[i,j] = t[i]
    rr = repeat(collect(rj)', ntb, 1)          # (ntb, nρ_b), rr[i,j] = ρ[j]

    load_state!(eos; temp=tt, rho=rr)
    cM = @view chiM[:,  1:ncell]               # reused (nλ, ncell) buffers, filled in place
    sM = @view scatM[:, 1:ncell]
    
    #t_op = @elapsed begin
    opacity_threaded!(cM, sM, eos, lam; opacity_kwargs...)   # threads over cells
    #end
    #println("Worker $(Distributed.myid()): opacity_threaded! (OpenMP region) took $(round(t_op, digits=2)) seconds.")

    rho_blk = HOP.rho(eos)
    E2d[i0:i1, j0:j1]       .= HOP.energy(eos)
    pg2d[i0:i1, j0:j1]      .= HOP.pg(eos)
    ne2d[i0:i1, j0:j1]      .= HOP.ne(eos)
    chi500_2d[i0:i1, j0:j1] .= HOP.chi500(eos) ./ rho_blk

    # (λ, cell) -> (T, ρ, λ) transpose directly into the cube, then [1/cm] -> [cm²/g]
    dchi = @view chi[i0:i1,  j0:j1, :]
    dsca = @view scat[i0:i1, j0:j1, :]
    permutedims!(dchi, reshape(cM, nlam, ntb, nrb), (2, 3, 1))
    permutedims!(dsca, reshape(sM, nlam, ntb, nrb), (2, 3, 1))
    @inbounds @views for li in 1:nlam
        dchi[:, :, li] ./= rho_blk
        dsca[:, :, li] ./= rho_blk
    end
    nothing
end

function _compute_grid_chunk(t, rho, lam; HOP_kwargs, opacity_kwargs,
                             t_block, rho_block, mmap, Teltype, verbose)
    ntemp, nrho, nlam = length(t), length(rho), length(lam)

    if t_block <= 0 || rho_block <= 0
        auto_cells = max(1, HOP_TILE_ELEMS_AUTO ÷ nlam)
        t_block    = t_block   > 0 ? t_block   : ntemp
        rho_block  = rho_block > 0 ? rho_block : clamp(auto_cells ÷ min(t_block, ntemp), 1, nrho)
    end

    t_block   = clamp(t_block,   1, ntemp)
    rho_block = clamp(rho_block, 1, nrho)
    ceil_cells = max(1, HOP_TILE_ELEMS_CEIL ÷ nlam)

    if t_block * rho_block > ceil_cells
        rho_block = clamp(ceil_cells ÷ t_block, 1, rho_block)
        t_block   = clamp(ceil_cells ÷ rho_block, 1, t_block)
        @warn "HOP tiling: tile exceeds the 32-bit-safe size for nλ=$nlam; reduced to " *
              "t_block=$t_block, rho_block=$rho_block (cells/tile=$(t_block*rho_block))."
    end

    chi, scat = if mmap
        opa_file  = tempname(pwd(), suffix=".bin")
        scat_file = tempname(pwd(), suffix=".bin")
        Mmap.mmap(open(opa_file,  "w+"), Array{Teltype, 3}, (ntemp, nrho, nlam)),
        Mmap.mmap(open(scat_file, "w+"), Array{Teltype, 3}, (ntemp, nrho, nlam))
    else
        Array{Teltype, 3}(undef, ntemp, nrho, nlam),
        Array{Teltype, 3}(undef, ntemp, nrho, nlam)
    end

    E2d  = Array{Float64}(undef, ntemp, nrho); pg2d      = similar(E2d)
    ne2d = similar(E2d);                       chi500_2d = similar(E2d)

    tkw = (; dlam=opacity_kwargs.dlam, line_scattering=opacity_kwargs.line_scattering,
             vmic=opacity_kwargs.vmic, nthreads=opacity_kwargs.nthreads)

    maxtile = t_block * rho_block
    chiM  = Matrix{Float32}(undef, nlam, maxtile)
    scatM = Matrix{Float32}(undef, nlam, maxtile)

    eos = HopEOS(; HOP_kwargs...)
    try
        #t_load = @elapsed begin
        isnothing(opacity_kwargs.linelist) || load_linelist!(eos, opacity_kwargs.linelist; dlam=opacity_kwargs.dlam, lam=lam)
        #end
        #println("Worker $(Distributed.myid()): load_linelist! took $(round(t_load, digits=2)) seconds.")
        
        for i0 in 1:t_block:ntemp
            i1 = min(i0 + t_block - 1, ntemp)
            for j0 in 1:rho_block:nrho
                j1 = min(j0 + rho_block - 1, nrho)
                _hop_fill_block!(eos, t, rho, i0, i1, j0, j1, lam,
                                 chi, scat, E2d, pg2d, ne2d, chi500_2d,
                                 chiM, scatM; tkw...)
                verbose && @info "HOP tile T[$i0:$i1] ρ[$j0:$j1] of ($ntemp, $nrho) done"
            end
        end
    finally
        close!(eos)
    end
    (chi, scat, E2d, pg2d, ne2d, chi500_2d)
end

function _compute_grid_tiled(
    T, rho, lam;
    HOP_kwargs=Dict{Symbol,Any}(:feh => 0.0),
    mmap=true,
    linelist=nothing,
    dlam=100.0, include_absmet=false, line_scattering=true, vmic=nothing,
    nthreads=nothing, nprocs=nothing,
    t_block=length(T), rho_block=length(rho),
    energy_offset=ENERGY_OFFSET,
    Teltype=Float32,
    verbose=false,
)
    opacity_kwargs = (; 
        linelist=linelist, dlam=dlam, include_absmet=include_absmet,line_scattering=line_scattering, vmic=vmic, nthreads=nthreads
    )

    np = isnothing(nprocs) ? 1 : max(1, nprocs)
    spawned_workers = Int[]

    try
        if np > 1 && length(lam) >= np
            verbose && @info "compute_eos: Running 2D parallelism ($np wavelength processes)"
            
            # Always spawn exactly np workers for this execution
            spawned_workers = setup_workers!(np)
            pool = Distributed.WorkerPool(spawned_workers)

            jobs = [(k, freq_split(length(lam), k, np)...) for k in 1:np]

            results = Distributed.pmap(pool, jobs) do (k, fs, fe)
                lam_k = lam[fs:fe]
                _compute_grid_chunk(T, rho, lam_k;
                    HOP_kwargs=HOP_kwargs, opacity_kwargs=opacity_kwargs,
                    t_block=t_block, rho_block=rho_block, mmap=false, Teltype=Teltype, verbose=false)
            end

            E2d_res, pg2d_res, ne2d_res, chi500_res = results[1][3], results[1][4], results[1][5], results[1][6]

            chi_full  = Array{Teltype, 3}(undef, length(T), length(rho), length(lam))
            scat_full = Array{Teltype, 3}(undef, length(T), length(rho), length(lam))

            for (k, fs, fe) in jobs
                chi_full[:, :, fs:fe]  .= results[k][1]
                scat_full[:, :, fs:fe] .= results[k][2]
            end

            E2d_res .+= energy_offset
            return (chi_full, scat_full, E2d_res, pg2d_res, ne2d_res, chi500_res)
        else
            res = _compute_grid_chunk(T, rho, lam;
                HOP_kwargs=HOP_kwargs, opacity_kwargs=opacity_kwargs,
                t_block=t_block, rho_block=rho_block, mmap=mmap, Teltype=Teltype, verbose=verbose)
            res[3] .+= energy_offset
            return res
        end
    finally
        remove_workers!(spawned_workers)
    end
end

"""
    compute_eos(T, rho, lam; nprocs=nothing, ...)

Compute the full EoS grid using 2D parallelism:
  - Wavelengths are split across `nprocs` multiprocessing workers.
  - T-ρ grid is tiled sequentially to minimize memory usage.
  - T-ρ points within each tile are threaded locally.
"""
compute_eos(args...; kwargs...) = _compute_grid_tiled(args...; kwargs...)

end # module
