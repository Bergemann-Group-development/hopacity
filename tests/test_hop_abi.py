"""Keep the generated libhop ABI in sync with the library and the Python wrapper.

``standalone/libhop.h`` is written by ``standalone/gen_bindings.py`` from the
``BIND(c, name=...)`` declarations in ``hop_c_api.f90`` and committed (``make -C
standalone bindings`` rewrites it after an interface change); the build does not
regenerate it, so these tests do, and fail on a diff.
The header is also what ``hopacity`` hands to ``ffi.cdef``, so every name and signature
the wrapper calls comes from here; the calls the other HOP tests make through cffi are
what checks the generator's type mapping.
"""
import importlib.util
import re
import subprocess
import sys

import pytest

from hop_reference import TESTS_DIR, hop_library

HOP_DIR = TESTS_DIR.parent
HOP_LIB = hop_library()
STANDALONE = HOP_DIR / "standalone"
HEADER = STANDALONE / "libhop.h"
SHIM_F90 = STANDALONE / "hop_c_api.f90"
# the Python wrapper, the only Python code that calls libhop
PY_HOSTS = sorted((HOP_DIR / "hopacity").glob("*.py"))

pytestmark = pytest.mark.skipif(
    HOP_LIB is None,
    reason="libhop.so not built - run 'make -C standalone -j8'")


def _header_names(text):
    return sorted(re.findall(r"\b(hop_[a-z0-9_]+)\s*\(", text))


def _regenerate_header():
    """Run gen_bindings.py on the shim, as the Makefile's bindings target does."""
    spec = importlib.util.spec_from_file_location("gen_bindings", STANDALONE / "gen_bindings.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.prototypes(SHIM_F90.read_text())


def test_header_matches_the_fortran():
    """The committed header is what the generator writes today."""
    assert HEADER.read_text() == _regenerate_header(), (
        "standalone/libhop.h is stale - regenerate with "
        "'make -C standalone bindings' and commit it")


def test_header_declares_every_bind_c_entry_point():
    """Header names == the BIND(c, name='...') names in the shim."""
    declared = sorted(re.findall(r"BIND\(c,\s*name='([a-z0-9_]+)'\)", SHIM_F90.read_text()))
    assert _header_names(HEADER.read_text()) == declared


def test_every_bind_c_entry_point_is_public():
    """A private BIND(C) procedure is exported by gfortran 15 but not by gfortran 16, so
    the export test above passes on one compiler and the wheel fails on the other."""
    src = SHIM_F90.read_text()
    declared = set(re.findall(r"BIND\(c,\s*name='([a-z0-9_]+)'\)", src))
    public = set()
    for stmt in re.findall(r"^\s*public\s*::((?:[^\n]*&[ \t]*\n)*[^\n]*)", src, re.M | re.I):
        public |= set(re.findall(r"\bhop_[a-z0-9_]+", stmt))
    assert not declared - public, f"not public in hop_c_api: {sorted(declared - public)}"


def test_header_matches_the_exported_symbols():
    """Header names == what the built library actually exports."""
    # a Linux wheel's library is stripped, so its exports are read from the dynamic table
    flags = ["-gU"] if sys.platform == "darwin" else ["-D", "--defined-only"]
    out = subprocess.run(["nm", *flags, str(HOP_LIB)], capture_output=True, text=True).stdout
    # macOS prefixes global symbols with an underscore; Linux does not
    exported = sorted({m for m in re.findall(r"^\S+ [TDS] _?(hop_[a-z0-9_]+)$", out, re.M)})
    assert _header_names(HEADER.read_text()) == exported


def _fortran_array_dummies():
    """Per entry point, the dummies declared with a single dimension.

    ``out(*)``, ``temp(n)``, ...  Multi-dimensional locals (``t3(:,:,:)``) and deferred
    shapes (``lam8(:)``) do not match, so what is left is the set of names that must
    cross the C boundary as pointers.  Kept per-function because the same name can be an
    array in one entry point and a scalar in another - ``value`` is ``value(*)`` in
    hop_set_option_str but a by-value scalar in hop_set_option_int.
    """
    out, cur = {}, None
    for line in SHIM_F90.read_text().splitlines():
        m = re.search(r"BIND\(c,\s*name='([a-z0-9_]+)'\)", line)
        if m:                                   # follows the FUNCTION line, before decls
            cur = m.group(1)
            out[cur] = set()
            continue
        if re.match(r"\s*END\s+FUNCTION", line, re.I):
            cur = None
            continue
        if cur is None or "::" not in line:
            continue
        lhs, rhs = line.split("::", 1)
        if not re.search(r"\b(real|integer|character|logical)\b", lhs, re.I):
            continue
        out[cur].update(re.findall(r"(\w+)\s*\(\s*(?:\*|\w+)\s*\)", rhs))
    return out


def test_no_array_argument_is_passed_by_value():
    """Every array dummy must reach the header as a pointer.

    gfortran renders an explicit-shape dummy (``temp(n)``) as a by-value scalar
    (``const double temp``) and only gets ``(*)`` right, so an explicit-shape declaration
    silently produces a header that passes arrays by value.  Deriving the expected set
    from the Fortran rather than whitelisting scalar names means a regression back to
    ``temp(n)`` still fails here.
    """
    arrays = _fortran_array_dummies()
    assert arrays.get("hop_load_state") >= {"temp", "rho", "ne"}, "decl scan broken?"
    for proto in HEADER.read_text().splitlines():
        fn = re.search(r"\b(hop_[a-z0-9_]+)\s*\(", proto).group(1)
        for param in proto[proto.index("(") + 1:proto.rindex(")")].split(","):
            param = param.strip()
            if not param:
                continue
            name = param.split()[-1].lstrip("*")
            if name in arrays.get(fn, ()) and "*" not in param:
                pytest.fail(f"{proto.strip()}: {name!r} is passed by value - declare it "
                            f"assumed-size '(*)' in hop_c_api.f90")


def test_python_wrapper_binds_only_exported_symbols():
    """Every hop_* name a hopacity module calls must exist in the library.

    Scans every module, since a module nothing imports can otherwise bind a name the
    library never exported.
    """
    exported = set(_header_names(HEADER.read_text()))
    used = {}
    for path in PY_HOSTS:
        used.setdefault(path, set()).update(
            re.findall(r"\b_?lib\.(hop_[a-z0-9_]+)", path.read_text()))

    bogus = {p.name: sorted(names - exported) for p, names in used.items()
             if names - exported}
    assert not bogus, f"these modules call libhop symbols that do not exist: {bogus}"


def test_option_names_match_the_select_case():
    """hop_option_names must list exactly what hop_set_option's select case handles.

    The two sit next to each other in hop_c_api.f90 and are maintained by hand, so this
    is what keeps an option from being settable but unlisted, or listed but rejected.
    """
    src = SHIM_F90.read_text()
    body = src[src.index("FUNCTION hop_set_option("):src.index("END FUNCTION hop_set_option")]
    handled = set(re.findall(r"case \('([a-z0-9_]+)'\)", body))
    listed = set(re.findall(r"'([a-z0-9_ ]+)'",
                            src[src.index("option_names = "):src.index("CONTAINS")]))
    listed = {w for chunk in listed for w in chunk.split()}
    assert handled == listed, (f"only in select case: {sorted(handled - listed)}; "
                               f"only in option_names: {sorted(listed - handled)}")


# Exported deliberately without an external caller.  Keep this list to one line per
# entry, with the reason - anything else added to the shim must be wrapped by a host.
INTERNAL_ONLY = {
    # released by hop_free already; exported so a host can hand the OpenMP threads
    # back without tearing down the loaded state.  No host needs that yet.
    "hop_omp_release",
}


def test_hosts_cover_every_entry_point():
    """Every exported entry point is reached from a host wrapper - no silent gaps.

    This is what keeps an entry point from being added to the library and then never
    wrapped.  There are three hosts: the Python wrapper (hopacity/), the Julia one
    (HOP.jl), which owns the threaded entry points HopEOS has no use for, and the
    hop_prune.x program next to the shim, which USEs the entry points it drives.
    INTERNAL_ONLY carries the documented exceptions.
    """
    exported = set(_header_names(HEADER.read_text()))
    wrapped = set(INTERNAL_ONLY)
    prune_src = (SHIM_F90.parent / "hop_prune.f90").read_text()
    wrapped |= set(re.findall(r"\b(hop_[a-z0-9_]+)\b",
                              prune_src[prune_src.index("USE hop_c_api"):
                                        prune_src.index("USE hop_env_mod")]))
    for path in PY_HOSTS:
        src = path.read_text()
        # attribute access (lib.hop_x), or a quoted name - how set_option picks its setter
        wrapped |= set(re.findall(r"\b\w*lib\.(hop_[a-z0-9_]+)", src))
        wrapped |= set(re.findall(r"['\"](hop_[a-z0-9_]+)['\"]", src))
    for path in (HOP_DIR / "HOP.jl" / "src").glob("*.jl"):
        # ccall targets are symbols: _fp(:hop_x), or :hop_x straight
        wrapped |= set(re.findall(r":(hop_[a-z0-9_]+)", path.read_text()))
    assert not exported - wrapped, f"no host caller for: {sorted(exported - wrapped)}"
