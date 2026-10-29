#!/usr/bin/env python3
"""Write the C prototypes of the BIND(C) functions in hop_c_api.f90.

The prototypes are what ``hopacity`` hands to ``cffi.cdef``; the Makefile's
``bindings`` target commits them as ``libhop.h``. The mapping is the C interoperability
of Fortran 2003 for the small vocabulary the shim uses: ``integer(c_int)``,
``real(c_float)``, ``real(c_double)`` and ``character(kind=c_char)`` dummies that are
either ``value`` scalars or assumed-size arrays, and ``c_int``/``c_double`` results.
``intent(in)`` arrays become ``const`` pointers. Anything else in a BIND(C) interface is
an error, so a new kind of argument is added here deliberately.

Usage::

    python3 gen_bindings.py hop_c_api.f90 > libhop.h
"""
from __future__ import annotations

from pathlib import Path
import re
import sys

C_TYPES = {"integer(c_int)": "int", "real(c_float)": "float", "real(c_double)": "double",
           "character(kind=c_char)": "char"}

HEAD = re.compile(r"^\s*function\s+(\w+)\s*\((.*?)\)\s*result\s*\((\w+)\)\s*"
                  r"bind\s*\(\s*c\s*,\s*name\s*=\s*'(\w+)'\s*\)", re.I)
DECL = re.compile(r"^\s*(integer|real|character)\s*\(([^)]*)\)\s*(.*?)\s*::\s*(.*)$", re.I)
END = re.compile(r"^\s*end\s+function\b", re.I)


def logical_lines(text):
    """Source lines with comments removed and ``&`` continuations joined."""
    lines, pending = [], ""
    for raw in text.splitlines():
        code = re.sub(r"!.*", "", raw) if "'" not in raw else _strip_comment(raw)
        code = code.strip()
        if pending:
            code = pending + code.lstrip("&")
            pending = ""
        if code.endswith("&"):
            pending = code[:-1]
            continue
        lines.append(code)
    return lines


def _strip_comment(raw):
    quoted = False
    for i, ch in enumerate(raw):
        if ch == "'":
            quoted = not quoted
        elif ch == "!" and not quoted:
            return raw[:i]
    return raw


def parse_declarations(lines):
    """``{name: (ctype, attributes, is_array)}`` from a function's declaration lines."""
    out = {}
    for line in lines:
        m = DECL.match(line)
        if not m:
            continue
        kind = re.sub(r"\s+", "", m.group(2)).lower()
        key = f"{m.group(1).lower()}({kind})"
        if key not in C_TYPES:
            continue                      # a local variable of a non-interoperable kind
        attrs = {a.strip().lower() for a in m.group(3).strip(", ").split(",") if a.strip()}
        for item in m.group(4).split(","):
            name, _, dims = item.strip().partition("(")
            out[name.lower()] = (C_TYPES[key], attrs, bool(dims))
    return out


def prototype(name, args, result, decls, c_name):
    parts = []
    for arg in args:
        try:
            ctype, attrs, is_array = decls[arg]
        except KeyError:
            raise ValueError(f"{c_name}: no interoperable declaration for argument {arg}") from None
        if "value" in attrs:
            if is_array:
                raise ValueError(f"{c_name}: {arg} is an array passed by value")
            parts.append(f"{ctype} {arg}")
        else:
            const = "const " if "intent(in)" in attrs else ""
            parts.append(f"{const}{ctype} *{arg}")
    try:
        rtype, rattrs, rarray = decls[result]
    except KeyError:
        raise ValueError(f"{c_name}: result {result} has no interoperable declaration") from None
    if rarray or rtype == "char":
        raise ValueError(f"{c_name}: result {result} is not a C scalar")
    return f"{rtype} {c_name} ({', '.join(parts)});"


def prototypes(source):
    lines = logical_lines(source)
    found = {}
    i = 0
    while i < len(lines):
        m = HEAD.match(lines[i])
        if not m:
            i += 1
            continue
        fname, arglist, result, c_name = m.groups()
        args = [a.strip().lower() for a in arglist.split(",") if a.strip()]
        body = []
        i += 1
        while i < len(lines) and not END.match(lines[i]):
            body.append(lines[i])
            i += 1
        found[c_name] = prototype(fname, args, result.lower(), parse_declarations(body), c_name)
    return "".join(found[k] + "\n" for k in sorted(found))


def main(argv):
    if len(argv) != 2:
        sys.exit(__doc__)
    sys.stdout.write(prototypes(Path(argv[1]).read_text()))


if __name__ == "__main__":
    main(sys.argv)
