#!/usr/bin/env python3
"""
Generate Fortran module dependency information for Make.

Given a list of .F90 source files, scan for:
- provided modules ("module <name>")
- used modules ("use [,...] :: <name>" or "use <name>")

Output make-style variable definitions and dependency rules:
  OBJS = path/to/file1.o path/to/file2.o ...
  MODS = module1.mod module2.mod ...
  path/to/target.o: path/to/provider.o
so that object files depending on modules build after the objects that provide them.

Notes:
- Ignores "module procedure", "end module", and "submodule" statements.
- Basic handling of line continuations with '&'.
- Case-insensitive parsing; emitted file paths retain original case.
"""
from __future__ import annotations
import os
import re
import sys
from typing import Dict, List, Set, Tuple

USE_RE_1 = re.compile(r"^\s*use\s*[,\s][^:]*::\s*([a-zA-Z0-9_]+)", re.IGNORECASE)
USE_RE_2 = re.compile(r"^\s*use\s+([a-zA-Z0-9_]+)", re.IGNORECASE)
MOD_RE = re.compile(r"^\s*module\s+([a-zA-Z0-9_]+)\b", re.IGNORECASE)
SUBMOD_RE = re.compile(r"^\s*submodule\b", re.IGNORECASE)
ENDMOD_RE = re.compile(r"^\s*end\s*module\b", re.IGNORECASE)
MODPROC_RE = re.compile(r"^\s*module\s+procedure\b", re.IGNORECASE)


def read_statements(path: str) -> List[str]:
    stmts: List[str] = []
    buf: str = ""
    try:
        with open(path, "r", encoding="utf-8", errors="ignore") as f:
            for raw in f:
                # strip comments
                line = re.sub(r"!.*", "", raw)
                if not line.strip():
                    continue
                # handle continuations
                s = line.rstrip().strip()
                if buf:
                    s = (buf + " " + s).strip()
                if s.endswith("&"):
                    buf = s[:-1].strip()
                    continue
                if s.startswith("&"):
                    # continuation from previous (unlikely after above), keep
                    s = s[1:].strip()
                    if buf:
                        s = (buf + " " + s).strip()
                stmts.append(s)
                buf = ""
        if buf:
            stmts.append(buf)
    except FileNotFoundError:
        pass
    return stmts


def to_obj(path: str) -> str:
    return os.path.splitext(path)[0] + ".o"


def to_src(path: str) -> str:
    """Normalize a path that may end in .o or .F90 to its .F90 source form."""
    return os.path.splitext(path)[0] + ".F90"


def discover_modules_and_uses(srcs: List[str]):
    provides: Dict[str, Set[str]] = {}
    uses_by_src: Dict[str, Set[str]] = {}

    for src in srcs:
        stmts = read_statements(src)
        uses: Set[str] = set()
        for st in stmts:
            st_l = st.lower()
            if SUBMOD_RE.match(st_l):
                continue
            if ENDMOD_RE.match(st_l):
                continue
            if MODPROC_RE.match(st_l):
                continue
            # provided module
            m = MOD_RE.match(st_l)
            if m:
                modname = m.group(1).lower()
                provides.setdefault(modname, set()).add(src)
                continue
            # use statements
            m1 = USE_RE_1.match(st_l)
            if m1:
                uses.add(m1.group(1).lower())
                continue
            m2 = USE_RE_2.match(st_l)
            if m2:
                uses.add(m2.group(1).lower())
                continue
        if uses:
            uses_by_src[src] = uses
    return provides, uses_by_src


def fmt_make_var(name: str, values: List[str], indent: str = "\t") -> str:
    """Format a Make variable assignment with continuation lines."""
    if not values:
        return f"{name} ="
    lines = [f"{name} = {values[0]}"]
    for v in values[1:]:
        lines[-1] += " \\"
        lines.append(f"{indent}{v}")
    return "\n".join(lines)


def main(argv: List[str]) -> int:
    if len(argv) < 2:
        print("Usage: gen_deps.py <file1.F90> [file2.F90 ...]", file=sys.stderr)
        return 1
    # Accept both .F90 and .o paths; normalize to .F90 for scanning
    srcs = list(dict.fromkeys(to_src(a) for a in argv[1:]))
    provides, uses_by_src = discover_modules_and_uses(srcs)

    # OBJS: object files corresponding to the given sources
    objs = [to_obj(s) for s in srcs]
    print(fmt_make_var("OBJS", objs))
    print()

    # MODS: .mod files for all modules defined in the given sources
    mods = sorted(f"{m}.mod" for m in provides)
    print(fmt_make_var("MODS", mods))
    print()

    # Dependency rules
    emitted: Set[Tuple[str, str]] = set()

    for src, used_mods in uses_by_src.items():
        for mod in sorted(used_mods):
            provs = provides.get(mod)
            if not provs:
                # likely an intrinsic or external library module; ignore
                continue
            for prov_src in provs:
                if prov_src == src:
                    continue
                dep = (to_obj(src), to_obj(prov_src))
                if dep not in emitted:
                    print(f"{dep[0]}: {dep[1]}")
                    emitted.add(dep)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
