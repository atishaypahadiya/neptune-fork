"""A6d leakage, static check (docs/PHASE1_SPEC.md, M6 acceptance tests): no code path in
M6a–M6d reads Neptune's real ephemeris or the M6e precovery file.

Pure AST, nothing imported or run. The runtime half is
``tests/test_m6c.py::test_dataset_is_table77_with_the_frozen_error_model_and_never_reads_neptune``.

1. The modules: everything under ``fork`` reachable by import (top-level or inside a function)
   from the M6a–M6d entry modules.
2. The fit path: a call graph over those modules, from the M6a transcription and loading, the
   M6b calibration and the M6c fit's ``dataset``, ``job_space``, ``do_job`` and ``run_all``.
   Calls are resolved by bare name to every function, method or class of that name in any of the
   modules, so the graph over-approximates the real one: a name clash can only add functions.
3. Neptune: a string constant containing "neptune" (docstrings aside), the ``NEPTUNE`` body or
   ``true_neptune_geocentric_longitude``. Every Horizons read takes a body name and looks it up
   in ``bodies.BODIES``, so a read of Neptune needs one of these in the path, or an iteration
   over ``BODIES`` (only ``horizons.default_fetch_plan`` does that, outside the fit path).
4. M6d's comparison constants (``JPL_*`` in ``fork.export.m6c``) are applied after sampling, to
   the posterior; the check also asserts that no fit-path function names one.
"""

from __future__ import annotations

import ast
from pathlib import Path

SRC = Path(__file__).resolve().parents[1] / "src"
ENTRY_MODULES = ("fork.data.transcription", "fork.data.observations", "fork.data.calibration",
                 "fork.export.m6c")  # fmt: skip
M6C_FIT = ("dataset", "job_space", "do_job", "run_all")
NEPTUNE_NAMES = {"NEPTUNE", "true_neptune_geocentric_longitude"}


def _path(mod: str) -> Path | None:
    p = SRC.joinpath(*mod.split("."))
    for q in (p.with_suffix(".py"), p / "__init__.py"):
        if q.exists():
            return q
    return None


def _tree(mod: str) -> ast.Module:
    return ast.parse(_path(mod).read_text(encoding="utf-8"))


def _imports(mod: str) -> set[str]:
    out = set()
    for n in ast.walk(_tree(mod)):
        if isinstance(n, ast.ImportFrom) and n.module and n.module.split(".")[0] == "fork":
            out.add(n.module)
            out |= {f"{n.module}.{a.name}" for a in n.names}
        elif isinstance(n, ast.Import):
            out |= {a.name for a in n.names if a.name.split(".")[0] == "fork"}
    return {m for m in out if _path(m)}


def reachable_modules() -> set[str]:
    seen: set[str] = set()
    todo = list(ENTRY_MODULES)
    while todo:
        m = todo.pop()
        if m not in seen:
            seen.add(m)
            todo += _imports(m)
    return seen


def _docstrings(tree: ast.Module) -> set[int]:
    kinds = (ast.Module, ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)
    return {id(b.body[0].value) for b in ast.walk(tree)
            if isinstance(b, kinds) and b.body and isinstance(b.body[0], ast.Expr)
            and isinstance(b.body[0].value, ast.Constant)}  # fmt: skip


class Index:
    """Functions of the reachable modules: names used, Neptune references, BODIES iteration."""

    def __init__(self, modules: set[str]):
        self.used: dict[tuple[str, str], set[str]] = {}
        self.neptune: dict[tuple[str, str], list[str]] = {}
        self.iterates_bodies: set[tuple[str, str]] = set()
        self.by_name: dict[str, set[tuple[str, str]]] = {}
        self.module_level: list[str] = []
        self.precovery: list[str] = []
        for m in sorted(modules):
            tree = _tree(m)
            docs = _docstrings(tree)
            for n in ast.walk(tree):
                if (isinstance(n, ast.Constant) and isinstance(n.value, str) and id(n) not in docs
                        and "precover" in n.value.lower()):  # fmt: skip
                    self.precovery.append(f"{m}:{n.lineno}")
            self._visit(m, tree, "", None, docs)

    def _visit(self, m, node, prefix, cls, docs):
        for ch in ast.iter_child_nodes(node):
            if isinstance(ch, (ast.FunctionDef, ast.AsyncFunctionDef)):
                q = f"{prefix}.{ch.name}" if prefix else ch.name
                key = (m, q)
                self.by_name.setdefault(ch.name, set()).add(key)
                if cls:
                    self.by_name.setdefault(cls, set()).add(key)  # a class name reaches its methods
                self.used[key] = set()
                self.neptune[key] = []
                for n in ast.walk(ch):
                    if isinstance(n, ast.Name):
                        self.used[key].add(n.id)
                    elif isinstance(n, ast.Attribute):
                        self.used[key].add(n.attr)
                    hit = _neptune_hit(n, docs)
                    if hit:
                        self.neptune[key].append(f"line {n.lineno}: {hit}")
                    if isinstance(n, (ast.For, ast.comprehension)) and _names_bodies(n.iter):
                        self.iterates_bodies.add(key)
                # Nested functions are part of the enclosing one (already walked above).
            elif isinstance(ch, ast.ClassDef):
                self._visit(m, ch, f"{prefix}.{ch.name}" if prefix else ch.name, ch.name, docs)

    def reach(self, entries: set[tuple[str, str]]) -> set[tuple[str, str]]:
        seen: set[tuple[str, str]] = set()
        todo = list(entries)
        while todo:
            f = todo.pop()
            if f in seen:
                continue
            seen.add(f)
            for name in self.used[f]:
                todo += self.by_name.get(name, ())
        return seen


def _neptune_hit(n: ast.AST, docs: set[int]) -> str | None:
    if isinstance(n, ast.Constant) and isinstance(n.value, str) and id(n) not in docs:
        if "neptune" in n.value.lower():
            return repr(n.value[:40])
    if isinstance(n, ast.Name) and n.id in NEPTUNE_NAMES:
        return n.id
    if isinstance(n, ast.Attribute) and n.attr in NEPTUNE_NAMES:
        return n.attr
    return None


def _names_bodies(expr: ast.AST) -> bool:
    """`BODIES`, `bodies.BODIES` or a call on either (`.values()`, `.items()`, `sorted(...)`)."""
    return any(
        (isinstance(n, ast.Name) and n.id == "BODIES")
        or (isinstance(n, ast.Attribute) and n.attr == "BODIES")
        for n in ast.walk(expr)
    )


def fit_path() -> tuple[Index, set[tuple[str, str]]]:
    idx = Index(reachable_modules())
    m6a_m6b = {k for k in idx.used if k[0] in ENTRY_MODULES[:3] and "." not in k[1]}
    m6c = {("fork.export.m6c", f) for f in M6C_FIT}
    assert m6c <= set(idx.used), "an M6c fit entry point was renamed; update M6C_FIT"
    return idx, idx.reach(m6a_m6b | m6c)


def test_no_fit_path_function_names_neptune_or_iterates_bodies():
    idx, path = fit_path()
    hits = {f"{m}.{q}": idx.neptune[(m, q)] for m, q in sorted(path) if idx.neptune[(m, q)]}
    assert not hits, f"Neptune referenced on the M6a–M6d fit path: {hits}"
    loops = sorted(f"{m}.{q}" for m, q in path & idx.iterates_bodies)
    assert not loops, f"fit-path functions iterate over bodies.BODIES (includes Neptune): {loops}"
    # The graph is not vacuous: it reaches the Horizons reader and the likelihood.
    names = {q for _, q in path}
    assert {"state_vectors_ecl_j2000", "dataset", "do_job"} <= names


def test_neptune_readers_exist_and_are_off_the_fit_path():
    """The check can see a Neptune read: the known ones are found, and none is on the path."""
    idx, path = fit_path()
    readers = {k for k, v in idx.neptune.items() if v}
    expected = {("fork.model.hidden_mass", "true_neptune_geocentric_longitude"),
                ("fork.data.horizons", "default_fetch_plan")}  # fmt: skip
    assert expected <= readers
    assert ("fork.data.horizons", "default_fetch_plan") in idx.iterates_bodies
    assert not readers & path


def test_the_graph_carries_a_neptune_read_to_its_caller():
    """From an M5a export that compares with the real Neptune, the graph reaches a reader: the
    same walk from the fit's entry points would have found one if it were there."""
    idx, _ = fit_path()
    readers = {k for k, v in idx.neptune.items() if v}
    assert idx.reach({("fork.export.to_json", "m5a_profile_data")}) & readers


def test_m6d_comparison_constants_stay_off_the_fit_path():
    idx, path = fit_path()
    named = sorted(
        f"{m}.{q}" for m, q in path if any(u.startswith("JPL_") for u in idx.used[(m, q)])
    )
    assert not named, f"a fit-path function names a JPL comparison constant: {named}"


def test_no_reference_to_a_precovery_file():
    idx, _ = fit_path()
    assert not idx.precovery, f"M6a–M6d modules mention a precovery file: {idx.precovery}"
