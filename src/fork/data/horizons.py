"""Heliocentric state vectors from JPL Horizons, with a local cache.

All states are **Sun-centered** (Horizons center ``500@10``), **ecliptic J2000** (Horizons
``refplane='ecliptic'``, i.e. ``ecl_j2000`` in ``fork.data.frames``), **geometric** (no light-time
or aberration), in AU and AU/day, at epochs given as JD **TDB** (Horizons interprets vector-table
epochs as TDB).

Cache: one ``data/cache/horizons_<body>.npz`` per body holding every epoch ever fetched. A request
fetches only epochs that are not already cached, so nothing is ever refetched.

CLI::

    uv run python -m fork.data.horizons --fetch
"""

from __future__ import annotations

import argparse
import os
import tempfile
from collections.abc import Callable
from pathlib import Path

import numpy as np
from numpy.typing import ArrayLike

from fork.model.bodies import BODIES, INNER_BODIES, M_CENTRAL_MSUN

# repo_root/src/fork/data/horizons.py → repo_root/data/cache
DEFAULT_CACHE_DIR = Path(__file__).resolve().parents[3] / "data" / "cache"

HORIZONS_CENTER = "500@10"  # Sun body center
HORIZONS_REFPLANE = "ecliptic"  # ecliptic & mean equinox of J2000 (IAU76 obliquity)
CHUNK = 50  # epochs per Horizons request (TLIST is sent in the URL)
_JD_MATCH_TOL_DAYS = 1e-6  # ~0.09 s, only for checking Horizons echoed the epochs we asked for

# fetcher(horizons_id, jd_tdb) -> state_ecl_j2000 (N, 6), rows aligned with jd_tdb
Fetcher = Callable[[str, np.ndarray], np.ndarray]


def cache_dir() -> Path:
    """Cache directory; override with env var ``FORK_CACHE_DIR``."""
    env = os.environ.get("FORK_CACHE_DIR")
    return Path(env) if env else DEFAULT_CACHE_DIR


CACHE_ONLY_ENV = "FORK_HORIZONS_CACHE_ONLY"


def cache_only() -> bool:
    """True when ``FORK_HORIZONS_CACHE_ONLY`` is set (to anything but "" or "0"): a cache miss
    then raises instead of going to Horizons. Pre-registered runs set it (M10.md, pre-run
    change (b)), so a run can never fetch an ephemeris silently."""
    return os.environ.get(CACHE_ONLY_ENV, "") not in ("", "0")


def cache_only_fetcher(body: str) -> Fetcher:
    """A fetcher that raises on any miss, naming the body and the missing epochs."""

    def fetch(horizons_id: str, jd_tdb: np.ndarray) -> np.ndarray:
        jd = np.asarray(jd_tdb, dtype=np.float64)
        shown = ", ".join(f"{v:.6f}" for v in jd[:5]) + (", …" if jd.size > 5 else "")
        raise RuntimeError(
            f"{CACHE_ONLY_ENV} is set and the cache misses {jd.size} epoch(s) for {body} "
            f"(Horizons id {horizons_id}): JD {shown}. Run `python -m fork.data.horizons "
            "--fetch` without the variable, or pass the epochs, before a cache-only run."
        )

    return fetch


def _cache_path(body: str, directory: Path) -> Path:
    return directory / f"horizons_{body}.npz"


def fetch_from_horizons(horizons_id: str, jd_tdb: np.ndarray) -> np.ndarray:
    """Query Horizons (network) for Sun-relative ecliptic-J2000 states at the given TDB JDs."""
    from astroquery.jplhorizons import Horizons  # deferred: network-only dependency path

    out = np.empty((jd_tdb.size, 6))
    for start in range(0, jd_tdb.size, CHUNK):
        chunk = jd_tdb[start : start + CHUNK]
        obj = Horizons(id=horizons_id, location=HORIZONS_CENTER, epochs=chunk.tolist())
        tab = obj.vectors(refplane=HORIZONS_REFPLANE, aberrations="geometric", cache=False)
        got_jd = np.asarray(tab["datetime_jd"], dtype=np.float64)
        try:
            rows = _match_rows(chunk, got_jd)
        except RuntimeError as err:
            raise RuntimeError(f"{err} for id={horizons_id}") from None
        cols = np.column_stack(
            [np.asarray(tab[c], dtype=np.float64) for c in ("x", "y", "z", "vx", "vy", "vz")]
        )
        out[start : start + chunk.size] = cols[rows]
    return out


def _match_rows(requested_jd: np.ndarray, got_jd: np.ndarray) -> np.ndarray:
    """Row of ``got_jd`` for each requested epoch, by nearest echoed JD.

    Horizons may reorder a TLIST, and it echoes a JD that is not exactly representable a few
    1e-10 d either side of the request, so the match is to the nearest echo, within
    ``_JD_MATCH_TOL_DAYS``. Raises if any epoch is missing, extra or unmatched.
    """
    req = np.asarray(requested_jd, dtype=np.float64)
    got = np.asarray(got_jd, dtype=np.float64)
    if got.size != req.size:
        raise RuntimeError("Horizons returned epochs that do not match the request")
    order = np.argsort(got)
    srt = got[order]
    hi = np.clip(np.searchsorted(srt, req), 0, srt.size - 1)
    lo = np.clip(hi - 1, 0, srt.size - 1)
    near = np.where(np.abs(srt[lo] - req) < np.abs(srt[hi] - req), lo, hi)
    if np.any(np.abs(srt[near] - req) > _JD_MATCH_TOL_DAYS) or np.unique(near).size != req.size:
        raise RuntimeError("Horizons returned epochs that do not match the request")
    return order[near]


def _load_cache(path: Path) -> tuple[np.ndarray, np.ndarray]:
    if not path.exists():
        return np.empty(0), np.empty((0, 6))
    with np.load(path) as f:
        return f["jd_tdb"], f["state_ecl_j2000"]


def _save_cache(path: Path, jd_tdb: np.ndarray, state: np.ndarray, horizons_id: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=path.parent, suffix=".npz.tmp")
    os.close(fd)
    try:
        with open(tmp, "wb") as fh:
            np.savez(
                fh,
                jd_tdb=jd_tdb,
                state_ecl_j2000=state,
                horizons_id=np.array(horizons_id),
                center=np.array(HORIZONS_CENTER),
                refplane=np.array(HORIZONS_REFPLANE),
            )
        os.replace(tmp, path)  # atomic: a crash never leaves a half-written cache
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)


def state_vectors_ecl_j2000(
    body: str,
    jd_tdb: ArrayLike,
    *,
    directory: Path | None = None,
    fetcher: Fetcher | None = None,
) -> np.ndarray:
    """Sun-relative state of ``body`` at each TDB JD: array ``(N, 6)`` = (x, y, z, vx, vy, vz).

    Frame ``ecl_j2000``; AU and AU/day. ``body`` is a key of ``fork.model.bodies.BODIES``
    (jupiter, saturn, uranus, neptune, earth). Epochs missing from the cache are fetched with
    ``fetcher`` (default: live Horizons) and merged into the cache; cached epochs never are.
    Epochs are matched exactly, so pass the same float JDs you will ask for later. With
    ``FORK_HORIZONS_CACHE_ONLY`` set and no ``fetcher``, a miss raises (``cache_only_fetcher``).
    """
    if body not in BODIES:
        raise KeyError(f"unknown body {body!r}; choose from {sorted(BODIES)}")
    directory = cache_dir() if directory is None else directory
    if fetcher is None:
        fetcher = cache_only_fetcher(body) if cache_only() else fetch_from_horizons
    horizons_id = BODIES[body].horizons_id

    jd_req = np.atleast_1d(np.asarray(jd_tdb, dtype=np.float64))
    if jd_req.ndim != 1:
        raise ValueError("jd_tdb must be a scalar or 1-D")

    path = _cache_path(body, directory)
    jd_cached, st_cached = _load_cache(path)

    missing = np.setdiff1d(jd_req, jd_cached)  # sorted, unique
    if missing.size:
        st_new = np.asarray(fetcher(horizons_id, missing), dtype=np.float64)
        if st_new.shape != (missing.size, 6):
            raise RuntimeError(
                f"fetcher returned shape {st_new.shape}, expected {(missing.size, 6)}"
            )
        jd_all = np.concatenate([jd_cached, missing])
        st_all = np.concatenate([st_cached, st_new])
        order = np.argsort(jd_all)
        jd_cached, st_cached = jd_all[order], st_all[order]
        _save_cache(path, jd_cached, st_cached, horizons_id)

    idx = np.searchsorted(jd_cached, jd_req)
    return st_cached[idx]


def inner_barycenter_state_ecl_j2000(
    jd_tdb: ArrayLike, *, directory: Path | None = None, fetcher: Fetcher | None = None
) -> np.ndarray:
    """State (N, 6) of the Sun + inner-planet barycenter relative to the Sun's center.

    This is the point the folded central mass (``M_CENTRAL_MSUN``) must sit at. It is ~3e-6 AU
    and ~7 cm/s from the Sun's center: tiny in position, but the velocity offset matters.
    """
    acc = None
    for body in INNER_BODIES:
        st = state_vectors_ecl_j2000(body.name, jd_tdb, directory=directory, fetcher=fetcher)
        acc = body.mass_msun * st if acc is None else acc + body.mass_msun * st
    return acc / M_CENTRAL_MSUN


def states_about_inner_barycenter_ecl_j2000(
    body: str, jd_tdb: ArrayLike, *, directory: Path | None = None, fetcher: Fetcher | None = None
) -> np.ndarray:
    """State (N, 6) of ``body`` relative to the Sun + inner-planet barycenter, ecl_j2000.

    Use these, not Sun-centred states, as initial conditions for any integration whose central
    mass is ``M_CENTRAL_MSUN``. Positions differ from heliocentric by ≤ 3e-6 AU (≤ 0.03″ in
    Uranus's longitude), so model outputs may be compared with heliocentric observations.
    """
    helio = state_vectors_ecl_j2000(body, jd_tdb, directory=directory, fetcher=fetcher)
    return helio - inner_barycenter_state_ecl_j2000(jd_tdb, directory=directory, fetcher=fetcher)


# --------------------------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------------------------

A1_EPOCH_ISO_TDB = "1781-03-13T00:00:00"  # Herschel's discovery of Uranus; A1 initial conditions
# Galle's discovery night, Berlin: the search ran around local midnight (≈ 23 h UT). Given in TDB;
# the ~6 s ΔT is irrelevant (Neptune moves ~1e-4″ in 6 s).
GALLE_EPOCH_ISO_TDB = "1846-09-23T23:00:00"


def default_fetch_plan() -> dict[str, np.ndarray]:
    """Epochs fetched by ``--fetch``: the A1 epoch for every body, a yearly 1690–1846 grid
    (Jan 1.0 TDB) for Uranus, Neptune and Earth, and Galle's discovery epoch (M5) for Earth,
    Neptune and the inner planets."""
    from fork.data.frames import jd_tdb_from_iso

    a1 = np.array([jd_tdb_from_iso(A1_EPOCH_ISO_TDB)])
    yearly = np.asarray(jd_tdb_from_iso([f"{y}-01-01T00:00:00" for y in range(1690, 1847)]))
    plan = {name: a1 for name in BODIES}
    for name in ("uranus", "neptune", "earth"):
        plan[name] = np.union1d(a1, yearly)
    galle = jd_tdb_from_iso(GALLE_EPOCH_ISO_TDB)
    for name in ("earth", "neptune", *(b.name for b in INNER_BODIES)):
        plan[name] = np.union1d(plan[name], [galle])
    return plan


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--fetch", action="store_true", help="fetch the default plan into cache")
    args = parser.parse_args(argv)
    if not args.fetch:
        parser.print_help()
        return
    for body, jds in default_fetch_plan().items():
        st = state_vectors_ecl_j2000(body, jds)
        r = np.linalg.norm(st[:, :3], axis=1)
        print(f"{body:8s} {jds.size:4d} epochs  |r| {r.min():7.3f}–{r.max():7.3f} AU")
    print(f"cache: {cache_dir()}")


if __name__ == "__main__":
    main()
