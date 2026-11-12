"""JPL Horizons reference states for the Phase 2 validation (V4): the eight planets at one
epoch and Mercury on a grid, 1800–2000, Sun-centred, ecliptic of J2000.

Cached under ``data/cache/phase2/`` (ignored), apart from the Phase 1 cache, so that nothing
Phase 1 reads is touched. A small subset is committed as a test fixture
(``tests/fixtures/phase2/v4_reference.npz``, built by ``build_fixture``).

    uv run python -m fork.phase2.reference --fetch      # network: fills the cache
    uv run python -m fork.phase2.reference --fixture    # writes the committed subset
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from fork.data import horizons as hz

ROOT = Path(__file__).resolve().parents[3]
CACHE = ROOT / "data" / "cache" / "phase2"
FIXTURE = ROOT / "tests" / "fixtures" / "phase2" / "v4_reference.npz"

PLANETS = ("mercury", "venus", "earth_moon", "mars", "jupiter", "saturn", "uranus", "neptune")
EPOCH_JD = 2415020.5  # 1900 January 1.0 TDB
JD_1800 = 2378496.5  # 1800 January 1.0
JD_2000 = 2451544.5  # 2000 January 1.0
GRID_STEP_DAYS = 20.3  # not commensurate with Mercury's 87.97 d
FIXTURE_STEP_DAYS = 1826.3  # about every 5 years: 41 epochs


def grid(step: float = GRID_STEP_DAYS) -> np.ndarray:
    """Julian dates (TDB) from 1800 to 2000 at ``step`` days, as offsets from the epoch that
    are exact multiples of the step."""
    k_lo = int(np.ceil((JD_1800 - EPOCH_JD) / step))
    k_hi = int(np.floor((JD_2000 - EPOCH_JD) / step))
    return EPOCH_JD + step * np.arange(k_lo, k_hi + 1)


def epoch_states(directory: Path = CACHE) -> np.ndarray:
    """(8, 6) Sun-centred states of the planets (system barycentres from Venus outward
    except as ``fork.model.bodies`` defines them) at ``EPOCH_JD``."""
    return np.vstack(
        [hz.state_vectors_ecl_j2000(b, EPOCH_JD, directory=directory) for b in PLANETS]
    )


def mercury_states(jd: np.ndarray, directory: Path = CACHE) -> np.ndarray:
    return hz.state_vectors_ecl_j2000("mercury", jd, directory=directory)


def load(path: Path | None = None) -> dict[str, np.ndarray]:
    """Reference set: ``epoch_states`` (8, 6), ``jd`` (T,), ``mercury`` (T, 6). From the
    cache by default (a miss goes to Horizons unless ``FORK_HORIZONS_CACHE_ONLY`` is set); from
    an ``.npz`` written by ``build_fixture`` when ``path`` is given."""
    if path is not None:
        with np.load(path) as d:
            return {k: d[k] for k in ("epoch_states", "jd", "mercury")}
    jd = grid()
    return {"epoch_states": epoch_states(), "jd": jd, "mercury": mercury_states(jd)}


def build_fixture() -> None:
    jd = grid(FIXTURE_STEP_DAYS)
    FIXTURE.parent.mkdir(parents=True, exist_ok=True)
    np.savez(FIXTURE, epoch_states=epoch_states(), jd=jd, mercury=mercury_states(jd))
    print(f"wrote {FIXTURE}: {jd.size} epochs")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--fetch", action="store_true")
    ap.add_argument("--fixture", action="store_true")
    a = ap.parse_args()
    CACHE.mkdir(parents=True, exist_ok=True)
    if a.fetch:
        d = load()
        print(f"epoch states {d['epoch_states'].shape}, Mercury at {d['jd'].size} epochs")
    if a.fixture:
        build_fixture()


if __name__ == "__main__":
    main()
