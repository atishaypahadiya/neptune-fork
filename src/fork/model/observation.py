"""M6 observation model: apparent RA/Dec measurements with per-group errors and offsets.

One measurement is an RA or a Dec at one epoch. The RA-only observations (Dec printed
"» » »") contribute an RA measurement and nothing else.

**Residuals**, in arcseconds on the sky:
- RA: wrap(α_obs − α_model) · cos δ_model;
- Dec: δ_obs − δ_model.

**Error groups** (``group_of``): observer series × quantity. The spec asks for observer ×
instrument × era, falling back to coarser groups when counts are small. Table 77 names no
instruments, and a group needs at least ``MIN_GROUP`` measurements of its quantity. The
groups are:

  - ``ancient`` (1690–1771): Flamsteed, Bradley, Mayer and Le Monnier together. Only Le
    Monnier reaches 10 alone (12 RA, 12 Dec); the others total 7 RA and 6 Dec, so the era is
    one group.
  - ``greenwich`` (1781–1800, Jan 1807–1809, 1828–1830), ``paris_1801_1827``,
    ``paris_1835_1845``: the series of p. 127.

**Likelihood.** Each measurement in group g has Gaussian error σ_g (fixed, from the error
model) about model + b_g. The group offset b_g ~ N(0, τ_g²) is a nuisance parameter,
marginalised analytically:

  −2 ln L = Σ_g [ Σ_{i∈g} r_i²/σ_g² − (Σ_{i∈g} r_i/σ_g²)² / (n_g/σ_g² + 1/τ_g²) ] + const.

This equals min over b of Σ (r_i − b_g)²/σ_g² + Σ b_g²/τ_g², so ``whitened`` returns that
augmented residual vector, which least-squares fits use directly. The constant
Σ_g ln(1 + n_g τ_g²/σ_g²) + Σ_i ln(2π σ_i²) depends only on the error model. It is dropped, as
M5b dropped its normalisation, so ln Z is comparable between models fitted to the same data and
error model only.
"""

from __future__ import annotations

import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from fork.model.apparent import ApparentGeometry, apparent_ra_dec
from fork.model.forward import ForwardContext, predict_states_ecl_j2000

ARCSEC_PER_RAD = 180.0 * 3600.0 / math.pi
TWO_PI = 2.0 * math.pi
I_URANUS = 2  # body index in predict_states_ecl_j2000 output: [J, S, U, (P)]
MIN_GROUP = 10
SERIES = ("ancient", "greenwich", "paris_1801_1827", "paris_1835_1845")
QUANTITIES = ("ra", "dec")
GROUPS = tuple(f"{s}:{q}" for s in SERIES for q in QUANTITIES)
ROOT = Path(__file__).resolve().parents[3]
ERROR_MODEL_FILE = ROOT / "data" / "historical" / "m6_error_model.json"


def series_of(row: dict[str, str]) -> str:
    """Observer series of a ``leverrier_1846.csv`` row (see module docstring)."""
    if row["era"] == "ancient":
        return "ancient"
    if row["observatory"] == "Greenwich":
        return "greenwich"
    if row["observatory"] != "Paris":
        raise ValueError(f"{row['id']}: modern row at {row['observatory']!r}")
    return "paris_1801_1827" if int(row["epoch"][:4]) <= 1830 else "paris_1835_1845"


def group_of(row: dict[str, str]) -> str:
    return f"{series_of(row)}:{row['quantity']}"


@dataclass(frozen=True)
class ErrorModel:
    """σ_g and τ_g in arcsec per group, with where they came from."""

    sigma_arcsec: dict[str, float]
    offset_prior_arcsec: dict[str, float]
    basis: str

    def __post_init__(self) -> None:
        for name in ("sigma_arcsec", "offset_prior_arcsec"):
            d = getattr(self, name)
            if set(d) != set(GROUPS) or not all(v > 0 for v in d.values()):
                raise ValueError(f"{name} must give a positive value for each of {GROUPS}")

    @classmethod
    def load(cls, path: Path = ERROR_MODEL_FILE) -> ErrorModel:
        d = json.loads(path.read_text(encoding="utf-8"))
        return cls(d["sigma_arcsec"], d["offset_prior_arcsec"], d["basis"])


@dataclass(frozen=True)
class RaDecObservations:
    """Measurements at ``geometry.jd_tt`` epochs (sorted, unique). Never carries the truth."""

    geometry: ApparentGeometry
    obs_id: np.ndarray  # (n_epoch,) str
    epoch_index: np.ndarray  # (m,) int: which epoch each measurement belongs to
    is_dec: np.ndarray  # (m,) bool
    value_rad: np.ndarray  # (m,) α or δ, true equator and equinox of date
    group_index: np.ndarray  # (m,) int into GROUPS
    sigma_arcsec: np.ndarray  # (len(GROUPS),)
    offset_prior_arcsec: np.ndarray  # (len(GROUPS),)

    @property
    def jd_tdb(self) -> np.ndarray:
        """Epochs (TT ≈ TDB to < 2 ms), for code that only needs the time span."""
        return self.geometry.jd_tt

    @property
    def n_measurements(self) -> int:
        return int(self.value_rad.size)

    def group_counts(self) -> np.ndarray:
        return np.bincount(self.group_index, minlength=len(GROUPS))

    def predict(self, theta: np.ndarray, context: ForwardContext) -> tuple[np.ndarray, np.ndarray]:
        """Model apparent (α, δ) of Uranus at every epoch."""
        pos, vel = predict_states_ecl_j2000(theta, self.geometry.jd_tt, context)
        return apparent_ra_dec(pos[:, I_URANUS], vel[:, I_URANUS], self.geometry)

    def residuals_arcsec(self, ra_model: np.ndarray, dec_model: np.ndarray) -> np.ndarray:
        """Observed − model for every measurement, arcsec on the sky (wrap-safe RA)."""
        ra_m = ra_model[self.epoch_index]
        dec_m = dec_model[self.epoch_index]
        d_ra = ((self.value_rad - ra_m + math.pi) % TWO_PI - math.pi) * np.cos(dec_m)
        d_dec = self.value_rad - dec_m
        return np.where(self.is_dec, d_dec, d_ra) * ARCSEC_PER_RAD

    def best_offsets_arcsec(self, resid_arcsec: np.ndarray) -> np.ndarray:
        """b_g minimising the augmented χ² (the posterior mean of the offsets given θ)."""
        g = self.group_index
        w = 1.0 / self.sigma_arcsec**2
        s = np.bincount(g, weights=resid_arcsec, minlength=len(GROUPS)) * w
        n = self.group_counts() * w + 1.0 / self.offset_prior_arcsec**2
        return s / n

    def whitened_from_residuals(self, resid_arcsec: np.ndarray) -> np.ndarray:
        b = self.best_offsets_arcsec(resid_arcsec)
        g = self.group_index
        return np.concatenate([(resid_arcsec - b[g]) / self.sigma_arcsec[g],
                               b / self.offset_prior_arcsec])  # fmt: skip

    def whitened(self, theta: np.ndarray, context: ForwardContext) -> np.ndarray:
        """Augmented residual vector; its squared sum is the marginal −2 ln L (+ const)."""
        return self.whitened_from_residuals(self.residuals_arcsec(*self.predict(theta, context)))

    def chi2(self, theta: np.ndarray, context: ForwardContext) -> float:
        return float(np.sum(self.whitened(theta, context) ** 2))

    def with_values(self, value_rad: np.ndarray) -> RaDecObservations:
        return RaDecObservations(self.geometry, self.obs_id, self.epoch_index, self.is_dec,
                                 np.asarray(value_rad, dtype=np.float64), self.group_index,
                                 self.sigma_arcsec, self.offset_prior_arcsec)  # fmt: skip


def read_rows(path: Path | None = None) -> list[dict[str, str]]:
    from fork.data.observations import OUT_FILE

    with (OUT_FILE if path is None else path).open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def structure_from_rows(
    rows: list[dict[str, str]],
    error_model: ErrorModel,
    geometry: ApparentGeometry | None = None,
    allow_empty_groups: bool = False,
) -> RaDecObservations:
    """The measurement layout of ``rows`` (epochs, RA/Dec presence, groups) with values left
    as NaN: what a synthetic data set is built on. Reads no observed value.

    ``allow_empty_groups`` (off by default, so every existing caller is unchanged) lets a
    data choice drop a whole series, e.g. the modern observations alone; a group that is
    present must still have ``MIN_GROUP`` measurements."""
    jd = np.array([float(r["epoch_tt_jd"]) for r in rows])
    uniq, inverse = np.unique(jd, return_inverse=True)
    ids = {}
    for r, k in zip(rows, inverse, strict=True):
        ids.setdefault(int(k), r["obs_id"])
    if len(ids) != uniq.size or len({r["obs_id"] for r in rows}) != uniq.size:
        raise ValueError("each epoch must belong to exactly one observation")
    geom = ApparentGeometry.from_horizons(uniq) if geometry is None else geometry
    if geom.jd_tt.shape != uniq.shape or np.any(geom.jd_tt != uniq):
        raise ValueError("geometry epochs do not match the rows")
    groups = np.array([GROUPS.index(group_of(r)) for r in rows])
    counts = np.bincount(groups, minlength=len(GROUPS))
    small = (counts < MIN_GROUP) & ((counts > 0) | (not allow_empty_groups))
    if np.any(small):
        raise ValueError(
            f"groups below {MIN_GROUP}: {[g for g, s in zip(GROUPS, small, strict=True) if s]}"
        )
    return RaDecObservations(
        geom, np.array([ids[k] for k in range(uniq.size)]), inverse.astype(np.int64),
        np.array([r["quantity"] == "dec" for r in rows]), np.full(len(rows), np.nan), groups,
        np.array([error_model.sigma_arcsec[g] for g in GROUPS], dtype=np.float64),
        np.array([error_model.offset_prior_arcsec[g] for g in GROUPS], dtype=np.float64),
    )  # fmt: skip


def true_uranus_ra_dec(structure: RaDecObservations) -> tuple[np.ndarray, np.ndarray]:
    """Noise-free apparent (α, δ) of the real Uranus (Horizons, includes the real Neptune)."""
    from fork.data.horizons import states_about_inner_barycenter_ecl_j2000 as about

    st = about("uranus", structure.geometry.jd_tt)
    return apparent_ra_dec(st[:, :3], st[:, 3:], structure.geometry)


def synthetic_values(
    structure: RaDecObservations,
    truth: tuple[np.ndarray, np.ndarray],
    rng: np.random.Generator,
    *,
    noise: bool = True,
) -> tuple[np.ndarray, np.ndarray]:
    """(values_rad, offsets_arcsec): truth + a per-group offset drawn from N(0, τ_g²) + per-
    measurement noise N(0, σ_g²), on the sky (RA displaced by (b + ε)/cos δ). ``noise=False``
    gives the bare truth (no offsets either)."""
    ra, dec = truth
    k = structure.epoch_index
    if noise:
        offsets = rng.normal(0.0, 1.0, len(GROUPS)) * structure.offset_prior_arcsec
        eps = rng.normal(0.0, 1.0, structure.n_measurements) * structure.sigma_arcsec[
            structure.group_index]  # fmt: skip
        shift = (offsets[structure.group_index] + eps) / ARCSEC_PER_RAD
    else:
        offsets = np.zeros(len(GROUPS))
        shift = np.zeros(structure.n_measurements)
    val = np.where(structure.is_dec, dec[k] + shift,
                   np.mod(ra[k] + shift / np.cos(dec[k]), TWO_PI))  # fmt: skip
    return val, offsets
