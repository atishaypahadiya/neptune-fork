"""Phase 1b step 3: the Adams-variant fit (docs/notes/P1B_adams_fit.md, pre-registered).

Data set D: Table 77 minus Flamsteed 1690, minus the Paris 1801–1827 series, minus every row
after 1840; the rows after 1840 are the test rows (the memoir's §§ 6 and 53 use the later
observations as a test). Three fits, a held at 38.4 AU (H1, Hypothesis I) and 37.25 AU (H2,
Hypothesis II) and a free (F), each under Gaussian (G, the headline's σ) and Student-t (T, ν = 4,
S3's scales) errors; everything else the M6c headline's. Jobs are variant × cell in M6c's
structure; ``summarize`` writes docs/notes/P1B_adams_fit.json and .png with every quantity of the
plan's § 5 and the reading rule of § 6 applied.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from fork.export import m6c
from fork.export.to_json import A1_EPOCH_JD_TDB

ROOT = Path(__file__).resolve().parents[3]
CACHE = ROOT / "data" / "cache" / "p1b_adams"
PROGRESS_LOG = CACHE / "progress.log"
DRIVER_LOG = CACHE / "driver.log"
NOTES = ROOT / "docs" / "notes"
OUT_JSON = NOTES / "P1B_adams_fit.json"
OUT_PNG = NOTES / "P1B_adams_fit.png"

FITS = ("h1", "h2", "f")
ERRORS = ("g", "t")
VARIANTS = ("h1_g", "h2_g", "f_g", "h1_t", "h2_t", "f_t")  # the launch order (§ 8)
HELD_A_AU = {"h1": 38.4, "h2": 37.25}  # Hypothesis I (a/a′ = 0.5), Hypothesis II (0.515)
E_MAX = 0.3
NU = 4.0
E_CELLS_HELD = ((0.0, 0.12), (0.12, 0.3))  # the headline's e split for a ≥ 28 AU
NLIVE_DEFAULT = 1000
NLIVE_FREE_CELL0 = 2000
WORKERS_HELD = 6  # two cells → 12
WORKERS_FREE = {0: 6, 1: 2, 2: 2, 3: 2}
SEED_NAMESPACE = (2026, 11, 3)
SUMMARY_RNG = (7, 1846, 11, 3)
DERIVED_MAX_SAMPLES = 6000
LAST_FIT_YEAR = 1840
TAU_ARCSEC = {"ancient": 5.0, "modern": 2.0}
A_BIN_HALF_WIDTH_AU = 0.5

FROZEN = {
    "table77_csv": ("data/historical/leverrier_1846.csv",
                    "2f546b9979f613c62d9c2dc77cd92b8269a9dbea4d3f87eecf39d2145509ebc4"),
    "error_model": ("data/historical/m6_error_model.json",
                    "0856fa30d1675bf926d3921580de919b45fb72474e2f41ba58a5d153a6bec277"),
    "adams_solutions": ("data/historical/adams_challis/adams_solutions.csv",
                        "9952542fa6283d7e0f67db58714173a0a17d746f41aa3cae7c2cd199088347ea"),
    "challis_sweeps": ("data/historical/adams_challis/challis_sweeps.csv",
                       "f9061de02545168636832645cca164260063ccf50ec1e1498b06950fd450fafa"),
}  # fmt: skip

# § 1's table, computed from the frozen CSV at pre-registration: (observations, measurements).
EXPECTED_FIT = {"total": (145, 278), "ancient": (18, 35), "greenwich": (87, 167),
                "paris_1801_1827": (0, 0), "paris_1835_1845": (40, 76)}  # fmt: skip
EXPECTED_TEST = (26, 47)

# § 4: report epochs; the strings are TDB (fork.data.horizons.P1B_REPORT_EPOCHS_ISO_TDB).
REPORT_EPOCHS = {
    "1845-09-30": "1845-09-30T23:30:00",  # S4's own date (end of September 1845)
    "1846-07-29": "1846-07-29T23:30:00",
    "1846-07-30": "1846-07-30T23:30:00",
    "1846-08-04": "1846-08-04T23:30:00",
    "1846-08-12": "1846-08-12T23:30:00",
    "1846-09-23": "1846-09-23T23:00:00",  # Galle's night, the headline's epoch
    "1846-09-29": "1846-09-29T23:30:00",
    "1846-10-01": "1846-10-01T23:30:00",  # Adams's date in the 2 September letter
}
CHALLIS_NIGHTS = ("1846-07-29", "1846-07-30", "1846-08-04", "1846-08-12", "1846-09-29")
READING_EPOCH = "1846-10-01"
GALLE_EPOCH = "1846-09-23"
S4_EPOCH = "1845-09-30"
S4_PRINTED_GEOCENTRIC_DEG = 320.5  # "320° 30′ nearly", end of September 1845 (S4)
S7_MEAN_LON_DEG = 315.0 + 20.0 / 60.0  # about 315° 20′ on 1 October 1846, the extrapolation
ADAMS_S53_ERRORS_ARCSEC = {"H1": {"1843": 7.11, "1844": 8.79, "1845": 12.40},
                           "H2": {"1843": 5.77, "1844": 7.05, "1845": 10.18}}  # fmt: skip


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_frozen() -> None:
    from fork.inference.nested import DLOGZ
    from fork.inference.priors import INC_PRIOR_MAX_DEG
    from fork.model.observation import ErrorModel

    for name, (rel, digest) in FROZEN.items():
        if sha256(ROOT / rel) != digest:
            raise RuntimeError(f"{name} ({rel}) differs from the pre-registered file")
    em = ErrorModel.load(ROOT / FROZEN["error_model"][0])
    for g, t in em.offset_prior_arcsec.items():
        if abs(t - TAU_ARCSEC["ancient" if g.startswith("ancient") else "modern"]) > 1e-12:
            raise RuntimeError(f"tau for {g} is {t}, not the pre-registered value")
    if abs(DLOGZ - 0.1) > 1e-12 or abs(INC_PRIOR_MAX_DEG - 10.0) > 1e-12:
        raise RuntimeError("sampler or prior settings differ from the pre-registration")


# ------------------------------------------------------------------ data set D


def _year(row: dict[str, str]) -> int:
    return int(row["epoch"][:4])


def adams_rows(rows: list[dict[str, str]]) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    """(fit rows, test rows) of Table 77 under Adams's choices (plan § 1, D1–D6)."""
    from fork.model.observation import series_of

    fit, test = [], []
    for r in rows:
        if r["observer"] == "Flamsteed" and _year(r) == 1690:  # D2: set aside, memoir § 25
            continue
        if series_of(r) == "paris_1801_1827":  # D4: Greenwich only, to 1830
            continue
        (fit if _year(r) <= LAST_FIT_YEAR else test).append(r)  # D5: the series ends at 1840
    return fit, test


def row_counts(rows: list[dict[str, str]]) -> dict[str, tuple[int, int]]:
    """(observations, measurements) in all and per series."""
    from fork.model.observation import SERIES, series_of

    out = {"total": (len({r["obs_id"] for r in rows}), len(rows))}
    for s in SERIES:
        sub = [r for r in rows if series_of(r) == s]
        out[s] = (len({r["obs_id"] for r in sub}), len(sub))
    return out


def split_variant(variant: str) -> tuple[str, str]:
    fit, err = variant.split("_")
    if fit not in FITS or err not in ERRORS:
        raise ValueError(f"unknown variant {variant!r}")
    return fit, err


def dataset(variant: str = "f_g"):
    """(obs, ctx, Uranus (i, Ω), (p, q), Uranus box, test obs): data set D with the frozen error
    model; the box by the M5 rule on D under the Gaussian model (as S3 used the headline's);
    the test rows in the same layout. Reads no Neptune state."""
    from fork.data.calibration import real_observations
    from fork.export.m9 import drop_measurements
    from fork.inference.fit import uranus_theta_to_x
    from fork.inference.fit_radec import uranus_box_radec
    from fork.model import observation as om
    from fork.model.forward import ForwardContext, uranus_inc_raan_at_epoch

    _, err = split_variant(variant)
    em = om.ErrorModel.load(ROOT / FROZEN["error_model"][0])
    rows = om.read_rows()
    fit_rows, test_rows = adams_rows(rows)
    # Table 77 in full, then the rows Adams did not use dropped (M9's way: an empty group keeps
    # its σ and τ, its offset term is 0), since the layout check wants every group populated.
    full = real_observations(rows, em)
    key = {(r["obs_id"], r["quantity"] == "dec") for r in fit_rows}
    tkey = {(r["obs_id"], r["quantity"] == "dec") for r in test_rows}
    ids = full.obs_id[full.epoch_index]
    keep = np.array([(i, bool(d)) in key for i, d in zip(ids, full.is_dec, strict=True)])
    tkeep = np.array([(i, bool(d)) in tkey for i, d in zip(ids, full.is_dec, strict=True)])
    obs = drop_measurements(full, keep)
    test = drop_measurements(full, tkeep)
    ctx = ForwardContext.from_horizons(A1_EPOCH_JD_TDB)
    inc_raan = uranus_inc_raan_at_epoch(A1_EPOCH_JD_TDB)
    box = uranus_box_radec(obs, ctx, inc_raan)
    pq = tuple(uranus_theta_to_x(np.array([19.2, 0.0, *inc_raan, 0.0, 0.0]))[4:])
    if err == "t":
        from fork.model.student_t import StudentTObservations

        obs = StudentTObservations.from_gaussian(obs, NU, m6c.t_scales_arcsec())
    return obs, ctx, inc_raan, pq, box, test


# ------------------------------------------------------------------ spaces and jobs


def variant_cells(variant: str) -> tuple:
    """F: the headline's four cells. H1, H2: two e cells, no a slicing (a is held)."""
    from fork.inference.nested import cells_for

    fit, _ = split_variant(variant)
    if fit == "f":
        return cells_for(E_MAX)
    return tuple((None, e) for e in E_CELLS_HELD)


def space_for(variant: str, cell: int, pq, box):
    """F: the headline's free-a space on the cell. H1/H2: ln a held (``PosteriorSpace.held``,
    the mechanism ``held_shape_space`` uses), the other four perturber entries and the two
    inclination entries free beside Uranus's four; the cell an e slice."""
    from fork.inference.nested import cell_space
    from fork.inference.priors import PosteriorSpace, free_a_space

    fit, _ = split_variant(variant)
    if fit == "f":
        return cell_space(free_a_space(pq, E_MAX, box, inclined=True), variant_cells(variant)[cell])
    held = PosteriorSpace((None, math.log(HELD_A_AU[fit]), None, None, None), e_max=E_MAX,
                          pq=pq, uranus_box=box, inclined=True)  # fmt: skip
    return held.with_e_slice(*E_CELLS_HELD[cell])


@dataclass(frozen=True)
class Job:
    variant: str
    cell: int

    @property
    def id(self) -> str:
        return f"{self.variant}/c{self.cell}"

    @staticmethod
    def from_id(job_id: str) -> Job:
        variant, part = job_id.split("/")
        return Job(variant, int(part.removeprefix("c")))

    @property
    def nlive(self) -> int:
        fit, _ = split_variant(self.variant)
        return NLIVE_FREE_CELL0 if fit == "f" and self.cell == 0 else NLIVE_DEFAULT

    @property
    def workers(self) -> int:
        fit, _ = split_variant(self.variant)
        return WORKERS_FREE[self.cell] if fit == "f" else WORKERS_HELD

    def path(self, suffix: str) -> Path:
        return CACHE / (self.id.replace("/", "__") + suffix)

    def rng_seed(self) -> list[int]:
        fit, err = split_variant(self.variant)
        return [*SEED_NAMESPACE, FITS.index(fit), ERRORS.index(err), self.cell + 1]

    def rng(self) -> np.random.Generator:
        return np.random.default_rng(self.rng_seed())


def all_jobs(variants: tuple[str, ...] = VARIANTS) -> list[Job]:
    return [Job(v, k) for v in variants for k in range(len(variant_cells(v)))]


def do_job(job: Job, n_workers: int, maxiter: int | None = None) -> None:
    from fork import ledger
    from fork.inference.nested import NestedResult, run_nested

    check_frozen()
    out = job.path(".pkl")
    if out.exists():
        print(f"{job.id}: cached")
        return
    obs, ctx, _, pq, box, _ = dataset(job.variant)
    space = space_for(job.variant, job.cell, pq, box)
    resumed = job.path(".ckpt").exists()
    t0, c0 = time.time(), ledger.cpu_now()
    res = run_nested(obs, ctx, space, job.rng(), nlive=job.nlive, n_workers=n_workers,
                     maxiter=maxiter, progress_log=PROGRESS_LOG,
                     checkpoint_file=None if maxiter else job.path(".ckpt"),
                     label=f"p1b_adams/{job.id}")  # fmt: skip
    if maxiter is not None:  # smoke runs are not results
        print(f"{job.id}: smoke run, {res.niter} it, {res.ncall} calls, {res.wall_s:.1f} s, "
              f"ndim {space.ndim}, names {space.names}")  # fmt: skip
        return
    NestedResult(res.space, res.results, res.nlive, res.dlogz, res.wall_s).save(out)
    yrs = ledger.sim_years_per_integration(obs.jd_tdb, A1_EPOCH_JD_TDB)
    ledger.append(f"p1b_adams/{job.id}", "nested", n_likelihood_calls=res.ncall,
                  n_integrations=res.ncall, sim_years=res.ncall * yrs,
                  wall_s=time.time() - t0, cpu_s=ledger.cpu_now() - c0, n_workers=n_workers,
                  estimated=resumed, resumed=resumed, nlive=job.nlive,
                  niter=res.niter)  # fmt: skip
    print(f"{job.id}: ln Z {res.logz:.3f} ± {res.logzerr:.3f}, {res.niter} it, "
          f"{res.wall_s / 60:.1f} min")  # fmt: skip


def _log(msg: str) -> None:
    CACHE.mkdir(parents=True, exist_ok=True)
    line = f"{time.strftime('%Y-%m-%d %H:%M:%S')} {msg}"
    with DRIVER_LOG.open("a") as f:
        f.write(line + "\n")
    print(line, flush=True)


def run_variant(variant: str, n_workers: int) -> None:
    """One variant's cells in parallel within ``n_workers``, largest first; one retry per failed
    job (it resumes from its checkpoint). M6c's scheduler, one variant at a time (§ 8)."""
    check_frozen()
    pending = sorted([j for j in all_jobs((variant,)) if not j.path(".pkl").exists()],
                     key=lambda j: (-j.nlive, j.cell))  # fmt: skip
    _log(f"run {variant}: {len(pending)} jobs pending, budget {n_workers} workers")
    running: dict[str, tuple[subprocess.Popen, float, int]] = {}
    failures: dict[str, int] = {}
    while pending or running:
        for jid, (proc, t0, _) in list(running.items()):
            rc = proc.poll()
            if rc is None:
                continue
            del running[jid]
            job = Job.from_id(jid)
            if rc == 0 and job.path(".pkl").exists():
                _log(f"done {jid} ({(time.time() - t0) / 60:.1f} min)")
            else:
                failures[jid] = failures.get(jid, 0) + 1
                _log(f"FAILED {jid} (exit {rc}, attempt {failures[jid]})")
                if failures[jid] < 2:
                    pending.insert(0, job)
        busy = sum(w for *_, w in running.values())
        while pending and (not running or busy + pending[0].workers <= n_workers):
            job = pending.pop(0)
            w = min(job.workers, n_workers)
            log = job.path(".log")
            log.parent.mkdir(parents=True, exist_ok=True)
            with log.open("a") as fh:
                proc = subprocess.Popen(
                    [sys.executable, "-m", "fork.export.p1b_adams", "job", job.id,
                     "--n-workers", str(w)],
                    stdout=fh, stderr=subprocess.STDOUT, cwd=ROOT)  # fmt: skip
            running[job.id] = (proc, time.time(), w)
            busy += w
            _log(f"start {job.id} (nlive {job.nlive}, {w} workers, pid {proc.pid})")
        time.sleep(30)
    _log(f"run {variant} finished; failures: {failures or 'none'}")


def plan() -> None:
    jobs = all_jobs()
    for j in jobs:
        state = ("done" if j.path(".pkl").exists() else
                 "resumable" if j.path(".ckpt").exists() else "pending")  # fmt: skip
        print(f"{j.id:8s} nlive {j.nlive:5d}  workers {j.workers:2d}  seed {j.rng_seed()}  {state}")
    print(f"{len(jobs)} jobs, {sum(j.path('.pkl').exists() for j in jobs)} done")


# ------------------------------------------------------------------ geometry at the report epochs


def mean_obliquity_deg(jd_tdb: float) -> float:
    """Mean obliquity of the ecliptic of date (Lieske et al. 1977, as the IAU 1976 precession
    this repository's frames use)."""
    t = (jd_tdb - 2451545.0) / 36525.0
    return (84381.448 - 46.8150 * t - 0.00059 * t * t + 0.001813 * t**3) / 3600.0


@dataclass(frozen=True)
class EpochGeometry:
    """What the derived quantities need at one report epoch: Earth there (about the inner
    barycenter, as the model's states), the precession ecl_j2000 → ecliptic of date and the
    mean obliquity of date. ``galle`` is the headline's GalleGeometry at this epoch."""

    label: str
    jd_tdb: float
    earth_pos_ecl_j2000_au: np.ndarray
    r_ecl_j2000_to_of_date: np.ndarray
    obliquity_deg: float

    @classmethod
    def from_horizons(cls, label: str, iso_tdb: str) -> EpochGeometry:
        from fork.data.frames import jd_tdb_from_iso, precession_matrix_ecl_j2000_to_ecl_of_date
        from fork.model.hidden_mass import earth_pos_about_inner_barycenter

        jd = float(jd_tdb_from_iso(iso_tdb))
        return cls(label, jd, earth_pos_about_inner_barycenter(jd),
                   np.asarray(precession_matrix_ecl_j2000_to_ecl_of_date(jd)),
                   mean_obliquity_deg(jd))  # fmt: skip

    @property
    def galle(self):
        from fork.inference.profile import GalleGeometry

        return GalleGeometry(self.jd_tdb, self.earth_pos_ecl_j2000_au, self.r_ecl_j2000_to_of_date)

    def derived(self, pos: np.ndarray, vel: np.ndarray, mu: float) -> dict[str, float]:
        """For a heliocentric state in ecl_j2000: geocentric ecliptic longitude of date with light
        time (the headline's quantity, ``GalleGeometry.lon_ecl_of_date_rad``); heliocentric true
        longitude of date; heliocentric mean longitude of date (Ω + ω + M of the osculating
        ellipse in the ecliptic of date); geocentric RA (hours) and Dec (degrees) referred to the
        mean equator and equinox of date (no nutation or aberration: both under 20″, far below
        the 9′ zone breadth they are compared with)."""
        from fork.data.frames import geocentric_lon_lat_ecl_j2000_light_time, state_to_elements

        geo_lon = math.degrees(self.galle.lon_ecl_of_date_rad(pos, vel))
        lon, lat, _ = geocentric_lon_lat_ecl_j2000_light_time(pos, vel, self.earth_pos_ecl_j2000_au)
        u = np.array([math.cos(lat) * math.cos(lon), math.cos(lat) * math.sin(lon), math.sin(lat)])
        w = self.r_ecl_j2000_to_of_date @ u
        eps = math.radians(self.obliquity_deg)
        y_eq = w[1] * math.cos(eps) - w[2] * math.sin(eps)
        z_eq = w[1] * math.sin(eps) + w[2] * math.cos(eps)
        ra_h = (math.degrees(math.atan2(y_eq, w[0])) % 360.0) / 15.0
        dec = math.degrees(math.atan2(z_eq, math.hypot(w[0], y_eq)))
        h = self.r_ecl_j2000_to_of_date @ pos
        helio_true = math.degrees(math.atan2(h[1], h[0])) % 360.0
        _, _, _, raan, argp, mean_anom = state_to_elements(h, self.r_ecl_j2000_to_of_date @ vel, mu)
        helio_mean = math.degrees(raan + argp + mean_anom) % 360.0
        return {"geo_lon": geo_lon, "helio_true": helio_true, "helio_mean": helio_mean,
                "ra_h": ra_h, "dec_deg": dec}  # fmt: skip


def report_geometries() -> dict[str, EpochGeometry]:
    return {k: EpochGeometry.from_horizons(k, v) for k, v in REPORT_EPOCHS.items()}


# ------------------------------------------------------------------ comparison targets (§ 3)


def _jd(iso: str) -> float:
    from fork.data.frames import jd_tdb_from_iso

    return float(jd_tdb_from_iso(iso))


def _published(label, a, e, mass, varpi_deg, mean_lon_deg, epoch_iso, sources):
    """A PublishedSolution from printed elements; its true longitude at the epoch is derived."""
    from fork.model.hidden_mass import PublishedSolution

    mean_anom = math.radians(mean_lon_deg - varpi_deg)
    ecc = mean_anom
    for _ in range(50):
        ecc -= (ecc - e * math.sin(ecc) - mean_anom) / (1.0 - e * math.cos(ecc))
    nu = 2.0 * math.atan2(math.sqrt(1 + e) * math.sin(ecc / 2),
                          math.sqrt(1 - e) * math.cos(ecc / 2))  # fmt: skip
    true_lon = (math.degrees(nu) + varpi_deg) % 360.0
    return PublishedSolution(label, a, e, mass, varpi_deg, mean_lon_deg, true_lon,
                             _jd(epoch_iso), sources)  # fmt: skip


LV_JUNE_SOURCES = (
    'Comptes rendus 22 (séance du 1er juin 1846), 907–918, at 917: "en assignant 325 degrés de '
    "longitude héliocentrique à la planète, au 1er janvier 1847, on ne commet pas une erreur "
    'de 10 degrés" (archive.org comptesrendusheb22acad, PDF page 931; P)',
    "Airy, MNRAS 7, 131–132 (summary; the fallback of amendment 2, not needed)",
)


def targets() -> dict:
    """The fixed targets of § 3. S4, S5, S6 and LV-June as Kepler orbits from the printed
    elements; LV-June has no elements (a circular orbit of 38.4 AU through its printed
    longitude at 1847.0; no mass); LV-Aug the repository's ``LE_VERRIER_1846``; S7 the
    extrapolation, a mean longitude only."""
    from fork.model.hidden_mass import LE_VERRIER_1846

    src = "adams_solutions.csv (step 1, P)"
    return {
        "S4": {"class": "communicated before the search", "date": "1845-09-18",
               "orbit": _published("S4 (18 Sep 1845, to Challis)", 38.4, 0.1428, 0.000173,
                                   320.5, 321.0 + 40.0 / 60.0, "1845-09-30T00:00:00", (src,)),
               "printed_geocentric_deg_1845_09_30": S4_PRINTED_GEOCENTRIC_DEG},
        "S5": {"class": "communicated before the search", "date": "1845-10-21",
               "orbit": _published("S5 (Hypothesis I)", 38.4, 0.16103, 0.0001656,
                                   315.0 + 57.0 / 60.0, 325.0 + 7.0 / 60.0,
                                   "1846-10-06T00:00:00", (src,))},
        "S6": {"class": "final computed", "date": "1846-09-02",
               "orbit": _published("S6 (Hypothesis II)", 37.25, 0.120615, 0.00015003,
                                   299.0 + 11.0 / 60.0, 323.0 + 2.0 / 60.0,
                                   "1846-10-06T00:00:00", (src,))},
        "S7": {"class": "extrapolation, not a solution", "date": "1846-09-02",
               "mean_lon_1846_10_01_deg": S7_MEAN_LON_DEG, "a_au_derived": 33.7,
               "note": "reported separately; enters no statement"},
        "LV-June": {"class": "comparison", "date": "1846-06-01",
                    "orbit": _published("Le Verrier, 1 June 1846", 38.4, 0.0, 0.0, 0.0, 325.0,
                                        "1847-01-01T00:00:00", LV_JUNE_SOURCES),
                    "printed": "325° heliocentric at 1 January 1847, 'not an error of 10°' (P)",
                    "note": "no elements printed: circular orbit of 38.4 AU through 325° at "
                            "1847.0 (0.38° of motion over the quarter)"},
        "LV-Aug": {"class": "comparison", "date": "1846-08-31", "orbit": LE_VERRIER_1846},
    }  # fmt: skip


def published_state_ecl_j2000(sol, jd_tdb: float, inc_raan: tuple[float, float]
                              ) -> tuple[np.ndarray, np.ndarray]:  # fmt: skip
    """A published orbit as a heliocentric Kepler ellipse, angles in the ecliptic and mean
    equinox of its epoch (assumed, as M5 and M6c), carried to ``jd_tdb``, rotated to ecl_j2000
    and laid in Uranus's plane (§ 3; the model's convention for a perturber with i_rel = 0:
    ``uranus_frame_matrix``). ``m6c.le_verrier_state_ecl_j2000`` generalised."""
    from fork.data.frames import precession_matrix_ecl_j2000_to_ecl_of_date
    from fork.model.bodies import M_CENTRAL_MSUN, G
    from fork.model.forward import uranus_frame_matrix

    mu = G * (M_CENTRAL_MSUN + sol.mass_msun)
    a, e = sol.a_au, sol.e
    varpi = math.radians(sol.varpi_ecl_of_date_deg)
    mean_anom = math.radians(sol.mean_lon_ecl_of_date_deg) - varpi + math.sqrt(mu / a**3) * (
        jd_tdb - sol.epoch_jd)  # fmt: skip
    ecc = mean_anom
    for _ in range(50):
        ecc -= (ecc - e * math.sin(ecc) - mean_anom) / (1.0 - e * math.cos(ecc))
    half = ecc / 2.0
    nu = 2.0 * math.atan2(math.sqrt(1 + e) * math.sin(half), math.sqrt(1 - e) * math.cos(half))
    r, lon, p = a * (1 - e * math.cos(ecc)), nu + varpi, a * (1 - e * e)
    vr, vt = math.sqrt(mu / p) * e * math.sin(nu), math.sqrt(mu * p) / r
    pos = np.array([r * math.cos(lon), r * math.sin(lon), 0.0])
    vel = np.array([vr * math.cos(lon) - vt * math.sin(lon),
                    vr * math.sin(lon) + vt * math.cos(lon), 0.0])  # fmt: skip
    rot = precession_matrix_ecl_j2000_to_ecl_of_date(sol.epoch_jd)
    f = np.asarray(uranus_frame_matrix(*inc_raan))
    return f @ (rot.T @ pos), f @ (rot.T @ vel)


def target_places(geos: dict[str, EpochGeometry], inc_raan: tuple[float, float]) -> dict:
    """Every target's derived quantities at every report epoch, and Neptune's from the cache."""
    from fork.data.horizons import states_about_inner_barycenter_ecl_j2000
    from fork.model.bodies import M_CENTRAL_MSUN, NEPTUNE, G

    out: dict = {}
    for name, t in targets().items():
        if "orbit" not in t:
            continue
        sol = t["orbit"]
        mu = G * (M_CENTRAL_MSUN + sol.mass_msun)
        out[name] = {k: g.derived(*published_state_ecl_j2000(sol, g.jd_tdb, inc_raan), mu)
                     for k, g in geos.items()}  # fmt: skip
    mu_n = G * (M_CENTRAL_MSUN + NEPTUNE.mass_msun)
    out["neptune"] = {}
    for k, g in geos.items():
        st = states_about_inner_barycenter_ecl_j2000("neptune", [g.jd_tdb])[0]
        out["neptune"][k] = g.derived(st[:3], st[3:], mu_n)
    return out


# ------------------------------------------------------------------ summary (§ 5)


def _circ(x: np.ndarray, nd: int = 3) -> dict:
    from fork.inference.mcmc import circular_interval_deg

    lo95, lo, med, hi, hi95 = circular_interval_deg(x, (2.5, 16.0, 50.0, 84.0, 97.5))
    return {"median": round(med, nd), "lo68": round(lo, nd), "hi68": round(hi, nd),
            "lo95": round(lo95, nd), "hi95": round(hi95, nd)}  # fmt: skip


def status_of(vs: dict) -> str:
    """§ 6: inside (68 %), marginal (68–95 %) or outside (beyond 95 %)."""
    return "inside" if vs["in_68"] else "marginal" if vs["in_95"] else "outside"


def derived_samples(space, ctx, eq: np.ndarray, geos: dict[str, EpochGeometry]) -> dict:
    """Per sample and report epoch, the five derived quantities, one integration per sample."""
    from fork.model.bodies import M_CENTRAL_MSUN, M_EARTH_MSUN, G
    from fork.model.forward import predict_states_ecl_j2000

    labels = list(geos)
    jds = [geos[k].jd_tdb for k in labels]
    keys = ("geo_lon", "helio_true", "helio_mean", "ra_h", "dec_deg")
    out = {k: {q: np.empty(len(eq)) for q in keys} for k in labels}
    for i, v in enumerate(eq):
        theta = space.theta(v)
        mu = G * (M_CENTRAL_MSUN + 10.0 ** theta[6] * M_EARTH_MSUN)
        pos, vel = predict_states_ecl_j2000(theta, jds, ctx)
        for j, k in enumerate(labels):
            d = geos[k].derived(pos[j, 3], vel[j, 3], mu)
            for q in keys:
                out[k][q][i] = d[q]
    return out


def _mass_vs(q: dict, tm: float) -> dict:
    half = (q["hi68"] - q["median"]) if tm > q["median"] else (q["median"] - q["lo68"])
    off = q["median"] - tm
    return {"target_mass_msun": tm, "median_minus_target_msun": round(off, 8),
            "offset_half_widths": round(off / half, 3) if half > 0 else None,
            "in_68": bool(q["lo68"] <= tm <= q["hi68"]),
            "in_95": bool(q["lo95"] <= tm <= q["hi95"])}  # fmt: skip


def a_density(a: np.ndarray) -> dict:
    """§ 5.8: F's posterior density in a at 38.4 and 37.25 AU (bins of ±0.5 AU, sample counts
    given), relative to the fullest 1 AU bin."""
    grid = np.arange(20.5, 60.0, 1.0)
    counts = np.array([np.sum(np.abs(a - c) <= A_BIN_HALF_WIDTH_AU) for c in grid])
    at = {}
    for name, c in (("H1", HELD_A_AU["h1"]), ("H2", HELD_A_AU["h2"]), ("neptune", 30.07)):
        n = int(np.sum(np.abs(a - c) <= A_BIN_HALF_WIDTH_AU))
        at[name] = {"a_au": c, "n_in_bin": n, "fraction_in_bin": round(n / len(a), 5),
                    "relative_to_max_bin": round(n / max(int(counts.max()), 1), 4)}  # fmt: skip
    return {"bin_half_width_au": A_BIN_HALF_WIDTH_AU, "n_samples": len(a),
            "max_bin_centre_au": float(grid[int(np.argmax(counts))]), "at": at}  # fmt: skip


def test_row_report(test, space, ctx, eq: np.ndarray, years: dict[str, str]) -> dict:
    """§ 5.6: posterior-predictive O − C of the Paris 1841–1845 measurements (median and 68 %
    per measurement), and per year the mean of the medians, beside Adams's § 53 errors."""
    tr = np.array([test.residuals_arcsec(*test.predict(space.theta(v), ctx)) for v in eq])
    rows = []
    for m in range(test.n_measurements):
        q = np.percentile(tr[:, m], [16, 50, 84])
        oid = str(test.obs_id[test.epoch_index[m]])
        rows.append({"obs_id": oid, "year": years[oid],
                     "quantity": "dec" if test.is_dec[m] else "ra",
                     "o_minus_c_arcsec": {"lo68": round(float(q[0]), 2),
                                          "median": round(float(q[1]), 2),
                                          "hi68": round(float(q[2]), 2)}})  # fmt: skip
    by_year: dict[str, dict] = {}
    for qty in ("ra", "dec"):
        for r in rows:
            if r["quantity"] == qty:
                by_year.setdefault(r["year"], {}).setdefault(qty, []).append(
                    r["o_minus_c_arcsec"]["median"]
                )
    summary = {}
    for y, v in sorted(by_year.items()):
        summary[y] = {}
        for qty in ("ra", "dec"):
            vals = v.get(qty, [])
            summary[y][f"{qty}_n"] = len(vals)
            summary[y][f"{qty}_mean_median_arcsec"] = (round(float(np.mean(vals)), 2)
                                                        if vals else None)  # fmt: skip
    return {"n_measurements": test.n_measurements, "n_samples": len(eq),
            "convention": "O − C in arcsec; RA residuals as the observation model scales them; "
                          "group offsets not subtracted (the fit's best-fit offsets are listed "
                          "under best_fit)",
            "by_year": summary, "rows": rows,
            "adams_s53_errors_arcsec": ADAMS_S53_ERRORS_ARCSEC}  # fmt: skip


def summarize_variant(variant: str, geos: dict[str, EpochGeometry], tplaces: dict,
                      years: dict[str, str]) -> tuple[dict, dict] | None:  # fmt: skip
    """(§ 5 record, raw arrays for the figure) for one variant, or None if a job is missing."""
    from fork.export.m5b import _e_summary
    from fork.inference.nested import ERR_INFLATION, NestedResult, PartitionedResult
    from fork.model.bodies import M_EARTH_MSUN
    from fork.model.observation import GROUPS

    jobs = all_jobs((variant,))
    if not all(j.path(".pkl").exists() for j in jobs):
        return None
    res = [NestedResult.load(j.path(".pkl")) for j in jobs]
    fit, err = split_variant(variant)
    rng = np.random.default_rng([*SUMMARY_RNG, VARIANTS.index(variant)])
    part = PartitionedResult(tuple(res), err_scale=ERR_INFLATION)
    eq = part.equal_weight(rng)
    obs, ctx, _, _, _, test = dataset(variant)
    space = res[0].space
    idx = rng.choice(len(eq), size=min(DERIVED_MAX_SAMPLES, len(eq)), replace=False)
    der = derived_samples(space, ctx, eq[idx], geos)
    thetas = np.array([space.theta(v) for v in eq])
    best = max(res, key=lambda r: r.max_logl)
    th_best = space.theta(best.samples[int(np.argmax(best.results["logl"]))])
    resid = obs.residuals_arcsec(*obs.predict(th_best, ctx))

    epochs = {}
    for k in geos:
        tg_geo = {n: tplaces[n][k]["geo_lon"] for n in tplaces}
        tg_true = {n: tplaces[n][k]["helio_true"] for n in tplaces}
        tg_mean = {n: tplaces[n][k]["helio_mean"] for n in tplaces}
        if k == READING_EPOCH:
            tg_mean["S7"] = S7_MEAN_LON_DEG
        if k == S4_EPOCH:
            tg_geo["S4_printed"] = S4_PRINTED_GEOCENTRIC_DEG
            tg_mean["S4_printed"] = targets()["S4"]["orbit"].mean_lon_ecl_of_date_deg
        rec = {"jd_tdb": geos[k].jd_tdb,
               "geocentric_lon_of_date_deg": m6c.lon_summary(der[k]["geo_lon"], tg_geo),
               "helio_true_lon_of_date_deg": m6c.lon_summary(der[k]["helio_true"], tg_true),
               "helio_mean_lon_of_date_deg": m6c.lon_summary(der[k]["helio_mean"],
                                                             tg_mean)}  # fmt: skip
        if k in CHALLIS_NIGHTS:
            d16, d84 = np.percentile(der[k]["dec_deg"], [16, 84])
            rec["ra_dec_mean_of_date"] = {
                "ra_h": m6c._quantiles(der[k]["ra_h"], 4),
                "dec_deg": m6c._quantiles(der[k]["dec_deg"], 3),
                "dec_width68_arcmin": round(60.0 * float(d84 - d16), 1),
                "zone_breadth_arcmin": 9.0,
                "targets": {n: {"ra_h": round(tplaces[n][k]["ra_h"], 4),
                                "dec_deg": round(tplaces[n][k]["dec_deg"], 3)} for n in tplaces},
            }  # fmt: skip
        epochs[k] = rec

    out = {
        "variant": variant, "fit": fit, "error": err, "held_a_au": HELD_A_AU.get(fit),
        "ndim": space.ndim, "n_samples": len(eq), "n_derived": len(idx),
        "evidence": {"logz": round(part.logz, 4), "logzerr": round(part.logzerr, 4),
                     "logzerr_inflation": ERR_INFLATION,
                     "max_logl": round(float(best.max_logl), 4),
                     "cells": [{"a_au": None if c[0] is None else list(c[0]),
                                "e": None if c[1] is None else list(c[1]), "nlive": r.nlive,
                                "share": round(float(s), 5), "share_err": round(float(se), 5),
                                "logz_k": round(r.logz, 4), "logzerr_k": round(r.logzerr, 4),
                                "niter": int(r.niter), "ncall": int(r.ncall),
                                "wall_min": round(r.wall_s / 60, 1),
                                "minus2_lnl_min": round(-2 * float(r.max_logl), 3)}
                               for c, r, s, se in zip(variant_cells(variant), res, part.shares,
                                                      part.share_errs, strict=True)]},
        "epochs": epochs,
        "mass_earth": m6c._quantiles(10.0 ** thetas[:, 6], 3),
        "mass_msun": m6c._quantiles(10.0 ** thetas[:, 6] * M_EARTH_MSUN, 8),
        "a_au": m6c._quantiles(thetas[:, 7], 3),
        "e": _e_summary(thetas[:, 8], None, E_MAX),
        "varpi_uranus_frame_deg": _circ(np.degrees(thetas[:, 9])),
        "inc_rel_uranus_deg": m6c._quantiles(np.degrees(thetas[:, 11]), 3)
        | {"prior_16_50_84": m6c.inc_prior_quantiles_deg()},
        "best_fit": {"minus2_lnl": round(-2 * float(best.max_logl), 3),
                     "statistic": "chi2 (augmented, constant dropped)" if err == "g" else "Q",
                     "n_measurements": int(obs.n_measurements),
                     "group_offsets_arcsec": dict(zip(
                         GROUPS, (round(float(v), 3) for v in obs.best_offsets_arcsec(resid)),
                         strict=True))},
        "test_rows": test_row_report(test, space, ctx, eq[idx], years),
    }  # fmt: skip
    if fit in HELD_A_AU:
        tname = "S5" if fit == "h1" else "S6"
        out["mass_vs_target"] = {"target": tname} | _mass_vs(
            out["mass_msun"], targets()[tname]["orbit"].mass_msun
        )
    else:
        out["a_density"] = a_density(thetas[:, 7])
    raw = {"lon_1oct_deg": der[READING_EPOCH]["geo_lon"], "a_au": thetas[:, 7]}
    return out, raw


# ------------------------------------------------------------ reading rule (§ 6, amendments 3–4)


def band_of(delta_lnl: float) -> str:
    """Statement 3's bands on Δ ln L_max(held) − ln L_max(F): M9's thresholds, with amendment
    3's middle band."""
    if delta_lnl <= -2.0:
        return "disfavoured"
    if delta_lnl <= -0.5:
        return "weakly disfavoured; reported, no statement"
    return "not distinguished"


def statement_1(f: dict, held: dict[str, dict]) -> dict:
    s = status_of(f["epochs"][READING_EPOCH]["geocentric_lon_of_date_deg"]["vs"]["neptune"])
    out = {"F": {"status": s, "reading": {"inside": "yes", "marginal": "partly",
                                           "outside": "no"}[s]}}  # fmt: skip
    for name, h in held.items():
        out[name] = {"status": status_of(
            h["epochs"][READING_EPOCH]["geocentric_lon_of_date_deg"]["vs"]["neptune"]),
            "note": "under his assumed distance"}  # fmt: skip
    return out


def statement_2(held: dict[str, dict]) -> dict:
    out = {}
    for tname, hname in (("S4", "H1"), ("S5", "H1"), ("S6", "H2")):
        h = held.get(hname)
        if h is None:
            continue
        e = h["epochs"][READING_EPOCH]
        rec = {"under": hname,
               "true_lon": status_of(e["helio_true_lon_of_date_deg"]["vs"][tname]),
               "mean_lon": status_of(e["helio_mean_lon_of_date_deg"]["vs"][tname])}  # fmt: skip
        rec["differ"] = rec["true_lon"] != rec["mean_lon"]
        out[tname] = rec
    return out


def statement_3(f: dict, held: dict[str, dict]) -> dict:
    out = {"density_of_F": f["a_density"]["at"]}
    for name, h in held.items():
        d = h["evidence"]["max_logl"] - f["evidence"]["max_logl"]
        out[name] = {"delta_lnl_max_vs_F": round(d, 3), "band": band_of(d)}
    return out


def interpret(s1: dict, s2: dict, s3: dict) -> dict:
    """The interpretation table's row for a reading, with amendment 4's two cases. Statement 2
    is read on the true longitude (where the planet is); the mean-longitude reading is reported
    beside it by ``statement_2``."""
    if s1["F"]["reading"] != "yes":
        return {
            "row": 4,
            "text": "His data, as reconstructed here, did not locate the planet; "
            "the comparison of his solutions is reported as description and no statement "
            "about his method is made.",
        }
    if "S5" not in s2 or "S6" not in s2:
        return {"row": None, "text": "statement 2 incomplete (a held-a run is missing)"}
    s5, s6 = s2["S5"]["true_lon"], s2["S6"]["true_lon"]
    if s5 != s6 and "inside" in (s5, s6):
        return {
            "row": "amendment 4 (S5 and S6 disagree)",
            "text": "S5 and S6 disagree; the "
            "reading is reported per solution and no single statement about his solutions "
            "is made.",
        }
    if "marginal" in (s5, s6):
        return {
            "row": "amendment 4 (marginal)",
            "text": "at the edge of what his data allowed "
            "under his assumption; no row applied, no statement about his method.",
        }
    if s5 == "outside" and s6 == "outside":
        return {
            "row": 3,
            "text": "His data located the planet but his solutions did not "
            "represent its posterior under his own assumption: the gap is in his solution "
            "method or data handling (A1–A5 are the candidates), not in the data.",
        }
    bands = {k: v["band"] for k, v in s3.items() if k in ("H1", "H2")}
    if set(bands.values()) == {"disfavoured"}:
        return {
            "row": 1,
            "text": "His data located the planet; his solutions sat where his "
            "data put it under his assumed distance; the distance was the error, and the "
            "data could have told him.",
        }
    if set(bands.values()) == {"not distinguished"}:
        return {
            "row": 2,
            "text": "His data located the planet and his solutions represented "
            "them; the data of 1840 did not fix the distance, so his assumption was a free "
            "choice within what they allowed.",
        }
    return {
        "row": "amendment 3",
        "text": "statement 3 is in the weakly-disfavoured band or "
        f"mixed ({bands}); rows 1–2 not applied, no statement.",
    }


def reading(fits: dict) -> dict:
    """§ 6 on the Gaussian runs, the Student-t runs beside them; where the two disagree on a
    target's status on 1 October 1846, both are given and the Gaussian reading stands."""
    out: dict = {
        "epoch": READING_EPOCH,
        "rule": "inside ≤ 68 %, marginal 68–95 %, outside > "
        "95 %; statement 1 on the geocentric longitude, statement 2 on the true and "
        "mean heliocentric longitudes, statement 3 on Δ ln L_max and F's density in a",
    }
    for err in ERRORS:
        f = fits.get(f"f_{err}")
        held = {n: fits[f"{k}_{err}"] for n, k in (("H1", "h1"), ("H2", "h2"))
                if f"{k}_{err}" in fits}  # fmt: skip
        rec: dict = {}
        if f is not None:
            rec["statement_1"] = statement_1(f, held)
        if held:
            rec["statement_2"] = statement_2(held)
        if f is not None and held:
            rec["statement_3"] = statement_3(f, held)
        comp = {}
        for name, h in (("F", f), *held.items()):
            if h is not None:
                vs = h["epochs"][READING_EPOCH]["geocentric_lon_of_date_deg"]["vs"]
                comp[name] = {n: status_of(vs[n]) for n in ("LV-June", "LV-Aug")}
        rec["comparison_le_verrier"] = comp
        if {"statement_1", "statement_2", "statement_3"} <= rec.keys():
            rec["interpretation"] = interpret(rec["statement_1"], rec["statement_2"],
                                              rec["statement_3"])  # fmt: skip
        out[err] = rec
    disagree = []
    for v in VARIANTS[:3]:
        g, t = fits.get(v), fits.get(v.replace("_g", "_t"))
        if g is None or t is None:
            continue
        vg = g["epochs"][READING_EPOCH]["geocentric_lon_of_date_deg"]["vs"]
        vt = t["epochs"][READING_EPOCH]["geocentric_lon_of_date_deg"]["vs"]
        for n in vg:
            if status_of(vg[n]) != status_of(vt[n]):
                disagree.append({"fit": v[:-2], "target": n, "gaussian": status_of(vg[n]),
                                 "student_t": status_of(vt[n])})  # fmt: skip
    out["gaussian_vs_student_t_disagreements_1846_10_01"] = disagree
    out["note"] = "the Gaussian reading is the pre-registered one"
    return out


# ------------------------------------------------------------------ driver


def _json_default(o):
    if isinstance(o, np.floating | np.integer):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    raise TypeError(type(o))


def _headline_galle() -> dict | None:
    p = ROOT / "site" / "data" / "m6_fits.json"
    if not p.exists():
        return None
    g = json.loads(p.read_text(encoding="utf-8"))["fits"]["headline"]["galle_lon_ecl_of_date_deg"]
    return {k: g[k] for k in ("median", "lo68", "hi68", "width68")}


def _target_record(t: dict) -> dict:
    rec = {k: v for k, v in t.items() if k != "orbit"}
    if "orbit" in t:
        o = t["orbit"]
        rec["elements"] = {"a_au": o.a_au, "e": o.e, "mass_msun": o.mass_msun,
                           "varpi_deg": o.varpi_ecl_of_date_deg,
                           "mean_lon_deg": o.mean_lon_ecl_of_date_deg,
                           "true_lon_deg_at_epoch": round(o.true_lon_ecl_of_date_deg, 4),
                           "epoch_jd": o.epoch_jd, "sources": list(o.sources)}  # fmt: skip
    return rec


def summarize(png: bool = True) -> dict:
    from fork import ledger
    from fork.inference.nested import BOUND, DLOGZ, ERR_INFLATION, SAMPLE
    from fork.model import observation as om
    from fork.model.forward import uranus_inc_raan_at_epoch

    check_frozen()
    t0, c0 = time.time(), ledger.cpu_now()
    geos = report_geometries()
    inc_raan = uranus_inc_raan_at_epoch(A1_EPOCH_JD_TDB)
    tplaces = target_places(geos, inc_raan)
    years = {r["obs_id"]: r["epoch"][:4] for r in om.read_rows()}
    fits, raw = {}, {}
    for v in VARIANTS:
        s = summarize_variant(v, geos, tplaces, years)
        if s is not None:
            fits[v], raw[v] = s
    head = _headline_galle()
    for v in ("f_g", "f_t"):
        if v in fits and head is not None:
            g = fits[v]["epochs"][GALLE_EPOCH]["geocentric_lon_of_date_deg"]
            shift = (g["median"] - head["median"] + 180.0) % 360.0 - 180.0
            fits[v]["vs_headline"] = {
                "galle_median_shift_deg": round(shift, 3),
                "galle_median_shift_headline_half_widths": round(
                    shift / (0.5 * head["width68"]), 3),
                "width68_ratio": round(g["width68"] / head["width68"], 3),
                "headline_galle_lon_deg": head,
                "note": "no ln B against the headline: the data differ"}  # fmt: skip
    payload = {
        "description": "Phase 1b step 3: the Adams-variant fit (docs/notes/P1B_adams_fit.md). "
        "Data set D (Table 77 minus Flamsteed 1690, Paris 1801-1827 and rows after 1840), three "
        "fits (a held at 38.4 and 37.25 AU, a free) under Gaussian and Student-t (nu 4) errors; "
        "M6c's observation model, priors and sampler. ln Z errors x err_inflation.",
        "frozen_inputs": {k: {"path": p, "sha256": h} for k, (p, h) in FROZEN.items()},
        "data_set": {"fit": EXPECTED_FIT, "test": EXPECTED_TEST,
                     "note": "(observations, measurements); the test rows are Paris 1841–1845"},
        "sampler": {"bound": BOUND, "sample": SAMPLE, "dlogz": DLOGZ,
                    "nlive": {j.id: j.nlive for j in all_jobs()}, "err_inflation": ERR_INFLATION,
                    "rng": {j.id: j.rng_seed() for j in all_jobs()},
                    "cells": {v: [{"a_au": None if a is None else list(a),
                                   "e": None if e is None else list(e)}
                                  for a, e in variant_cells(v)] for v in VARIANTS}},
        "priors": {"log10_m_earth": [0.0, 2.0], "a_au_log_uniform": [20.0, 60.0],
                   "a_held_au": HELD_A_AU, "e_uniform_max": E_MAX,
                   "inc_rel_uranus_deg_max": 10.0, "inc": "uniform in cos i",
                   "angles": "uniform", "uranus": "box ±50 scaled σ around the Uranus-only "
                   "RA/Dec fit to D"},
        "report_epochs": {k: {"iso_tdb": v, "jd_tdb": geos[k].jd_tdb}
                          for k, v in REPORT_EPOCHS.items()},
        "targets": {n: _target_record(t) | {"places": tplaces.get(n)}
                    for n, t in targets().items()},
        "neptune_places": tplaces["neptune"],
        "fits": fits,
        "missing_variants": [v for v in VARIANTS if v not in fits],
        "reading": reading(fits),
    }  # fmt: skip
    OUT_JSON.write_text(json.dumps(payload, indent=1, ensure_ascii=False, default=_json_default)
                        + "\n", encoding="utf-8")  # fmt: skip
    if png and fits:
        from fork.export.p1b_adams_plot import draw

        draw(payload, raw, OUT_PNG)
    n = sum(d["n_derived"] for d in fits.values()) * 2
    ledger.append("p1b_adams/summarize", "summarize", n_likelihood_calls=0, n_integrations=n,
                  sim_years=n * (1847.0 - 1781.2), wall_s=time.time() - t0,
                  cpu_s=ledger.cpu_now() - c0, n_workers=1, estimated=True,
                  note="derived places at the report epochs and the test rows")  # fmt: skip
    print(json.dumps(payload["reading"], indent=1, default=_json_default, ensure_ascii=False))
    print(f"wrote {OUT_JSON}")
    return payload


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("plan")
    p_run = sub.add_parser("run")
    p_run.add_argument("--n-workers", type=int, default=min(12, os.cpu_count() or 1))
    p_run.add_argument("--variant", action="append", choices=VARIANTS,
                       help="repeatable; variants run one after another (default: all, in "
                            "the pre-registered order)")  # fmt: skip
    p_job = sub.add_parser("job")
    p_job.add_argument("job_id")
    p_job.add_argument("--n-workers", type=int, default=4)
    p_job.add_argument("--maxiter", type=int, default=None, help="smoke run: not saved")
    p_sum = sub.add_parser("summarize")
    p_sum.add_argument("--no-png", action="store_true")
    args = ap.parse_args(argv)
    if args.cmd == "plan":
        plan()
    elif args.cmd == "run":
        for v in args.variant or VARIANTS:
            run_variant(v, args.n_workers)
    elif args.cmd == "job":
        do_job(Job.from_id(args.job_id), args.n_workers, args.maxiter)
    else:
        summarize(png=not args.no_png)


if __name__ == "__main__":
    main()
