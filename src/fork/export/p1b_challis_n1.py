"""Phase 1b step 4, the chart case C2 under the named input N1 itself.

The frozen plan (docs/notes/P1B_challis_detection.md, § 5, C2) runs the Berlin Hora XXII case
under N1's two fallbacks "and under N1 itself if obtained". N1 was read on 12 November 2026
from the sheet (a published photograph of the Potsdam copy; the note's dated section gives the
source table). This module runs that third version on the same accepted draws, from the plan's
seed, and writes docs/notes/P1B_challis_detection_n1.json. The step 4 code, cache and JSON are
read and not changed.

    uv run python -m fork.export.p1b_challis_n1 run        # about 3 min, one process
"""

from __future__ import annotations

import argparse
import json
import math
import time

import numpy as np

from fork.export import p1b_challis as S
from fork.historical import challis_1846 as C

OUT_JSON = S.ROOT / "docs" / "notes" / "P1B_challis_detection_n1.json"
CACHE = S.CACHE / "n1"

# N1 as read from the sheet. Right ascension limits are the outermost lines of its grid, one
# 4-minute column beyond the lines lettered XXII and XXIII; the engraved title gives the
# declination range and the epoch ("15° südlicher bis 15° nördlicher Abweichung, für 1800").
N1 = {
    "sheet": "Akademische Sternkarten, Hora XXII (Argelander, Åbo; engraved by K. Kolbe)",
    "ra_limits_h": [21.0 + 56.0 / 60.0, 23.0 + 4.0 / 60.0],
    "ra_limits_text": ["21h 56m", "23h 04m"],
    "dec_limits_deg": [-15.0, 15.0],
    "equinox": 1800,
    "faintest_class_in_key": "9-10",
    "read_from": "the photograph of the Potsdam copy in Krajnović 2021, Fig. 6.10 (lower part "
    "of the sheet: declination −8° to −15°, the key and the title)",
    "date_read": "2026-11-12",
}


def geometry_n1(pre: dict) -> dict:
    """The step 4 geometry with the chart test of N1 itself in the slot C2_1800 reads."""
    geo = S.geometry(pre)
    n = pre["nights"]
    ra1800 = np.array([x["neptune"]["ra1800_h"] for x in n])
    dec1800 = np.array([x["neptune"]["dec1800_deg"] for x in n])
    lo, hi = N1["ra_limits_h"]
    dlo, dhi = N1["dec_limits_deg"]
    on = (ra1800 >= lo) & (ra1800 <= hi) & (dec1800 >= dlo) & (dec1800 <= dhi)
    geo = dict(geo)
    geo["on_chart_n1"] = on
    return geo


def run(n_accept: int = S.N_ACCEPT) -> dict:
    pre = json.loads(S.PRERUN_JSON.read_text())
    geo = geometry_n1(pre)
    same_as_first_fallback = bool(np.array_equal(geo["on_chart_n1"], geo["on_chart_1800"]))
    # C.groups reads the chart test of the case C2_1800 from geo["on_chart_1800"]; under N1
    # that slot holds N1's test. Everything else of the geometry is the step 4 geometry.
    geo_case = dict(geo)
    geo_case["on_chart_1800"] = geo["on_chart_n1"]
    CACHE.mkdir(parents=True, exist_ok=True)
    t0, c0 = time.time(), time.process_time()
    out = {}
    for tag, priors, n in (
        ("main", None, n_accept),
        ("halved", C.halved(C.PRIORS), min(S.N_ACCEPT_HALVED, n_accept)),
    ):
        acc, _, counts = S.accepted_draws(geo, n, priors)
        S._log(f"n1 {tag}: {n} accepted of {counts['n_drawn']} at {time.time() - t0:.0f} s")
        t1 = np.zeros((2, n, 2))
        t2 = np.zeros((2, n, 2))
        a0 = np.zeros(n)
        for k in range(math.ceil(n / S.CHUNK)):
            f = CACHE / f"{tag}_chunk_{k:03d}.npz"
            sl = slice(k * S.CHUNK, min((k + 1) * S.CHUNK, n))
            if not f.exists():
                p1 = np.zeros((2, sl.stop - sl.start, 2))
                p2 = np.zeros_like(p1)
                pa = np.zeros(sl.stop - sl.start)
                for i, d in enumerate(acc[sl]):
                    for ri, name in enumerate(S.READINGS):
                        f1, f2 = C.cdfs(d, C.groups(d, geo_case, "C2_1800", held=name == "held"))
                        p1[ri, i] = f1[C.PRIMARY_LAST], f1[C.SECONDARY_LAST]
                        p2[ri, i] = f2[C.PRIMARY_LAST], f2[C.SECONDARY_LAST]
                    fa, _ = C.cdfs(d, C.groups(d, geo, "A0", held=True))
                    pa[i] = fa[C.PRIMARY_LAST]
                tmp = CACHE / f"{tag}_chunk_{k:03d}.tmp.npz"
                np.savez_compressed(tmp, t1=p1, t2=p2, a0=pa)
                tmp.rename(f)
                S._log(f"n1 {tag}: chunk {k + 1} at {time.time() - t0:.0f} s")
            z = np.load(f)
            t1[:, sl], t2[:, sl], a0[sl] = z["t1"], z["t2"], z["a0"]
        out[tag] = {"t1": t1, "t2": t2, "a0": a0, "acceptance": counts}
    m, h = out["main"], out["halved"]
    # the same draws as step 4: its cached values of the first fallback and of A0, draw by draw
    step4 = S._load("main")
    ci, c0i = C.CASES.index("C2_1800"), C.CASES.index("A0")
    checks = {
        "n1_chart_test_equals_first_fallback_on_every_night": same_as_first_fallback,
        "acceptance_equals_step_4": m["acceptance"] == step4["acceptance"],
        "max_abs_A0_minus_step_4_A0": float(np.abs(m["a0"] - step4["t1"][0, c0i, :, 0]).max()),
        "max_abs_n1_minus_step_4_first_fallback": float(
            np.abs(m["t1"][:, :, 0] - step4["t1"][:, ci, :, 0]).max()
        ),
    }
    entry = {
        "label": "C2 chart, N1 itself (the sheet: limits for the equinox of 1800)",
        "conditional_on": "N1 as read from the sheet on 12 November 2026 (no fallback)",
    }
    for ri, name in enumerate(S.READINGS):
        entry[name] = {
            "primary": S._stats(m["t1"][ri, :, 0]),
            "by_1_october": S._stats(m["t1"][ri, :, 1]),
            "t2_primary": S._stats(m["t2"][ri, :, 0]),
            "t2_by_1_october": S._stats(m["t2"][ri, :, 1]),
        }
    entry["held"]["paired_minus_A0"] = S._stats(m["t1"][0, :, 0] - m["a0"])
    s = entry["held"]["primary"]
    row = C.reading(s["q05"], s["q95"])
    hs = S._stats(h["t1"][0, :, 0])
    entry["reading_row"] = row
    entry["halved_priors"] = {"primary": hs, "reading_row": C.reading(hs["q05"], hs["q95"])}
    entry["prior_dependent"] = entry["halved_priors"]["reading_row"] != row
    pd_ = entry["held"]["paired_minus_A0"]
    entry["paired_reading"] = C.paired_reading(pd_["q05"], pd_["q95"])
    committed = json.loads((S.ROOT / "docs" / "notes" / "P1B_challis_detection.json").read_text())
    rows = {v: committed["cases"][v]["reading_row"] for v in ("C2_1800", "C2_date")}
    rows["C2_N1"] = row
    nights = pre["nights"]
    result = {
        "plan": "docs/notes/P1B_challis_detection.md, § 5 C2: 'and under N1 itself if obtained'",
        "step_4_json_sha256": S.sha256(S.ROOT / "docs" / "notes" / "P1B_challis_detection.json"),
        "prerun_sha256": S.sha256(S.PRERUN_JSON),
        "n1": N1,
        "neptune_ra_1800_h": {
            x["date"]: x["neptune"]["ra1800_h"]
            for x in nights
            if x["date"] in ("1846-07-29", "1846-07-30", "1846-08-04", "1846-08-12", "1846-09-29")
        },
        "n_accepted": int(m["t1"].shape[1]),
        "seeds": {"draws": list(S.SEED_DRAWS)},
        "closed_form_only": "the event-sampling check of the closed form was made in step 4 on "
        "the first fallback, whose per-draw values these equal",
        "checks": checks,
        "case": entry,
        "reading_rows_by_version_of_n1": rows,
        "c2_combined_reading_row_in_step_4": committed["c2_combined_reading_row"],
        "compute": {"wall_s": time.time() - t0, "cpu_s": time.process_time() - c0},
    }
    S._dump(OUT_JSON, result)
    S._log(f"n1 done: wall {time.time() - t0:.0f} s; wrote {OUT_JSON}")
    return result


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("run")
    ap.parse_args(argv)
    r = run()
    print(
        json.dumps(
            {
                "checks": r["checks"],
                "held": r["case"]["held"]["primary"],
                "row": r["case"]["reading_row"],
            },
            indent=1,
        )
    )


if __name__ == "__main__":
    main()
