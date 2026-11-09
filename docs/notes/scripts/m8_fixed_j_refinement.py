"""After-the-fact diagnostic (not a spec mode): mode (b) with J₀ reused, on the A8a points.

    uv run python docs/notes/scripts/m8_fixed_j_refinement.py

A8a (docs/notes/M8.md) showed that mode (a)'s error is all in χ²: its errors in χ² and in
−2 ln L agree to 1e-5, so the fixed J₀ is not what fails. This evaluates "(b′)": mode (b)'s
second integration at u₁ = u*(a), but with the second Gauss–Newton step taken with J₀ instead
of a forward-difference J₁. That is 2 integrations per call instead of 6. It is measured here
against the same references as A8a (the stored full profile χ²_min and exact-Hessian Laplace),
and timed. It is a proposal for the A8b–A8f pre-registration, not the chosen mode.

Reads data/cache/m8/a8a/*.json; writes docs/notes/M8_fixed_j_refinement.json.
"""

import json
import math
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

from fork.export import m8
from fork.model import linearized as lz

OUT = Path(__file__).resolve().parents[1] / "M8_fixed_j_refinement.json"
_P: dict = {}


def _init() -> None:
    _P["headline"] = m8.problem("headline")
    _P["s1_e05"] = m8.problem("s1_e05")


def fixed_j_refined(prob: lz.LinearizedProblem, phi: np.ndarray) -> lz.GaussianMarginal:
    z = prob.z_fn(phi)
    m = lz.linear_marginal(z(prob.u_ref), prob.jac_ref, prob.u_ref, prob.log_volume)
    return lz.linear_marginal(z(m.u_star), prob.jac_ref, m.u_star, prob.log_volume)


def _eval(rec: dict) -> dict:
    _, _, _, prob = _P["s1_e05" if rec["source"] == "s1" else "headline"]
    m = fixed_j_refined(prob, np.array(rec["q"])[4:])
    return {"source": rec["source"], "class": rec["class"], "chi2_full_min": rec["chi2_full_min"],
            "d_chi2": m.chi2_min - rec["chi2_full_min"],
            "d_m2lnl": m.m2lnl - rec["m2lnl_laplace"]}  # fmt: skip


def main() -> None:
    recs = m8.load_a8a()
    best = m8.best_m2lnl()
    with ProcessPoolExecutor(8, initializer=_init) as ex:
        out = list(ex.map(_eval, recs, chunksize=25))
    summary = {"description": __doc__.strip().splitlines()[0], "best_m2lnl": round(best, 3),
               "tolerances": {"abs_d_chi2": m8.TOL_CHI2, "abs_d_m2lnl": m8.TOL_M2LNL},
               "groups": {}}  # fmt: skip
    for name in ("core", "ridge", "other_cells", "prior_tail"):

        def keep(r, name=name):
            if name == "prior_tail":
                return r["source"].startswith("prior")
            return r["source"] in ("r2", "s1") and r["class"] == name

        rs = [r for r in out if keep(r)]
        dom = [r for r in rs if (r["chi2_full_min"] - best) / 2 <= m8.DOMAIN_DLNL]
        for label, sub in (("in_domain", dom), ("all", rs)):
            if not sub:
                summary["groups"][f"{name}_{label}"] = {"n": 0}
                continue
            dc = np.array([r["d_chi2"] for r in sub])
            dm = np.array([r["d_m2lnl"] for r in sub])
            summary["groups"][f"{name}_{label}"] = {
                "n": len(sub), "d_chi2": m8._dist(dc), "d_m2lnl": m8._dist(dm),
                "n_over_chi2_tol": int(np.sum(np.abs(dc) >= m8.TOL_CHI2)),
                "n_over_m2lnl_tol": int(np.sum(np.abs(dm) >= m8.TOL_M2LNL))}  # fmt: skip
    spec = [r for r in out if r["source"] in ("r2", "s1")
            and (r["chi2_full_min"] - best) / 2 <= m8.DOMAIN_DLNL]  # fmt: skip
    dc = np.abs([r["d_chi2"] for r in spec])
    dm = np.abs([r["d_m2lnl"] for r in spec])
    summary["spec_points"] = {"n": len(spec), "max_abs_d_chi2": round(float(dc.max()), 6),
                              "max_abs_d_m2lnl": round(float(dm.max()), 6),
                              "would_pass": bool(dc.max() < m8.TOL_CHI2
                                                 and dm.max() < m8.TOL_M2LNL)}  # fmt: skip

    # timing, one thread, the A8e points
    _init()
    _, _, _, prob = _P["headline"]
    qs = [np.array(r["q"]) for r in recs if r["source"] == "r2"][:100]
    fixed_j_refined(prob, qs[0][4:])
    best_ms = math.inf
    for _ in range(3):
        t0 = time.process_time()
        for q in qs:
            fixed_j_refined(prob, q[4:])
        best_ms = min(best_ms, (time.process_time() - t0) / len(qs) * 1e3)
    summary["ms_per_call"] = round(best_ms, 3)
    summary["integrations_per_call"] = 2
    OUT.write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps({k: summary[k] for k in ("spec_points", "ms_per_call")}, indent=1))
    for k, v in summary["groups"].items():
        dc, dm = v.get("d_chi2", {}), v.get("d_m2lnl", {})
        print(k, v.get("n"), dc.get("abs_max"), dm.get("abs_max"), v.get("n_over_chi2_tol"),
              v.get("n_over_m2lnl_tol"))  # fmt: skip


if __name__ == "__main__":
    main()
