"""Diagnosis of M8's gate failures (docs/notes/M8.md, "Gate diagnosis"): a hand check, outside
``m8_runs.summarize``, of gate 1 and of the cell-2 ridge counts.

    uv run python docs/notes/scripts/m8_gate_diagnosis.py     # ~1 min, reads the run cache only

1. Gate 1 by hand, with plain numpy (not ``fork.inference.reweight``), for each fit:
   - each cell's Δ = ln L_b − ln L_fast, its ESS_k and ln Z_k = ln Z_fast + ln mean e^Δ;
   - the combined weights ln c_k + ln Z_fast,k + Δ_i − ln N_k, as in ``reweight.combine``;
   - the combined ESS/N, and the identity ESS = 1 / Σ_k (p_k² / ESS_k), p_k the cell shares;
   - the same at Δ ≡ 0, the best any fast likelihood could do with 6000 draws per cell;
   - the ratio of the two, the part of ESS/N due to reweighting.
   ln Z_fast is checked against the dynesty result itself.
2. The same Δ ≡ 0 bound for M6c's published cell ln Z_k (full likelihood, no reweighting).
3. Cell-2 ridge branches per run: counts among the 6000 draws, against the raw dynesty
   importance-weighted mass (no resampling), so a resampling bug would show as a mismatch.

Writes nothing.
"""

from __future__ import annotations

import math

import numpy as np
from scipy.special import logsumexp

from fork.export import m8_runs as R

M6C_PUBLISHED_LN_Z = np.array([-364.87, -365.18, -359.77, -364.87])  # M8.md, A8c row


def gate1(v: str, s: int) -> dict:
    lc, lz, ess_k, n_k, logw = [], [], [], [], []
    for k in range(R.n_cells(v)):
        d = R.load_reweighted(v, k, s)
        delta = d["lnl_b"] - d["lnl_fast"]
        w = np.exp(delta - delta.max())
        ess_k.append(w.sum() ** 2 / (w**2).sum())
        n_k.append(delta.size)
        c = math.log(R._cell_space(v, k).slice_weight)
        lc.append(c)
        lz.append(float(d["ln_z_fast"]) + float(logsumexp(delta)) - math.log(delta.size))
        logw.append(c + float(d["ln_z_fast"]) + delta - math.log(delta.size))
        if s == R.MERGED:
            dyn = R.load_posterior(v, k, s).logz
            assert abs(dyn - float(d["ln_z_fast"])) < 1e-9, (v, k, dyn)
    lw = np.concatenate(logw)
    W = np.exp(lw - lw.max())
    n = lw.size
    p = np.exp(np.array(lc) + np.array(lz) - logsumexp(np.array(lc) + np.array(lz)))
    ess = W.sum() ** 2 / (W**2).sum()
    ess_identity = 1.0 / np.sum(p**2 / np.array(ess_k))
    ess_delta0 = 1.0 / np.sum(p**2 / np.array(n_k))
    return {"ess_frac": ess / n, "identity": ess_identity / n, "delta0": ess_delta0 / n,
            "reweighting_part": ess / ess_delta0, "shares": p}  # fmt: skip


def main() -> None:
    print("gate 1: fit | ESS/N | identity | at Δ ≡ 0 | ESS/ESS(Δ≡0) | cell shares")
    for v in R.VARIANTS:
        for s in (R.MERGED, *R.SEEDS):
            g = gate1(v, s)
            name = f"{v}/{'merged' if s == R.MERGED else f's{s}'}"
            print(f"  {name:22s} {g['ess_frac']:.4f} {g['identity']:.4f} {g['delta0']:.4f} "
                  f"{g['reweighting_part']:.4f}  {np.round(g['shares'], 4)}")  # fmt: skip

    lc = np.array([math.log(R._cell_space("headline", k).slice_weight) for k in range(4)])
    p6 = np.exp(lc + M6C_PUBLISHED_LN_Z - logsumexp(lc + M6C_PUBLISHED_LN_Z))
    print(f"M6c published shares {np.round(p6, 4)}: Δ ≡ 0 bound {1 / (4 * np.sum(p6**2)):.4f}")

    print("cell 2 ridge: run | draws inner, outer (fraction) | raw dynesty mass inner, outer")
    for v in R.VARIANTS:
        space = R.m8.problem(v)[2]
        for s in (*R.SEEDS, R.MERGED):
            res = R.load_posterior(v, 2, s)
            pr = R._perturber(space, res.samples)
            wr = res.weights / res.weights.sum()
            pd = R._perturber(space, R.load_reweighted(v, 2, s)["phi"])
            out = []
            for pp, ww in ((pd, np.full(len(pd), 1.0)), (pr, wr)):
                a, e = pp[:, 1], pp[:, 2]
                out.append((ww[(e >= R.E_RIDGE) & (a < R.INNER_OUTER_A)].sum(),
                            ww[(e >= R.E_RIDGE) & (a >= R.INNER_OUTER_A)].sum()))  # fmt: skip
            (ni, no), (ri, ro) = out
            print(f"  {R.pid(v, 2, s):26s} {int(ni):4d} {int(no):4d} ({(ni + no) / len(pd):.3f})"
                  f" | {ri:.4f} {ro:.4f} ({ri + ro:.3f})")  # fmt: skip


if __name__ == "__main__":
    main()
