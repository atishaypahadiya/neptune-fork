"""Noise-free M5b posterior (free a, e ≤ 0.3): e limits as `m5b summarize` computes them, and
how the log m / ln a / 1846-longitude offsets line up with the mass–distance degeneracy
(correlations, principal axis, Mahalanobis distance). Numbers quoted in docs/notes/M5.md.

    uv run python docs/notes/scripts/m5b_noise_free_degeneracy.py
"""

import math

import numpy as np

from fork.export import m5b
from fork.export.to_json import A1_EPOCH_JD_TDB
from fork.inference.nested import ERR_INFLATION, NestedResult, PartitionedResult
from fork.inference.profile import GalleGeometry, galle_lon_of_date_rad
from fork.model.forward import ForwardContext

run = "free_a_e0.3"
jobs = [j for j in m5b.jobs_for(None) if j.run == run and j.rep == 1]
res = [NestedResult.load(j.path(".pkl")) for j in jobs]
rng = m5b.summary_rng(None, run)
p = PartitionedResult(tuple(res), err_scale=ERR_INFLATION)
eq = p.equal_weight(rng)
space = res[0].space
th = np.array([space.theta(v) for v in eq])
tp = m5b.truth_params()
e = m5b._e_summary(th[:, 8], tp["e"], m5b.E_MAX[run])
print("e:", e)
idx = rng.choice(len(eq), size=min(m5b.DERIVED_MAX_SAMPLES, len(eq)), replace=False)
ctx, geo = ForwardContext.from_horizons(A1_EPOCH_JD_TDB), GalleGeometry.from_horizons()
lon = np.array([math.degrees(galle_lon_of_date_rad(space.theta(v), ctx, geo)) for v in eq[idx]])
lon = (lon - tp["galle_lon_ecl_of_date_deg"] + 180) % 360 - 180 + tp["galle_lon_ecl_of_date_deg"]
X = np.column_stack([th[idx, 6], np.log(th[idx, 7]), lon])  # log10 m, ln a, lon
names = ["log10 m", "ln a", "lon"]
C = np.corrcoef(X.T)
print("n", len(X), "corr\n", np.round(C, 3))
t = np.array([math.log10(tp["mass_earth"]), math.log(tp["a_au"]), tp["galle_lon_ecl_of_date_deg"]])
mu, S = X.mean(0), np.cov(X.T)
d = t - mu
z1 = d / np.sqrt(np.diag(S))
maha = float(math.sqrt(d @ np.linalg.solve(S, d)))
print("1-D z (truth − mean)/sd:", np.round(z1, 3),
      " in quadrature:", round(float(np.linalg.norm(z1)), 3))  # fmt: skip
print("3-D Mahalanobis:", round(maha, 3))
sd = np.sqrt(np.diag(S))
w, V = np.linalg.eigh(S / np.outer(sd, sd))
print("eig of corr:", np.round(w, 3))
print("leading eigvec:", np.round(V[:, -1], 3))
dz = z1 / np.linalg.norm(z1)
print("cos(offset, leading eigvec):", round(float(abs(dz @ V[:, -1])), 3))
d2 = d[:2]  # (log m, ln a) alone
print("2-D (m, a) Mahalanobis:", round(float(math.sqrt(d2 @ np.linalg.solve(S[:2, :2], d2))), 3))
