"""examples/toy_fit.py: the demo runs end to end with a capped iteration count (CI-fast)."""

import importlib.util
import math
from pathlib import Path

TOY = Path(__file__).resolve().parents[1] / "examples" / "toy_fit.py"


def _load_toy():
    spec = importlib.util.spec_from_file_location("toy_fit", TOY)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_toy_fit_capped_run(offline_horizons):
    """A 15-iteration run on the demo's cell: the summary is finite, the posterior draws lie in
    the cell's a and e range, and the report is labelled as a demo, not a result."""
    toy = _load_toy()
    s = toy.run(0, nlive=12, maxiter=15, n_draws=5)
    assert s["niter"] >= 15 and s["ncall"] > 0
    assert math.isfinite(s["logz_k"]) and s["logzerr_k"] > 0
    a_lo, a_hi = s["cell"]["a_au"]
    e_lo, e_hi = s["cell"]["e"]
    assert a_lo <= s["a_au_median"] <= a_hi and e_lo <= s["e_median"] <= e_hi
    lon = s["lon_1846_deg"]
    assert all(0.0 <= lon[k] < 360.0 for k in ("median", "lo68", "hi68"))
    assert abs(s["lon_1846_true_deg"] - 325.873) < 0.01
    text = toy.report(s)
    assert text.startswith(toy.BANNER) and text.rstrip().endswith("docs/notes/M5.md.")
