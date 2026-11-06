"""A6c runner: synthetic data built from Table 77's layout only, reproducible per seed."""

import numpy as np

from fork.export import a6c
from fork.model import observation as om


def test_jobs_cover_ten_seeds_and_four_cells():
    jobs = a6c.all_jobs()
    assert len(jobs) == 40 and {j.seed for j in jobs} == set(range(10))
    assert all(a6c.Job.from_id(j.id) == j for j in jobs)
    assert [j.nlive for j in jobs[:4]] == [2000, 1000, 1000, 1000]


def test_dataset_reads_no_observed_value_and_is_reproducible():
    rows = om.read_rows()
    blank = [dict(r, value="", value_as_printed="", lv_residual="") for r in rows]
    obs_a, *_, box_a, off_a = a6c.dataset(3, rows)
    obs_b, *_, box_b, off_b = a6c.dataset(3, blank)
    assert np.array_equal(obs_a.value_rad, obs_b.value_rad)  # same seed, same machine
    assert np.array_equal(off_a, off_b) and box_a == box_b
    obs_c, *_ = a6c.dataset(4, rows)
    assert not np.allclose(obs_a.value_rad, obs_c.value_rad)
    assert obs_a.n_measurements == 529


def test_coverage_counts_and_flags_the_range():
    def seed(k, hit):
        x = {"median": 325.0, "minus": 1.0, "plus": 1.0, "width68": 2.0,
             "median_minus_truth": 0.5 if hit else 3.0, "truth_in_68": hit}  # fmt: skip
        return {"seed": k, "galle_lon_ecl_of_date_deg": x}

    c = a6c.coverage([seed(k, k < 7) for k in range(10)])
    assert c["n_truth_in_68"] == 7 and c["pass_5_to_8_of_10"]
    c = a6c.coverage([seed(k, k < 9) for k in range(10)])
    assert not c["pass_5_to_8_of_10"]
