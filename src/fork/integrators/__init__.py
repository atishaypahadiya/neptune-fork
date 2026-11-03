"""N-body integrators. REBOUND is used only in tests to validate these (project rule 1)."""

# Project-wide default step [day]. 12 d keeps Wisdom–Holman's 1000-yr energy error < 1e-8 with
# margin (14 d sat 0.7% under it) — see docs/notes/M3.md for the measured numbers.
DEFAULT_DT_DAYS: float = 12.0
