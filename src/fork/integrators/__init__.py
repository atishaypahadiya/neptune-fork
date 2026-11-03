"""N-body integrators. REBOUND is used only in tests to validate these (CLAUDE.md rule 1)."""

# Project-wide default step [day]. 14 d keeps Wisdom–Holman's 1000-yr energy error < 1e-8
# (measured 9.93e-9 with the default composition) and its Uranus error ~0.008″ — see
# docs/notes/M3.md. 65 yr of Sun + 3 planets costs ~1.6 ms.
DEFAULT_DT_DAYS: float = 14.0
