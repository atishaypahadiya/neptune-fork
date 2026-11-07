# Project rules

These rules apply to all code in this repository. Code, tests and notes cite them as
"docs/RULES.md rule N".

1. **Own integrator in `src/fork/integrators/`.** REBOUND is a *test dependency only*, used to validate.
   Never import rebound from `src/`.
2. **Uranus's in-plane elements (a, e, ϖ, λ₀) are always free parameters** in any fit. Its i and Ω
   are held at their epoch (1781) values: longitude-only data can't constrain them (a second-order,
   ~10″ effect) and, freed, they absorb the perturber signal. The leak this introduces is pinned
   by a test (Neptune's effect on i, Ω over 1781–1846 is 0.017″ in longitude). Never hold Uranus's
   modern Horizons state (or any in-plane element) fixed when Neptune is removed — that leaks the
   answer.
3. **Synthetic before historical.** No code touches historical observations until test A3
   (recover a planted planet) passes.
4. **Frames are explicit.** Every angle carries its frame in the variable name or type
   (`lon_ecl_j2000_rad`, `lon_ecl_of_date_rad`). Precession (~50″/yr) dwarfs the ~130″ signal;
   a frame bug will fake a planet.
5. **Units:** AU, days, solar masses. G = k², k = 0.01720209895 (Gaussian). Angles internally in
   radians; arcseconds only at I/O and plotting boundaries.
6. **Wrap-safe residuals:** always `((obs - model + π) % 2π) - π`.
7. **No paid services.** Everything runs locally or on free GitHub Actions / Pages.
   API credits are reserved for Phase 3; do not add LLM calls in Phase 1.
