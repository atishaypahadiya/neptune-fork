# Test fixtures

`horizons/`: the JPL Horizons states the default test job (`-m "not slow and not network"`)
reads, and nothing else. `tests/conftest.py` points `FORK_CACHE_DIR` at a temporary copy for
the whole session and makes any live fetch raise, so that job makes no network calls.

- **Source:** NASA/JPL Horizons API (https://ssd.jpl.nasa.gov/api/horizons.api), DE441, via
  astroquery 0.4.11 (`fork.data.horizons.fetch_from_horizons`).
- **Fetched:** into the local cache (`data/cache/`) between 2026-11-03 and 2026-11-08 (the cache
  files' modification times; the cache does not record when each epoch was fetched). The fixture
  was cut from that cache on 2026-11-09.
- **Query:** `EPHEM_TYPE=VECTORS`, `CENTER='500@10'` (the Sun's centre), `REF_PLANE=ECLIPTIC`
  (ecliptic and mean equinox of J2000), `VEC_CORR=NONE` (geometric), `OUT_UNITS=AU-D`, epochs
  as a `TLIST` of TDB Julian Dates. Bodies by Horizons ID: Jupiter 5, Saturn 6, Uranus 7,
  Neptune 8, Earth 399, Earth–Moon barycentre 3, Mercury 199, Venus 299, Mars 4.
- **Format:** one `horizons_<body>.npz` per body, the cache's own format (`jd_tdb`,
  `state_ecl_j2000` in AU and AU/day, `horizons_id`, `center`, `refplane`).
- **Regenerate:** `uv run python tests/fixtures/build_horizons_fixture.py`. It runs the default
  job in record mode (states served from `data/cache/`, fetched live only for an epoch the local
  cache lacks) and writes exactly the epochs the tests ask for. It then reruns the job from the
  new fixture with live fetches disabled, to check the fixture is complete.

Neptune is in the fixture (2 epochs) because tests compare with it. The fit path never reads it:
A6d rule 5 (`tests/test_a6d.py`, `tests/test_m6c.py`) checks that no `fork` code on the fit
path names this folder.

`tests/data/horizons_m4.json` is the older, smaller fixture used by the `offline_horizons`
modules; it is unchanged.
