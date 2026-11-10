"""Pre-registration timeline for the site's panel 4 (site/data/m7_timeline.json).

Every time comes from a primary source, never from the notes:
- commit time: git (committer date, the moment the text was fixed);
- push time: GitHub's repository activity API (`gh api repos/{owner}/{repo}/activity`), the
  first push whose new head is the commit or a descendant of it;
- launch, finish and first-summary times: the run drivers' logs under data/cache/ (gitignored,
  so this script runs on the machine that ran the jobs).

Run from the repo root: `uv run python docs/notes/scripts/m7_timeline.py` (needs `gh`, network).
"""

from __future__ import annotations

import datetime as dt
import json
import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
CACHE = ROOT / "data" / "cache"
OUT = ROOT / "site" / "data" / "m7_timeline.json"
PST = dt.timezone(dt.timedelta(hours=-8))

# The 2026-11-07 history rewrite (docs/notes/M6.md, "Note on commit hashes"): the M6c
# pre-registration was first pushed as ed1bf8a, which no longer exists locally.
PUSHED_AS = {"2c9ebb3": "ed1bf8a"}


def sh(*args: str) -> str:
    return subprocess.run(args, cwd=ROOT, check=True, capture_output=True, text=True).stdout


def iso(t: dt.datetime) -> str:
    return t.astimezone(PST).isoformat()


def committed(sha: str) -> str:
    return iso(dt.datetime.fromisoformat(sh("git", "log", "-1", "--format=%cI", sha).strip()))


def pushes() -> list[dict]:
    raw = sh("gh", "api", "--paginate", "repos/{owner}/{repo}/activity?per_page=100")
    events = json.loads("[" + raw.replace("][", ",")[1:-1] + "]")
    events = [e for e in events if e["activity_type"] in ("push", "force_push", "branch_creation")]
    return sorted(events, key=lambda e: e["timestamp"])


def first_push(sha: str, events: list[dict]) -> dict:
    target = PUSHED_AS.get(sha, sha)
    for e in events:
        after = e["after"]
        exact = after.startswith(target)
        known = subprocess.run(["git", "cat-file", "-e", after], cwd=ROOT).returncode == 0
        descends = (
            known
            and subprocess.run(
                ["git", "merge-base", "--is-ancestor", target, after], cwd=ROOT
            ).returncode
            == 0
        )
        if exact or descends:
            t = dt.datetime.fromisoformat(e["timestamp"].replace("Z", "+00:00"))
            return {"time": iso(t), "head": after[:7], "as": target if target != sha else None}
    return {"time": None, "head": None, "as": None}


LINE = re.compile(r"^(\d{4}-\d\d-\d\d \d\d:\d\d:\d\d) (.*)$")


def log_time(
    path: Path, pattern: str, *, last: bool = False, after: str | None = None
) -> dict | None:
    """Time of the first (or last) line matching `pattern`, optionally only after a line
    matching `after`. Log times are local (PST)."""
    if not path.exists():
        return None
    hit, armed = None, after is None
    for line in path.read_text().splitlines():
        m = LINE.match(line)
        if not m:
            continue
        if not armed and re.search(after, m.group(2)):
            armed = True
        if armed and re.search(pattern, m.group(2)):
            hit = m.group(1)
            if not last:
                break
    if hit is None:
        return None
    t = dt.datetime.strptime(hit, "%Y-%m-%d %H:%M:%S").replace(tzinfo=PST)
    return {"time": iso(t), "source": str(path.relative_to(ROOT))}


def mtime(path: Path) -> dict:
    t = dt.datetime.fromtimestamp(path.stat().st_mtime, tz=PST).replace(microsecond=0)
    return {"time": iso(t), "source": f"{path.relative_to(ROOT)} (file time)"}


def run_commit(path: Path, pattern: str) -> str | None:
    if not path.exists():
        return None
    m = re.search(pattern + r": commit ([0-9a-f]{7})", path.read_text())
    return m.group(1) if m else None


def main() -> None:
    events = pushes()
    m6c = CACHE / "nested" / "m6c"
    m8 = CACHE / "m8"

    def prereg(sha: str, role: str) -> dict:
        return {
            "sha": sha,
            "role": role,
            "committed": committed(sha),
            "pushed": first_push(sha, events),
        }

    entries = [
        {
            "id": "m6c_headline",
            "title": "Headline fit",
            "notes": "docs/notes/M6.md, M6c pre-registration",
            "preregistration": [prereg("2c9ebb3", "pre-registration")],
            "launched": log_time(m6c / "driver.log", r"^run: 4 jobs pending"),
            "finished": log_time(
                m6c / "driver.log", r"^done headline/", last=True, after=r"^run: 4 jobs pending"
            ),
            # summarize.log holds the first computed Galle-night longitude; it was not rewritten.
            "first_summary": mtime(m6c / "summarize.log"),
        },
        {
            "id": "m6c_sensitivity",
            "title": "Sensitivity runs",
            "notes": "docs/notes/M6.md, M6c sensitivity runs S1–S3: pre-registration",
            "preregistration": [prereg("2981aa1", "pre-registration")],
            "launched": log_time(m6c / "driver.log", r"^run s1_e05:"),
            "finished": log_time(m6c / "driver.log", r"^done s3_t4/", last=True),
            "first_summary": None,
        },
        {
            "id": "m6c_repeat",
            "title": "Repeatability runs",
            "notes": "docs/notes/M6.md, M6c cell-2 repeatability: pre-registration",
            "preregistration": [prereg("328350a", "pre-registration")],
            "launched": log_time(m6c / "repeat.log", r"^run: 6 jobs pending"),
            "finished": log_time(m6c / "repeat.log", r"^done ", last=True),
            "first_summary": None,
            # A 31-iteration capped smoke job, no result saved (M6.md says so in the
            # pre-registration); recorded so the earlier log line is accounted for.
            "smoke": log_time(
                CACHE / "nested" / "m6c_repeat" / "progress.log", r"\[m6c_repeat/r2_1\] start"
            ),
        },
        {
            "id": "m8_a8",
            "title": "Fast-likelihood validation",
            "notes": "docs/notes/M8.md, Pre-registration: A8b–A8f; Before launch",
            "preregistration": [
                prereg("2d2b1ff", "pre-registration"),
                prereg("01c344c", "the confirmed fallbacks and the choices the plan left open"),
            ],
            "launched": log_time(m8 / "runs" / "progress.log", r"^run: 48 jobs pending"),
            "finished": log_time(m8 / "runs" / "progress.log", r"^run finished"),
            "first_summary": log_time(m8 / "runs" / "progress.log", r"^summarize:"),
            "run_commit": run_commit(m8 / "runs" / "progress.log", r"pipeline"),
        },
        {
            "id": "m8_a8f_rerun",
            "title": "Fast-likelihood rerun",
            "notes": "docs/notes/M8.md, Follow-up pre-registration",
            "preregistration": [prereg("049e88f", "pre-registration")],
            "launched": log_time(m8 / "runs_a8f" / "progress.log", r"^run: "),
            "finished": log_time(m8 / "runs_a8f" / "progress.log", r"^run finished"),
            "first_summary": log_time(m8 / "runs_a8f" / "progress.log", r"^summarize"),
            "run_commit": run_commit(m8 / "runs_a8f" / "progress.log", r"pipeline a8f"),
        },
        {
            "id": "m8b",
            "title": "Higher-resolution reruns (decide the 95 % range)",
            "notes": "docs/notes/M8.md, Follow-up pre-registration",
            "preregistration": [prereg("049e88f", "pre-registration")],
            "launched": log_time(m8 / "runs_m8b" / "progress.log", r"^run: "),
            "finished": log_time(m8 / "runs_m8b" / "progress.log", r"^run finished"),
            "first_summary": log_time(m8 / "runs_m8b" / "progress.log", r"^summarize"),
            "run_commit": run_commit(m8 / "runs_m8b" / "progress.log", r"pipeline m8b"),
        },
        {
            "id": "m8_arbiter",
            "title": "Full-likelihood check of the eccentricity test",
            "notes": "docs/notes/M8.md, Follow-up pre-registration",
            "preregistration": [prereg("049e88f", "pre-registration")],
            "launched": log_time(CACHE / "nested" / "m8_arbiter" / "progress.log", r"start"),
            "finished": log_time(CACHE / "nested" / "m8_arbiter" / "progress.log", r"run finished"),
            "first_summary": None,
        },
    ]
    for e in entries:
        e["status"] = (
            "not launched"
            if e["launched"] is None
            else "running"
            if e["finished"] is None
            else "done"
        )

    out = {
        "description": "Pre-registration timeline: when each pre-registration was committed and "
        "first pushed (public), and when its runs launched, finished and were first summarized. "
        "Times are PST (UTC-8). Built by docs/notes/scripts/m7_timeline.py from git, GitHub's "
        "repository activity API and the run logs.",
        "generated": iso(dt.datetime.now(tz=PST).replace(microsecond=0)),
        "sources": {
            "committed": "git committer date",
            "pushed": "GitHub repository activity API: first push whose head is the commit or "
            "a descendant",
            "launched/finished/first_summary": "run driver logs under data/cache/ "
            "(named per entry)",
        },
        "rewritten": {
            "2c9ebb3": {"pushed_as": "ed1bf8a", "source": "docs/notes/M6.md, Note on commit hashes"}
        },
        "entries": entries,
        # The spec's M6c fits not run in v1.0.0 (docs/release_checklist.md, "Deferred to v1.1",
        # moved there 2026-11-09). Not pre-registered yet, so they have no timeline entry.
        "deferred": {
            "source": "docs/release_checklist.md, Deferred to v1.1 (2026-11-09)",
            "to_version": "1.1",
            "fitted": "Le Verrier's data choice on the modern reduction track",
            "m6c_fits": [
                {"key": "bouvard", "name": "Bouvard", "kind": "variant"},
                {"key": "adams", "name": "Adams", "kind": "variant"},
                {"key": "track_1846", "name": "1846 reduction track", "kind": "track"},
                {"key": "lv_shape", "name": "Le Verrier shape model", "kind": "model"},
            ],
            "status": "not yet pre-registered; to be planned and pushed before launch",
            "opens": {"adams": "the Adams–Challis study"},
        },
        # The pre-registered reading rules of the M6c sensitivity runs, as written in
        # docs/notes/M6.md before any of their sampling (commit 2981aa1). The site applies them
        # to the numbers in m6_fits.json; it does not restate a verdict.
        "sensitivity_reading_rules": {
            "source": "docs/notes/M6.md, pre-registration of the sensitivity runs (M6c, commit 2981aa1)",
            "s1_e05": {
                "quantity": "p_a_below_25_au",
                "moved_above": 0.5,
                "not_at_all_below": 0.05,
                "text": "the evidence has moved to the low-a peak if P(a < 25 AU) > 0.5. It has "
                "moved partly if P(a < 25 AU) is 0.05–0.5, and not at all if it is < 0.05.",
            },
            "s2_sigma_data": {
                "jpl_keeps": "in_68",
                "max_shift_headline_half_widths": 0.5,
                "text": "the result survives if JPL's 325.873° keeps its headline status (inside "
                "the 68 % interval) and the Galle-night median moves by < 0.5 headline "
                "68 % half-widths.",
            },
            "s3_t4": {
                "jpl_keeps": "in_68",
                "max_shift_headline_half_widths": 0.5,
                "text": "heavier tails change the headline if JPL's status in the 68 % interval "
                "changes or the median moves by ≥ 0.5 headline half-widths.",
            },
        },
    }
    OUT.write_text(json.dumps(out, indent=1, ensure_ascii=False) + "\n")
    for e in entries:
        p = e["preregistration"]
        print(
            e["id"],
            [(x["sha"], x["committed"][11:19], (x["pushed"]["time"] or "-")[11:19]) for x in p],
            "launch",
            (e["launched"] or {}).get("time", "-")[11:19],
            "finish",
            (e["finished"] or {}).get("time", "-")[11:19],
            "summary",
            (e["first_summary"] or {}).get("time", "-")[11:19],
            e["status"],
        )


if __name__ == "__main__":
    main()
