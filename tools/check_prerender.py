"""Acceptance check for the prerendered site (tools/prerender_site.py's output).

For each page of the site (the main paper and its companion, site/adams-challis/):

1. Loads the built index.html with JavaScript disabled (what a crawler, link preview, archive
   or reader mode gets) and runs the number audit on it: every number in the text and in the
   SVG labels must be a value in the site's data files, a number from the embedded registry of
   derived numbers (``drv``), or a section number. 0 unmatched.
2. Loads the same page with JavaScript (it redraws everything from site/data/) and runs the
   same audit. The two pages must match exactly: the same numbers in the same order, the same
   text and the same SVG labels (the JS-only controls, ``data-js``, left out of both).
3. Head tags: title = the paper's title (h1), description = the abstract's first sentence,
   canonical URL, Open Graph and Twitter card tags, and plate-i.png at 1200 x 630.

    uv run --group site python tools/check_prerender.py _site [--base-url URL]
"""

from __future__ import annotations

import argparse
import json
import struct
import sys
from pathlib import Path

from playwright.sync_api import sync_playwright

sys.path.insert(0, str(Path(__file__).resolve().parent))
from prerender_site import (
    BASE_URL,
    OG_SIZE,
    PAGES,
    STRUCTURE,
    VIEWPORT,
    launch,
    open_rendered,
    serve,
)

# Runs in either page (JS on or off): text, numbers and the audit. The audit is the site's
# acceptance rule (docs/DESIGN.md): each number is a JSON value, a derived number or structure.
AUDIT = r"""({json, derived, structure}) => {
  const NUM = /(?<![\p{L}\d.])\d+(?:\.\d+)?(?![\p{L}\d])/gu;
  const STRUCTURE = new Set(structure);
  const norm = (t) => {
    if (!t.includes(".")) return String(Number(t));
    const [i, f] = t.split(".");
    return `${Number(i)}.${f}`;
  };
  const variants = (v, set) => {
    if (typeof v !== "number" || !isFinite(v)) return;
    const a = Math.abs(v);
    for (let k = 0; k <= 4; k++) set.add(norm(a.toFixed(k)));
    set.add(norm(String(a)));
  };
  const walk = (x, set) => {
    if (typeof x === "number") variants(x, set);
    else if (typeof x === "string") for (const m of x.matchAll(NUM)) set.add(norm(m[0]));
    else if (Array.isArray(x)) x.forEach((y) => walk(y, set));
    else if (x && typeof x === "object") for (const [k, y] of Object.entries(x)) {
      for (const m of k.matchAll(/\d+/g)) set.add(norm(m[0]));
      walk(y, set);
    }
  };
  const allow = new Set(), der = new Set();
  json.forEach((v) => walk(v, allow));
  derived.forEach((x) => variants(x, der));

  // The JS-only controls are hidden in the static page; leave them out of both.
  const js = [...document.querySelectorAll("[data-js]")], was = js.map((e) => e.hidden);
  js.forEach((e) => { e.hidden = true; });
  const text = document.body.innerText.replace(/[ \t]+/g, " ").replace(/\n\s*\n+/g, "\n").trim();
  js.forEach((e, i) => { e.hidden = was[i]; });
  const svg = [...document.querySelectorAll("svg text")].map((t) => t.textContent);
  const all = [text, ...svg].join("\n");

  const numbers = [], unmatched = [], counts = { json: 0, derived: 0, structure: 0 };
  for (const m of all.matchAll(NUM)) {
    const t = norm(m[0]);
    numbers.push(m[0]);
    if (allow.has(t)) counts.json++;
    else if (der.has(t)) counts.derived++;
    else if (STRUCTURE.has(t)) counts.structure++;
    else {
      const ctx = all.slice(Math.max(0, m.index - 50), m.index + 20).replace(/\s+/g, " ");
      unmatched.push(`${m[0]}  ...${ctx}...`);
    }
  }
  return { text, svg, numbers, unmatched, counts };
}"""

HEAD = r"""() => {
  const q = (s, a = "content") => document.head.querySelector(s)?.getAttribute(a) ?? null;
  const abstract = document.querySelector('[data-v="abstract"]').textContent.trim();
  return {
    title: document.title, h1: document.querySelector("h1").textContent.trim(),
    abstract, description: q('meta[name="description"]'),
    canonical: q('link[rel="canonical"]', "href"),
    og: Object.fromEntries([...document.head.querySelectorAll('meta[property^="og:"]')]
      .map((m) => [m.getAttribute("property"), m.content])),
    twitter: Object.fromEntries([...document.head.querySelectorAll('meta[name^="twitter:"]')]
      .map((m) => [m.getAttribute("name"), m.content])),
    conj: document.querySelector("#s4-2").textContent.trim(),
    required: ["post hoc", "by construction", "prior-dependent"]
      .filter((w) => !document.body.textContent.includes(w)),
    derived: JSON.parse(document.getElementById("derived-numbers")?.textContent ?? "null"),
  };
}"""


def png_size(path: Path) -> tuple[int, int]:
    head = path.read_bytes()[:24]
    if head[:8] != b"\x89PNG\r\n\x1a\n":
        raise ValueError(f"{path} is not a PNG")
    return struct.unpack(">II", head[16:24])


def check(out: Path, base_url: str, path: str = "") -> list[str]:
    """One page of the site: ``path`` under the output directory ("" is the main paper)."""
    fails: list[str] = []
    image = base_url + "plate-i.png"
    base_url += path
    structure = STRUCTURE[path]
    problems: list[str] = []
    fetched: list[str] = []
    with serve(out) as root, sync_playwright() as p:
        url = root + path
        browser = launch(p)
        # 1. JavaScript off: the built HTML as served.
        ctx = browser.new_context(
            java_script_enabled=False, viewport=VIEWPORT, reduced_motion="reduce"
        )
        static = ctx.new_page()
        static.goto(url)
        static.wait_for_load_state("load")
        head = static.evaluate(HEAD)
        # 2. JavaScript on: redrawn from site/data/ over the prerendered page.
        live = open_rendered(browser, url, VIEWPORT, problems, fetched=fetched)
        data = [json.loads((out / "data" / name).read_text()) for name in sorted(set(fetched))]
        live_derived = live.evaluate("globalThis.__derived.map((d) => d.x)")
        if head["derived"] is None:
            fails.append("static page: no embedded #derived-numbers registry")
            static_derived: list = []
        else:
            static_derived = [d["x"] for d in head["derived"]]
        a_static = static.evaluate(
            AUDIT, {"json": data, "derived": static_derived, "structure": structure}
        )
        a_live = live.evaluate(
            AUDIT, {"json": data, "derived": live_derived, "structure": structure}
        )
        browser.close()

    fails += [f"JS page: {p}" for p in problems]
    for name, a in (("static (JS off)", a_static), ("JS", a_live)):
        c = a["counts"]
        print(f"{name}: {len(a['numbers'])} numbers (json {c['json']}, derived {c['derived']}, "
              f"structure {c['structure']}), {len(a['unmatched'])} unmatched")  # fmt: skip
        fails += [f"{name} unmatched: {u}" for u in a["unmatched"]]
    if a_static["numbers"] != a_live["numbers"]:
        sn, jn = a_static["numbers"], a_live["numbers"]
        pairs = enumerate(zip(sn, jn, strict=False))
        i = next((k for k, (x, y) in pairs if x != y), min(len(sn), len(jn)))
        fails.append(f"numbers differ at #{i}: static {a_static['numbers'][i:i + 5]} vs JS "
                     f"{a_live['numbers'][i:i + 5]} ({len(a_static['numbers'])} vs "
                     f"{len(a_live['numbers'])})")  # fmt: skip
    if a_static["text"] != a_live["text"]:
        s, j = a_static["text"].splitlines(), a_live["text"].splitlines()
        diff = [f"static: {x!r} | JS: {y!r}" for x, y in zip(s, j, strict=False) if x != y][:3]
        fails.append(f"text differs ({len(s)} vs {len(j)} lines): {diff}")
    if a_static["svg"] != a_live["svg"]:
        fails.append(f"SVG labels differ: {len(a_static['svg'])} vs {len(a_live['svg'])}")
    same = {
        k: "identical" if a_static[k] == a_live[k] else "DIFFER" for k in ("numbers", "text", "svg")
    }
    print(f"static vs JS: numbers {same['numbers']}, text {same['text']}, "
          f"SVG labels {same['svg']} "
          f"({len(a_static['svg'])})")  # fmt: skip

    # 3. Head tags.
    first = head["description"] or ""
    want = {
        "title = h1": head["title"] == head["h1"],
        "description is the abstract's first sentence": bool(first)
        and " ".join(head["abstract"].split()).startswith(first)
        and first.endswith((".", "!", "?")),
        "canonical": head["canonical"] == base_url,
        "og:title/description/url/image": head["og"].get("og:title") == head["h1"]
        and head["og"].get("og:description") == first
        and head["og"].get("og:url") == base_url
        and head["og"].get("og:image") == image,
        "twitter card": head["twitter"].get("twitter:card") == "summary_large_image"
        and head["twitter"].get("twitter:image") == image,
        "plate-i.png 1200x630": png_size(out / "plate-i.png")
        == (OG_SIZE["width"], OG_SIZE["height"]),
    }
    if path == "":
        want["section 4.2 heading complete"] = head["conj"].rstrip()[-4:].isdigit()
    else:  # the companion page states its deviations where it makes the claims
        want["the text says post hoc, by construction and prior-dependent"] = not head["required"]
    for k, ok in want.items():
        print(f"{'ok  ' if ok else 'FAIL'} {k}")
        if not ok:
            fails.append(f"head: {k}")
    print(f"description: {first}")
    return fails


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("out", type=Path)
    ap.add_argument("--base-url", default=BASE_URL)
    a = ap.parse_args()
    base = a.base_url if a.base_url.endswith("/") else a.base_url + "/"
    fails: list[str] = []
    for path in PAGES:
        print(f"--- {path or 'main page'}")
        fails += [f"[{path or 'main'}] {f}" for f in check(a.out.resolve(), base, path)]
    if fails:
        print("\nFAILED:\n  " + "\n  ".join(fails))
        raise SystemExit(1)
    print("\nall checks pass")


if __name__ == "__main__":
    main()
