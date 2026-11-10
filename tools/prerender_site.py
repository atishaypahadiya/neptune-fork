"""Build the Pages output: site/ prerendered, so the page reads without JavaScript.

The site renders its text, numbers, tables, references and figures client-side from
site/data/*.json. This script copies site/ to the output directory, loads it in headless Chrome
(Playwright) until the page reports ``data-ready``, and writes the rendered DOM back as
index.html: text, tables, captions, references and the static SVG of every plate and figure.
JavaScript still runs on the built page and redraws the same content, adding what needs JS
(hover links, Plate III's animation and controls, Copy buttons; those controls are
``data-js hidden`` in the HTML).

It also writes the head tags a link preview reads (title, description = the abstract's first
sentence, canonical URL, Open Graph and Twitter card) and plate-i.png, Plate I at 1200 x 630,
as their image. The number registry the audit needs (``drv`` in site/js/data.js) is embedded
as ``<script type="application/json" id="derived-numbers">``. Check the output with
``tools/check_prerender.py``.

    uv run --group site python tools/prerender_site.py --out _site [--base-url URL]
"""

from __future__ import annotations

import argparse
import contextlib
import functools
import http.server
import shutil
import threading
from collections.abc import Iterator
from pathlib import Path

from playwright.sync_api import Browser, BrowserContext, Page, sync_playwright

ROOT = Path(__file__).resolve().parents[1]
SITE = ROOT / "site"
BASE_URL = "https://atishaypahadiya.github.io/neptune-fork/"
VIEWPORT = {"width": 1440, "height": 900}  # the drawn width of the static SVG
OG_SIZE = {"width": 1200, "height": 630}
OPTIONAL = ("m8_arbiter.json", "m8b.json")  # loaded with { optional: true }: a 404 is "pending"


class _Server(http.server.ThreadingHTTPServer):
    request_queue_size = 128  # the default 5 drops connections when the page fetches 16 files


@contextlib.contextmanager
def serve(directory: Path) -> Iterator[str]:
    """Serve ``directory`` on a free localhost port; yields the base URL."""

    class Quiet(http.server.SimpleHTTPRequestHandler):
        def log_message(self, *_args: object) -> None:
            pass

    srv = _Server(("127.0.0.1", 0), functools.partial(Quiet, directory=str(directory)))
    t = threading.Thread(target=srv.serve_forever, daemon=True)
    t.start()
    try:
        yield f"http://127.0.0.1:{srv.server_address[1]}/"
    finally:
        srv.shutdown()
        srv.server_close()


def launch(p) -> Browser:
    """Installed Google Chrome if there is one (GitHub's Ubuntu runners have it), else
    Playwright's own Chromium (``playwright install chromium``)."""
    try:
        return p.chromium.launch(channel="chrome")
    except Exception:
        return p.chromium.launch()


def open_rendered(
    browser: Browser | BrowserContext,
    url: str,
    viewport: dict,
    problems: list[str],
    fetched: list[str] | None = None,
) -> Page:
    """Load the page with JS, reduced motion (Plate III holds its conjunction frame), and wait
    for data-ready. Page errors, failed requests and HTTP errors (bar optional files) go to
    ``problems``; the names of the data files it loaded go to ``fetched``. Given a browser
    context (a phone, say), the page opens in it with the context's own settings."""
    if isinstance(browser, BrowserContext):
        page = browser.new_page()
    else:
        page = browser.new_page(viewport=viewport, reduced_motion="reduce", device_scale_factor=1)
    page.on("pageerror", lambda e: problems.append(f"pageerror: {e}"))
    page.on("requestfailed", lambda r: problems.append(f"request failed: {r.url}"))
    page.on(
        "response",
        lambda r: (
            problems.append(f"HTTP {r.status}: {r.url}")
            if r.status >= 400 and not r.url.endswith(OPTIONAL)
            else None
        ),
    )
    if fetched is not None:
        page.on(
            "response",
            lambda r: (
                fetched.append(r.url.rsplit("/", 1)[1])
                if r.ok and "/data/" in r.url and r.url.endswith(".json")
                else None
            ),
        )
    page.goto(url)
    page.wait_for_function(
        "document.documentElement.dataset.ready === '1' || "
        "!!document.querySelector('.error:not([hidden])')",
        timeout=90_000,
    )
    page.evaluate("document.fonts.ready")
    errors = page.eval_on_selector_all(".error:not([hidden])", "es => es.map(e => e.textContent)")
    problems.extend(f"figure error: {e}" for e in errors)
    return page


# Runs in the rendered page: head tags, JS-only controls hidden, the number registry embedded.
SERIALIZE = """({base, og}) => {
  const html = document.documentElement, head = document.head;
  delete html.dataset.ready;
  html.dataset.prerendered = "1";
  for (const el of document.querySelectorAll("[data-js]")) el.hidden = true;

  const title = document.querySelector("h1").textContent.trim();
  const abstract = document.querySelector('[data-v="abstract"]').textContent.trim();
  // First sentence: up to the first full stop followed by a space and a capital letter.
  const m = abstract.match(/^.+?[.!?](?=\\s+[A-Z(])/s);
  const description = (m ? m[0] : abstract).replace(/\\s+/g, " ");
  const caption = document.querySelector("#cap-plate-i").textContent.trim().replace(/\\s+/g, " ");
  document.title = title;

  const set = (sel, make, attrs) => {
    let el = head.querySelector(sel);
    if (!el) { el = document.createElement(make); head.append(el); }
    for (const [k, v] of Object.entries(attrs)) el.setAttribute(k, v);
  };
  const meta = (key, value, attr = "property") =>
    set(`meta[${attr}="${key}"]`, "meta", { [attr]: key, content: value });
  meta("description", description, "name");
  set('link[rel="canonical"]', "link", { rel: "canonical", href: base });
  const image = new URL("plate-i.png", base).href;
  meta("og:type", "article");
  meta("og:site_name", "neptune-fork");
  meta("og:title", title);
  meta("og:description", description);
  meta("og:url", base);
  meta("og:image", image);
  meta("og:image:width", String(og.width));
  meta("og:image:height", String(og.height));
  meta("og:image:alt", caption);
  meta("twitter:card", "summary_large_image", "name");
  meta("twitter:title", title, "name");
  meta("twitter:description", description, "name");
  meta("twitter:image", image, "name");
  meta("twitter:image:alt", caption, "name");

  const reg = document.createElement("script");
  reg.type = "application/json";
  reg.id = "derived-numbers";
  reg.textContent = JSON.stringify(globalThis.__derived);
  document.body.append(reg);
  return "<!doctype html>\\n" + html.outerHTML + "\\n";
}"""


def build(out: Path, base_url: str) -> None:
    if out.exists():
        shutil.rmtree(out)
    shutil.copytree(SITE, out)
    problems: list[str] = []
    with serve(out) as url, sync_playwright() as p:
        browser = launch(p)
        page = open_rendered(browser, url, VIEWPORT, problems)
        html = page.evaluate(SERIALIZE, {"base": base_url, "og": OG_SIZE})
        og = open_rendered(browser, url, OG_SIZE, problems)
        og.screenshot(path=str(out / "plate-i.png"), clip={"x": 0, "y": 0, **OG_SIZE})
        browser.close()
    if problems:
        raise SystemExit("prerender failed:\n  " + "\n  ".join(problems))
    (out / "index.html").write_text(html, encoding="utf-8")
    print(f"wrote {out / 'index.html'} ({len(html) / 1024:.0f} KB) and {out / 'plate-i.png'}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", type=Path, default=ROOT / "_site")
    ap.add_argument("--base-url", default=BASE_URL, help="the deployed URL, ending in /")
    a = ap.parse_args()
    build(a.out.resolve(), a.base_url if a.base_url.endswith("/") else a.base_url + "/")


if __name__ == "__main__":
    main()
