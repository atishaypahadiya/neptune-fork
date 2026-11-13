"""Public-facing text: no internal labels, no "we", no "the project's".

Read here: README.md, CHANGELOG.md, CITATION.cff, .zenodo.json and the string literals of
the site's copy modules (the pages are filled from them; ``tools/check_prerender.py`` checks
the built pages). Paths, link targets and code spans are not prose and are removed first.
The identifier table that maps each label to its public name is in ``docs/reproduce.md``.
"""

import json
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]

# Internal labels. S4, S5 and S6 are the companion study's own defined labels and are allowed.
LABELS = re.compile(
    r"\bphase[ -]?(?:1b|1|2|one|two)\b|\bP1B\b|\bM\d{1,2}[a-e]?\b|\bS1 cell 0\b",
    re.IGNORECASE,
)
PLURAL = re.compile(r"\b(?:we|our|ours|us|ourselves)\b", re.IGNORECASE)
POSSESSIVE = re.compile(r"\bthe project[’']s\b", re.IGNORECASE)

# Quotations that may carry any of the words, verbatim. Add the quoted sentence, not a word.
QUOTATION_ALLOWLIST: tuple[str, ...] = ()
# Look-alikes that are not the words, with the reason.
NOT_PLURAL = ("US",)  # the country


def _markdown_prose(text: str) -> str:
    text = re.sub(r"```.*?```", " ", text, flags=re.DOTALL)  # code blocks
    text = re.sub(r"`[^`\n]*`", " ", text)  # code spans: paths, modules, commands
    text = re.sub(r"\]\([^)\s]*\)", "]", text)  # link targets
    text = re.sub(r"<[^>]+>", " ", text)  # raw HTML tags and autolinks
    return re.sub(r"https?://\S+", " ", text)


def _js_strings(text: str) -> str:
    text = re.sub(r"^\s*//.*$", " ", text, flags=re.MULTILINE)
    text = re.sub(r"/\*.*?\*/", " ", text, flags=re.DOTALL)
    out = []
    for m in re.finditer(r'`((?:\\.|[^`\\])*)`|"((?:\\.|[^"\\\n])*)"', text):
        s = m.group(1) if m.group(1) is not None else m.group(2)
        out.append(re.sub(r"\$\{[^}]*\}", " ", s))
    return "\n".join(out)


def _yaml_prose(text: str) -> str:
    return "\n".join(
        ln for ln in text.splitlines() if not re.match(r"\s*(url|repository|doi|value)", ln)
    )


def _zenodo_prose(text: str) -> str:
    meta = json.loads(text)
    return re.sub(
        r"<[^>]+>", " ", " ".join([meta["title"], meta["description"], meta.get("notes", "")])
    )


SURFACES = {
    "README.md": _markdown_prose,
    "CHANGELOG.md": _markdown_prose,
    "CITATION.cff": _yaml_prose,
    ".zenodo.json": _zenodo_prose,
    "site/js/copy.js": _js_strings,
    "site/js/copy_p1b.js": _js_strings,
}


def _prose(name: str) -> str:
    text = SURFACES[name]((ROOT / name).read_text(encoding="utf-8"))
    for quote in QUOTATION_ALLOWLIST:
        text = text.replace(quote, " ")
    return text


def _hits(pattern: re.Pattern, text: str, skip: tuple[str, ...] = ()) -> list[str]:
    found = []
    for m in pattern.finditer(text):
        if m.group(0) in skip:
            continue
        found.append(text[max(0, m.start() - 40) : m.end() + 40].replace("\n", " "))
    return found


@pytest.mark.parametrize("name", sorted(SURFACES))
def test_no_internal_labels(name):
    assert not _hits(LABELS, _prose(name))


@pytest.mark.parametrize("name", sorted(SURFACES))
def test_first_person_singular(name):
    assert not _hits(PLURAL, _prose(name), NOT_PLURAL)


@pytest.mark.parametrize("name", sorted(SURFACES))
def test_no_possessive_for_the_work(name):
    assert not _hits(POSSESSIVE, _prose(name))


def test_the_patterns_catch_what_they_should():
    bad = "In Phase 1b we ran M10 and the project's fit, see P1B and phase 2; our M6c."
    assert len(_hits(LABELS, bad)) == 5
    assert len(_hits(PLURAL, bad)) == 2
    assert len(_hits(POSSESSIVE, bad)) == 1
    ok = _markdown_prose(
        "See [the write-up](docs/notes/P1B_writeup.md) and `fork.export.m10`; S4 and S5."
    )
    assert not _hits(LABELS, ok)
    assert not _hits(PLURAL, "The US edition; a bonus; the focus.", NOT_PLURAL)
