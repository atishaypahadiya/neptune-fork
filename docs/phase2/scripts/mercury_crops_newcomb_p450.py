"""Page bands for the transcription of Newcomb 1882, p. 450: the passage on the factor k.

    uv run --with pymupdf --with pillow python -I docs/phase2/scripts/mercury_crops_newcomb_p450.py

Reads the cached scan (data/cache/mercury/astronomicalpap00newcgoog.pdf, SHA-256 in
data/historical/mercury/MANIFEST.md) and takes the embedded page scan (4854 × 6636 px) of PDF
page 482, printed p. 450: the page before the comparison table of pp. 451–452. The page is
running text with one displayed formula and no table. The bands are whole lines of text, cut
in the white between lines, and together cover every printed line of the page. Output:
data/cache/mercury/renders_newcomb_p450/ (ignored).
"""

from pathlib import Path

import numpy as np
import pymupdf
from PIL import Image

ROOT = Path(__file__).resolve().parents[3]
PDF = ROOT / "data" / "cache" / "mercury" / "astronomicalpap00newcgoog.pdf"
OUT = ROOT / "data" / "cache" / "mercury" / "renders_newcomb_p450"
PDF_PAGE = 482  # 1-based; printed p. 450
LEFT, RIGHT = 330, 4240

# band -> (top, bottom, what it holds) in pixels of the embedded scan
BANDS = {
    "b00": (480, 650, "running head"),
    "b01": (690, 1290, "lines 1-5"),
    "b02": (1285, 1880, "lines 6-10"),
    "b03": (1875, 2590, "lines 11-16"),
    "b04": (2590, 3080, "lines 17-20"),
    "b05": (3080, 3860, "lines 21-25 and the displayed formula"),
    "b06": (3860, 4520, "lines 26-30"),
    "b07": (4520, 4950, "section number and title"),
    "b08": (4950, 5920, "lines 31-38"),
}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    doc = pymupdf.open(PDF)
    pix = pymupdf.Pixmap(doc, doc[PDF_PAGE - 1].get_images()[0][0])
    a = np.frombuffer(pix.samples, dtype=np.uint8).reshape(pix.height, pix.width, pix.n)[:, :, 0]
    im = Image.fromarray(a)
    for name, (top, bottom, what) in BANDS.items():
        im.crop((LEFT, top, RIGHT, bottom)).save(OUT / f"p450_{name}.png")
        print(f"p450_{name}.png", (LEFT, top, RIGHT, bottom), what)


if __name__ == "__main__":
    main()
