"""Page images for the transcription of Newcomb 1882, pp. 451–452 (docs/phase2/transcription.md).

    uv run --with pymupdf --with pillow python -I docs/phase2/scripts/mercury_crops.py

Reads the cached scan (data/cache/mercury/astronomicalpap00newcgoog.pdf, SHA-256 in
data/historical/mercury/MANIFEST.md), takes the embedded page scans (4854 × 6636 px) of PDF
pages 483 and 484 (printed 451 and 452), and cuts horizontal bands that show only the columns
transcribed: contact, date, observed time with its k term, weight. The tabular times and the
symbolic corrections to the right are outside every band. Output: data/cache/mercury/renders/
(ignored).
"""

import io
from pathlib import Path

import pymupdf
from PIL import Image

ROOT = Path(__file__).resolve().parents[3]
PDF = ROOT / "data" / "cache" / "mercury" / "astronomicalpap00newcgoog.pdf"
OUT = ROOT / "data" / "cache" / "mercury" / "renders"

# band name -> (left, top, right, bottom) in pixels of the embedded scan
FX, FY = 4854 / 1500, 6636 / 1848
P451 = {"head": (655, 750), "nov_1": (728, 870), "nov_2": (855, 990), "nov_3": (975, 1110),
        "nov_4": (1095, 1232), "nov_5": (1195, 1262), "may_1": (1285, 1400),
        "may_2": (1385, 1535), "may_3": (1440, 1545), "may_4": (1515, 1600)}  # fmt: skip
P452 = {"head": (1080, 1400), "nov_1": (1390, 2010), "nov_2": (1950, 2550), "may": (2790, 3300)}
BANDS = {
    451: {
        k: (int(285 * FX), int(a * FY), int(700 * FX), int(b * FY)) for k, (a, b) in P451.items()
    },
    452: {k: (1390, a, 2790, b) for k, (a, b) in P452.items()},
}
PDF_PAGE = {451: 483, 452: 484}  # 1-based


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    doc = pymupdf.open(PDF)
    for printed, bands in BANDS.items():
        page = doc[PDF_PAGE[printed] - 1]
        raw = doc.extract_image(page.get_images()[0][0])["image"]
        im = Image.open(io.BytesIO(raw)).convert("L")
        for name, box in bands.items():
            im.crop(box).save(OUT / f"p{printed}_{name}.png")
            print(f"p{printed}_{name}.png", box)


if __name__ == "__main__":
    main()
