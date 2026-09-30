"""Splice the exported bundle into the artifact template to produce one self-contained page.

The page has no backend: the coverage model is ported to JavaScript and the data tables are
inlined, so the hosted atlas recomputes every number in the browser rather than trusting a
precomputed answer. Inlining rather than fetching also means the first frame is complete.
"""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TEMPLATE = ROOT / "web" / "artifact" / "template.html"
BUNDLE = ROOT / "web" / "public" / "venomgap-bundle.json"
OUT = ROOT / "web" / "artifact" / "venomgap-atlas.html"

PLACEHOLDER = "__BUNDLE__"


def main() -> None:
    template = TEMPLATE.read_text()
    if PLACEHOLDER not in template:
        raise SystemExit(f"{TEMPLATE} has no {PLACEHOLDER} placeholder")
    bundle = BUNDLE.read_text()

    # The payload sits inside <script type="application/json">, so the only sequence that can
    # break out is a literal "</script". Nothing else needs escaping.
    bundle = bundle.replace("</", "<\\/")
    json.loads(bundle.replace("<\\/", "</"))  # the escape must round-trip

    page = template.replace(PLACEHOLDER, bundle)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(page)
    size_mb = OUT.stat().st_size / 1024 / 1024
    print(f"wrote {OUT} ({size_mb:.2f} MB)")
    if size_mb > 15.0:
        raise SystemExit("page exceeds the 16 MB artifact limit")


if __name__ == "__main__":
    main()
