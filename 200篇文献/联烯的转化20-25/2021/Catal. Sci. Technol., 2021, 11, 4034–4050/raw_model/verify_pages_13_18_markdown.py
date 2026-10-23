#!/usr/bin/env python3
"""Check generated expert Markdown for page-13-18 quantities/conflicts."""

from __future__ import annotations

import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
OUT = HERE.parent
PREFIX = "doc_837adcd3d5d8cb7b__"

expected = {
    "synthesis_complex11": [
        '"number": 2.88, "unit": "g"',
        '"number": 10.59, "unit": "mmol"',
        '"number": 1.12, "unit": "g"',
        '"number": 12.71, "unit": "mmol"',
        '"number": 1.2, "unit": "equiv"',
        '"number": 100, "unit": "mL"',
        '"number": -78, "unit": "°C"',
        "room temperature",
        "overnight",
        '"number": 1.08, "unit": "g"',
        '"number": 32, "unit": "%"',
        "exclusion of light",
    ],
    "synthesis_li_amidoborane": [
        '"number": 400, "unit": "mg"',
        '"number": 6.8, "unit": "mmol"',
        '"number": 15, "unit": "mL"',
        '"number": 2.9, "unit": "mL"',
        '"number": 7.5, "unit": "mmol"',
        '"number": 2.5, "unit": "M"',
        '"number": -78, "unit": "°C"',
        "room temperature",
        "overnight",
        "2 × 5 mL",
        '"number": 350, "unit": "mg"',
        '"number": 87, "unit": "%"',
        "identity of the evolved gas is not stated",
    ],
    "synthesis_mu_aminodiborane": [
        '"number": 30, "unit": "mg"',
        '"number": 0.51, "unit": "mmol"',
        '"number": 0.6, "unit": "mL"',
        '"number": 48, "unit": "uL"',
        '"number": 80, "unit": "°C"',
        "overnight",
        "-19.6 ppm",
        "96.32 MHz",
        "no isolated mass or yield is given",
    ],
    "isolation_complex6": [
        "0.27 μL",
        "0.43 mmol at 1.6 M implies about 0.27 mL",
        '"data_status": "conflict"',
        "Printed MeLi volume/amount/concentration conflict remains unresolved",
    ],
    "synthesis_complex9_from7": [
        "25 mg, 0.83 mmol",
        "25 mg at the reported M = 626 implies about 0.040 mmol",
        '"data_status": "conflict"',
        "Printed complex-7 amount conflict remains unresolved",
    ],
}

results = {}
for slug, tokens in expected.items():
    path = OUT / "records" / f"{PREFIX}{slug}.md"
    text = path.read_text()
    results[slug] = {
        "path": str(path.relative_to(OUT)),
        "checked_tokens": tokens,
        "missing_tokens": [token for token in tokens if token not in text],
    }

payload = {
    "pages": list(range(13, 19)),
    "all_checks_passed": all(not item["missing_tokens"] for item in results.values()),
    "records": results,
}
(HERE / "verification_pages_13_18_markdown.json").write_text(
    json.dumps(payload, ensure_ascii=False, indent=2) + "\n"
)
print(json.dumps(payload, ensure_ascii=False, indent=2))
