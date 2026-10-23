#!/usr/bin/env python3
"""Check representative generated Markdown against source-derived record tokens."""

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RECORDS = ROOT / "records"
OUTPUT = Path(__file__).with_name("verification_summary_markdown.json")

CHECKS = {
    "synthesis_complex4": [
        '"number": 400, "unit": "mg"',
        '"number": 0.57, "unit": "mmol"',
        '"number": 58.6, "unit": "mg"',
        '"number": 1.15, "unit": "mmol"',
        '"number": 15, "unit": "mL"',
        '"number": 326, "unit": "mg"',
        '"number": 82, "unit": "%"',
        "overnight",
    ],
    "synthesis_complex5": [
        '"number": 200, "unit": "mg"',
        '"number": 0.26, "unit": "mmol"',
        '"number": 20, "unit": "mg"',
        '"number": 0.53, "unit": "mmol"',
        '"number": 15, "unit": "mL"',
        '"number": 40, "unit": "°C"',
        "one week",
        '"number": 152, "unit": "mg"',
        '"number": 76, "unit": "%"',
    ],
    "isolation_complex6": [
        '"number": 300, "unit": "mg"',
        '"number": 0.43, "unit": "mmol"',
        "0.27 μL as printed",
        "0.43 mmol at 1.6 M implies about 0.27 mL",
        '"number": 7, "unit": "mL"',
        '"number": 40, "unit": "°C"',
        "overnight",
        '"number": 122, "unit": "mg"',
        '"number": 40, "unit": "%"',
        '"raw_value": "89:11"',
        '"raw_value": "48:34:18"',
        "no chosen correction",
    ],
    "catalyst_loading_2_50": [
        '"number": 1.33, "unit": "mmol"',
        '"number": 4, "unit": "mL"',
        '"number": 0.33, "unit": "M"',
        '"number": 2.5, "unit": "mol%"',
        '"number": 0.85, "unit": "equiv"',
    ],
}


def main() -> None:
    results = []
    for suffix, tokens in CHECKS.items():
        path = RECORDS / f"doc_837adcd3d5d8cb7b__{suffix}.md"
        body = path.read_text(encoding="utf-8")
        missing = [token for token in tokens if token not in body]
        results.append(
            {
                "record_id": f"doc_837adcd3d5d8cb7b:{suffix}",
                "markdown_path": str(path.relative_to(ROOT)),
                "checked_tokens": tokens,
                "missing_tokens": missing,
            }
        )
    payload = {
        "purpose": "Representative post-merge Markdown token checks; not a chemistry-accuracy claim.",
        "checks": results,
        "all_missing_token_lists_empty": all(not row["missing_tokens"] for row in results),
    }
    OUTPUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
