#!/usr/bin/env python3
"""Verify representative page-7--12 Markdown quantities and conflicts."""

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RECORDS = ROOT / "records"
OUTPUT = Path(__file__).with_name("verification_pages_7_12_markdown.json")
DOC = "doc_837adcd3d5d8cb7b"

CHECKS = {
    "synthesis_complex7_method_a": [
        '"number": 300, "unit": "mg"', '"number": 0.43, "unit": "mmol"',
        '"number": 20, "unit": "mL"', '"number": -78, "unit": "°C"',
        '"number": 1.08, "unit": "mL"', '"number": 1.08, "unit": "mmol"',
        '"number": 10, "unit": "min"', "overnight",
        '"number": 78, "unit": "mg"', '"number": 26, "unit": "%"',
    ],
    "synthesis_complex8": [
        '"number": 64, "unit": "mg"', '"number": 0.18, "unit": "mmol"',
        '"number": 19, "unit": "mg"', '"number": 10, "unit": "mL"',
        '"number": -78, "unit": "°C"', "overnight",
        '"number": 53, "unit": "mg"', '"number": 83, "unit": "%"',
    ],
    "synthesis_complex7_method_b": [
        '"number": 16, "unit": "mg"', '"number": 0.04, "unit": "mmol"',
        '"number": 10.8, "unit": "mg"', '"number": 0.6, "unit": "mL"',
        '"number": 75, "unit": "°C"', '"number": 3, "unit": "d"',
        '"number": 12, "unit": "mg"', '"number": 75, "unit": "%"',
    ],
    "synthesis_complex10": [
        '"number": 240, "unit": "mg"', '"number": 0.74, "unit": "mmol"',
        '"number": 350, "unit": "mg"', '"number": 50, "unit": "mL"',
        '"number": 60, "unit": "°C"', "overnight", '"number": 0.1, "unit": "mL"',
        '"number": 130, "unit": "mg"', '"number": 33, "unit": "%"',
    ],
    "synthesis_complex9_from10": [
        '"number": 133, "unit": "mg"', '"number": 0.24, "unit": "mmol"',
        '"number": 600, "unit": "mg"', '"number": 10.18, "unit": "mmol"',
        '"number": 15, "unit": "mL"', "overnight", "4 x 3 mL",
        '"number": 97, "unit": "mg"', '"number": 72, "unit": "%"',
    ],
    "synthesis_complex9_from7": [
        "25 mg, 0.83 mmol", "0.040 mmol", "no correction selected",
        '"number": 45, "unit": "mg"', '"number": 0.75, "unit": "mmol"',
        '"number": 0.6, "unit": "mL"', '"number": 1, "unit": "d"',
        "closed tube; H2 pressure allowed to build",
    ],
    "catalysis_complex7": [
        "[H3B·NMe2H] = 1.33 M", "1.33 mmol in 4 mL", "approximately 0.33 M",
        '"data_status": "conflict"', "concentration conflict is unresolved",
    ],
}


def main() -> None:
    rows = []
    for suffix, tokens in CHECKS.items():
        path = RECORDS / f"{DOC}__{suffix}.md"
        body = path.read_text(encoding="utf-8")
        missing = [token for token in tokens if token not in body]
        rows.append({
            "record_id": f"{DOC}:{suffix}",
            "markdown_path": str(path.relative_to(ROOT)),
            "checked_tokens": tokens,
            "missing_tokens": missing,
        })
    payload = {
        "purpose": "Post-merge page-7--12 Markdown quantity/conflict checks; not a chemistry-accuracy claim.",
        "checks": rows,
        "all_missing_token_lists_empty": all(not row["missing_tokens"] for row in rows),
    }
    OUTPUT.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, ensure_ascii=False, indent=2))
    raise SystemExit(0 if payload["all_missing_token_lists_empty"] else 1)


if __name__ == "__main__":
    main()
