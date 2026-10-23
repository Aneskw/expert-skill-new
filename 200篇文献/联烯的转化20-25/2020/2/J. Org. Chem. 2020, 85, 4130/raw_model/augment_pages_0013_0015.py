#!/usr/bin/env python3
"""Add the visually reviewed bibliography pages 13--15 to the v3.1 bundle.

These pages contain only the article's reference list.  They are retained as
source context and page coverage, not converted into synthetic experiments.
"""
from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "runs/current"
DOC = "doc_5909d9a1e9253598"
CONTEXT_ID = f"{DOC}:context:bibliography_pages13_15"


def ev(page: int) -> list[str]:
    return [f"{DOC}:page{page}:image"]


bundle_path = OUT / "submitted.json"
bundle = json.loads(bundle_path.read_text())
# Formal submissions must reuse the prepared dataset's evidence registry
# exactly; pages 13--15 were staged after the preceding continuation merge.
bundle["documents"] = json.loads((OUT / "dataset.json").read_text())["documents"]

context_exists = any(
    c.get("context_id") == CONTEXT_ID for c in bundle.get("source_context", [])
)
if not context_exists:
    bundle["source_context"].append({
        "context_id": CONTEXT_ID,
        "document_id": DOC,
        "pages": [13, 14, 15],
        "kind": "bibliographic_references",
        "evidence_ids": ev(13) + ev(14) + ev(15),
        "note": (
            "The original rendered pages contain only the article bibliography "
            "(references 1--25, with entries continued across pages). Citation "
            "titles that mention reactions are prior-work metadata, not present-"
            "paper experimental procedures, reaction schemes, or independent "
            "literature-example arrows. No synthesis record was created."
        ),
    })

raw = {
    "batch_id": "allene-stereo-v31-200-20261017",
    "document_id": DOC,
    "visual_pages": [13, 14, 15],
    "page_findings": [
        {
            "page": 13,
            "content": "Two-column bibliography: continuation of reference 1, references 2--11, and the beginning of reference 12.",
            "reaction_figures": 0,
            "tables": 0,
            "reaction_arrows": 0,
            "footnotes": 0,
            "present_paper_experiments": 0,
        },
        {
            "page": 14,
            "content": "Two-column bibliography: continuation of reference 12, references 13--20, and the beginning of reference 21.",
            "reaction_figures": 0,
            "tables": 0,
            "reaction_arrows": 0,
            "footnotes": 0,
            "present_paper_experiments": 0,
        },
        {
            "page": 15,
            "content": "Final bibliography page: continuation of reference 21 and references 22--25; remainder of page is blank journal layout.",
            "reaction_figures": 0,
            "tables": 0,
            "reaction_arrows": 0,
            "footnotes": 0,
            "present_paper_experiments": 0,
        },
    ],
    "decision": "Retain as bibliographic source_context; create no molecule or reaction records.",
    "limitations": [
        "Citation titles are not evidence for detailed conditions or outcomes in the cited papers.",
        "No supporting-information association is inferred from bibliography entries.",
        "Visual review is restricted to PDF pages 13--15 in this continuation.",
    ],
}
(OUT / "raw_model/visual_transcription_pages_0013_0015.json").write_text(
    json.dumps(raw, ensure_ascii=False, indent=2) + "\n"
)

coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {
            "page": 13,
            "state": "attempted",
            "regions": [
                {
                    "region": "full page, two-column REFERENCES section",
                    "finding": "Bibliography only: continuation of reference 1, references 2--11, and start of reference 12; no present-paper experiment, reaction scheme, table, arrow, or footnote.",
                    "evidence_ids": ev(13),
                }
            ],
            "record_ids": [],
            "unresolved": [
                "Detailed experimental facts in the cited papers cannot be derived from citation titles and were not imported."
            ],
        },
        {
            "page": 14,
            "state": "attempted",
            "regions": [
                {
                    "region": "full page, two-column REFERENCES continuation",
                    "finding": "Bibliography only: continuation of reference 12, references 13--20, and start of reference 21; cited reaction titles are prior-work metadata, not experiments reported by this article.",
                    "evidence_ids": ev(14),
                }
            ],
            "record_ids": [],
            "unresolved": [
                "No cited-paper full text or explicitly associated source was supplied, so citation entries are not expanded into literature experiments."
            ],
        },
        {
            "page": 15,
            "state": "attempted",
            "regions": [
                {
                    "region": "upper-left REFERENCES continuation and remaining blank journal page",
                    "finding": "Bibliography ends with continuation of reference 21 and references 22--25; no synthesis experiment, mechanism figure, table, arrow, footnote, or analytical peak table is present.",
                    "evidence_ids": ev(15),
                }
            ],
            "record_ids": [],
            "unresolved": [],
        },
    ],
}
(OUT / "coverage/pages_0013_0015.json").write_text(
    json.dumps(coverage, ensure_ascii=False, indent=2) + "\n"
)

bundle_path.write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")

print(json.dumps({
    "pages": [13, 14, 15],
    "new_reactions": 0,
    "new_molecules": 0,
    "new_source_context": 0 if context_exists else 1,
    "total_reactions": len(bundle["reactions"]),
    "total_molecules": len(bundle["molecules"]),
}, ensure_ascii=False))
