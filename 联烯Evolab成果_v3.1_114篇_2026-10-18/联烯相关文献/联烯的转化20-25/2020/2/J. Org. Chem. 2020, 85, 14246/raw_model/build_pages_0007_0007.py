#!/usr/bin/env python3
"""Add the source-bounded page-7 review without creating synthetic records."""
from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
RUN = ROOT / "runs/current"
DOC = "doc_7ca2cfc676d4ba33"
CONTEXT_ID = f"{DOC}:context_page7_notes_acknowledgments_references"


def write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n")


def main() -> None:
    dataset = json.loads((RUN / "dataset.json").read_text())
    document = next(item for item in dataset["documents"] if item["document_id"] == DOC)
    evidence_ids = {item["evidence_id"] for item in document["evidence"]}
    required = {f"{DOC}:page7:text", f"{DOC}:page7:image"}
    missing = sorted(required - evidence_ids)
    if missing:
        raise SystemExit(f"page-7 evidence missing from dataset: {missing}")

    bundle_path = RUN / "submitted.json"
    bundle = json.loads(bundle_path.read_text())
    # The prepared dataset gained the rendered page-7 image after the earlier
    # page-1–6 increment. Merge requires the document/evidence block to be an
    # exact mechanical copy; this does not alter any extracted chemistry.
    bundle["documents"] = dataset["documents"]
    contexts = bundle.setdefault("source_context", [])
    contexts[:] = [item for item in contexts if item.get("context_id") != CONTEXT_ID]
    contexts.append(
        {
            "context_id": CONTEXT_ID,
            "document_id": DOC,
            "kind": "bibliography_and_back_matter",
            "evidence_ids": [f"{DOC}:page7:text", f"{DOC}:page7:image"],
            "note": (
                "Page 7 contains the no-competing-financial-interest note, acknowledgments, "
                "and references 1–22. The citations describe prior literature only; the page "
                "contains no reaction arrow, experimental table, procedure, structure drawing, "
                "or standalone characterization dataset, so no molecule or reaction record was created."
            ),
        }
    )
    write_json(bundle_path, bundle)

    coverage = {
        "document_id": DOC,
        "semantic_accuracy_verified": False,
        "pages": [
            {
                "page": 7,
                "state": "attempted",
                "regions": [
                    {
                        "region": "upper-left Notes and Acknowledgments",
                        "finding": (
                            "Back matter only: the authors declare no competing financial interest "
                            "and acknowledge NRF support; no synthetic experiment or chemical result."
                        ),
                        "evidence_ids": [f"{DOC}:page7:image", f"{DOC}:page7:text"],
                    },
                    {
                        "region": "two-column References list (entries 1–22)",
                        "finding": (
                            "Bibliographic citations to prior work, including earlier allenoate and "
                            "natural-product studies. No reaction arrows, structures, tables, footnotes "
                            "with new conditions, or experiments performed in this paper occur on this page."
                        ),
                        "evidence_ids": [f"{DOC}:page7:image", f"{DOC}:page7:text"],
                    },
                ],
                "record_ids": [],
                "unresolved": [
                    "The cited papers themselves were not supplied in this bounded source and were not treated as experiments from the present article."
                ],
            }
        ],
    }
    write_json(RUN / "coverage/pages_0007_0007.json", coverage)

    raw = {
        "batch_id": "allene-stereo-v31-200-20261017",
        "document_id": DOC,
        "page_range": [7, 7],
        "semantic_accuracy_verified": False,
        "chemical_accuracy_verified": False,
        "visual_review": {"7": f"{DOC}:page7:image"},
        "page_inventory": [
            {
                "page": 7,
                "observed_regions": [
                    "Notes: no competing financial interest",
                    "Acknowledgments: NRF support",
                    "References 1–22 in two columns",
                ],
                "reaction_arrows": 0,
                "reaction_tables": 0,
                "experimental_procedures": 0,
                "structure_drawings": 0,
                "standalone_characterization_tables": 0,
            }
        ],
        "decisions": [
            "No molecule or reaction was created because page 7 is entirely back matter and bibliography.",
            "References are prior literature citations, not experiments performed in the present paper.",
            "Existing page-1–6 records were preserved only to avoid destructive replacement and are not counted as page-7 results.",
        ],
        "new_molecule_ids": [],
        "new_reaction_ids": [],
        "counts": {"new_molecules": 0, "new_reactions": 0, "new_source_context_items": 1},
    }
    write_json(RUN / "raw_model/pages_0007_0007_raw_extraction.json", raw)

    workflow = ROOT / "runs/workflow.md"
    marker = "## Increment — original PDF page 7"
    text = workflow.read_text()
    if marker not in text:
        text += f"""

{marker}

- Bounded scope for this increment is exactly original PDF page **7–7**. The rendered source image `runs/current/assets/{DOC}/page7.jpg` was opened at original detail; pages 1–6 were not re-claimed as covered here.
- Page 7 contains only the no-competing-financial-interest note, acknowledgments, and references 1–22. It has no reaction arrow, experimental table, procedure, structure drawing, standalone spectrum/peak table, or new chemical footnote.
- Reference entries are classified as prior literature/back matter and were not converted into experiments from this paper. Current increment: **0 molecules, 0 reactions, 1 source-context item**.
- No stereochemical or bond-formation review was applicable to a new page-7 molecule or reaction because none is depicted. Existing machine candidates from pages 1–6 remain unchanged and `needs_review`.
- SI association remains `unestablished`; no SI was inferred from the bibliography or back matter.
- The untouched `standard init` template for this increment is `runs/current/raw_model/pages_0007_0007_template.json`; the visual decision inventory is `runs/current/raw_model/pages_0007_0007_raw_extraction.json`; exact coverage is `runs/current/coverage/pages_0007_0007.json`.
"""
        workflow.write_text(text)

    print(
        json.dumps(
            {
                "page_range": [7, 7],
                "new_molecules": 0,
                "new_reactions": 0,
                "new_source_context_items": 1,
                "bundle_molecules_preserved": len(bundle.get("molecules", [])),
                "bundle_reactions_preserved": len(bundle.get("reactions", [])),
            },
            ensure_ascii=False,
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
