#!/usr/bin/env python3
"""Merge source-derived findings from original PDF pages 25--28.

No new chemical experiment begins on these pages.  Page 25 completes the
ALARM-MSPS analysis method and is linked to the five existing assay records;
pages 26--28 contain acknowledgments, abbreviations, and cited references.
"""

from __future__ import annotations

import copy
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
RUN = ROOT / "runs" / "current"
DOC = "doc_832eb8cea3728513"
PAGE25 = [f"{DOC}:page25:text", f"{DOC}:page25:image"]

REACTION_IDS = [
    f"{DOC}:rxn:rxn_p9_p10_alarm_msps_S17",
    f"{DOC}:rxn:rxn_p9_p10_alarm_msps_CPM",
    f"{DOC}:rxn:rxn_p9_p10_alarm_msps_51",
    f"{DOC}:rxn:rxn_p9_p10_alarm_msps_fluconazole",
    f"{DOC}:rxn:rxn_p9_p10_alarm_msps_DMSO",
]

METHOD_NOTE = (
    "Page 25 reports a minimum precursor signal of 5000 and HCD analysis at "
    "normalized collision energy 35 V, isolation width 2, and first mass m/z "
    "100.00. Each sample was run in technical triplicate and analyzed "
    "independently; peptide identification required all or all but one expected "
    "y ions and at least two b ions."
)


def evidence(page: int) -> list[str]:
    return [f"{DOC}:page{page}:text", f"{DOC}:page{page}:image"]


def context(context_id: str, page: int, region: str, classification: str, summary: str) -> dict:
    return {
        "context_id": context_id,
        "document_id": DOC,
        "location": {"page": page, "region": region},
        "classification": classification,
        "summary": summary,
        "evidence_ids": evidence(page),
        "review_status": "needs_review",
    }


bundle = json.loads((RUN / "submitted.json").read_text())
dataset = json.loads((RUN / "dataset.json").read_text())
bundle["documents"] = copy.deepcopy(dataset["documents"])

# Idempotently replace only this batch's contextual entries.
bundle["source_context"] = [
    item
    for item in bundle.get("source_context", [])
    if not item.get("context_id", "").startswith(("ctx_p25_", "ctx_p26_", "ctx_p27_", "ctx_p28_"))
]

reaction_by_id = {item["reaction_id"]: item for item in bundle["reactions"]}
missing = [rid for rid in REACTION_IDS if rid not in reaction_by_id]
if missing:
    raise SystemExit(f"Missing existing ALARM-MSPS record(s): {missing}")

for reaction_id in REACTION_IDS:
    reaction = reaction_by_id[reaction_id]
    reaction["evidence_ids"] = list(dict.fromkeys(reaction.get("evidence_ids", []) + PAGE25))
    reaction["review"]["status"] = "needs_review"

    reaction["notes"] = [
        note for note in reaction.get("notes", []) if note.get("raw_value") != METHOD_NOTE
    ]
    reaction["notes"].append(
        {
            "raw_value": METHOD_NOTE,
            "normalized_value": METHOD_NOTE,
            "acquisition_method": "original_extraction",
            "evidence_ids": PAGE25,
            "data_status": "pending_review",
            "verification": "unreviewed",
            "derivation": None,
            "alternatives": [],
        }
    )

    unresolved = reaction["review"].setdefault("unresolved", [])
    unresolved_note = (
        "Page 25 reports C232-containing peptide behavior as supporting context "
        "and selects the FSGDLDDQTC245R peptide for the study; no discrete C232 "
        "or C245 compound-adduct structure was isolated."
    )
    if unresolved_note not in unresolved:
        unresolved.append(unresolved_note)

    for outcome in reaction.get("outcomes", []):
        basis = outcome.get("measurement_basis")
        if isinstance(basis, dict):
            basis["evidence_ids"] = list(dict.fromkeys(basis.get("evidence_ids", []) + PAGE25))

bundle["source_context"].extend(
    [
        context(
            "ctx_p25_alarm_ms_analysis",
            25,
            "top-left continuation of La Antigen MS Analysis",
            "analytical_method",
            "ALARM-MSPS acquisition and manual-identification criteria: minimum precursor signal 5000; HCD normalized collision energy 35 V; isolation width 2; first mass m/z 100.00; technical triplicates analyzed independently. This completes the method for the five existing ALARM-MSPS test/control records and is not a new synthesis.",
        ),
        context(
            "ctx_p25_c232_interpretation",
            25,
            "left column C232 peptide interpretation",
            "mechanistic_interpretation",
            "C232-containing peptides showed slower/incomplete IAD labeling in DMSO and fluconazole controls. The authors therefore used FSGDLDDQTC245R for the principal analysis. This is assay interpretation/characterization, not an isolated chemical product.",
        ),
        context(
            "ctx_p25_associated_content_authors",
            25,
            "Associated Content, author information, funding and notes",
            "administrative_metadata",
            "Supporting-information and data-repository items are listed, followed by author/funding metadata. The listed SI was not supplied or explicitly associated in source_mapping.json and was not imported.",
        ),
        context(
            "ctx_p26_acknowledgments_abbreviations",
            26,
            "Acknowledgments and abbreviations",
            "administrative_metadata",
            "Acknowledgments and abbreviation definitions only; no current-paper reaction arrow, table row, or experimental run is reported.",
        ),
        context(
            "ctx_p26_references",
            26,
            "References 1-30",
            "cited_literature",
            "Cited literature only. Reference-list footnote 29 gives reagent-handling advice for a cited/currently discussed transformation, but it is not a new independently specified run on page 26.",
        ),
        context(
            "ctx_p27_references",
            27,
            "References 30-67",
            "cited_literature",
            "Cited literature only; no current-paper experimental record is created from bibliographic entries.",
        ),
        context(
            "ctx_p28_references",
            28,
            "References 68-70",
            "cited_literature",
            "Cited literature only; no current-paper experimental record is created from bibliographic entries.",
        ),
    ]
)

raw_extraction = {
    "batch_id": "allene-stereo-v31-200-20261017",
    "document_id": DOC,
    "source_kind": "original_article_pdf",
    "presentation_conversion": False,
    "si_association": "unestablished",
    "pages": [
        {
            "page": 25,
            "observations": [
                "ALARM-MSPS mass-spectrometric acquisition and peptide-identification method continues from page 24.",
                "C232 peptide behavior is interpreted as less reactive than C245; FSGDLDDQTC245R is selected for the main analysis.",
                "Associated Content, authors, funding, and notes follow; no reaction arrow or new experimental run appears.",
            ],
            "decision": "Update five existing ALARM-MSPS records and add source context; create no new reaction or molecule.",
            "linked_record_ids": REACTION_IDS,
            "evidence_ids": evidence(25),
        },
        {
            "page": 26,
            "observations": ["Acknowledgments, abbreviations, and references 1-30."],
            "decision": "Administrative/cited-literature context only; create no experiment.",
            "linked_record_ids": [],
            "evidence_ids": evidence(26),
        },
        {
            "page": 27,
            "observations": ["References 30-67."],
            "decision": "Cited-literature context only; create no experiment.",
            "linked_record_ids": [],
            "evidence_ids": evidence(27),
        },
        {
            "page": 28,
            "observations": ["References 68-70."],
            "decision": "Cited-literature context only; create no experiment.",
            "linked_record_ids": [],
            "evidence_ids": evidence(28),
        },
    ],
}

coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {
            "page": 25,
            "state": "attempted",
            "regions": [
                {
                    "region": "top-left continuation of La Antigen MS Analysis",
                    "finding": "ALARM-MSPS acquisition and peptide-identification criteria complete the five existing test/control records; analytical method, not a new synthesis.",
                    "evidence_ids": evidence(25),
                },
                {
                    "region": "left column C232 peptide interpretation",
                    "finding": "C232 peptides are described as less reactive than C245, and FSGDLDDQTC245R is selected for the main analysis; mechanistic/characterization context without an isolated adduct product.",
                    "evidence_ids": evidence(25),
                },
                {
                    "region": "Associated Content, author information, funding and notes",
                    "finding": "Administrative and repository/SI listings; no new current-paper experiment. Unassociated SI was not imported.",
                    "evidence_ids": evidence(25),
                },
            ],
            "record_ids": REACTION_IDS,
            "unresolved": [
                "Page 25 does not resolve the Figure-6 50 mM versus page-24 dilution-derived 0.5 mM test-compound concentration conflict.",
                "Listed supporting-information files were not supplied or explicitly associated and were not used.",
                "No discrete C232 or C245 compound-adduct structure is reported on this page.",
            ],
        },
        {
            "page": 26,
            "state": "attempted",
            "regions": [
                {
                    "region": "Acknowledgments and abbreviations",
                    "finding": "Administrative definitions and acknowledgments; no reaction arrow, table row, or experimental run.",
                    "evidence_ids": evidence(26),
                },
                {
                    "region": "References 1-30",
                    "finding": "Cited prior literature, including a reference-list reagent-handling note; not promoted to current-paper experiments.",
                    "evidence_ids": evidence(26),
                },
            ],
            "record_ids": [],
            "unresolved": [
                "Cited works were not externally inspected and do not establish new current-paper records on page 26."
            ],
        },
        {
            "page": 27,
            "state": "attempted",
            "regions": [
                {
                    "region": "References 30-67",
                    "finding": "Cited prior literature only; no current-paper reaction arrow, table row, or experimental run.",
                    "evidence_ids": evidence(27),
                }
            ],
            "record_ids": [],
            "unresolved": [
                "Cited works were not externally inspected and were not treated as experiments in this document."
            ],
        },
        {
            "page": 28,
            "state": "attempted",
            "regions": [
                {
                    "region": "References 68-70",
                    "finding": "Cited prior literature only; no current-paper reaction arrow, table row, or experimental run.",
                    "evidence_ids": evidence(28),
                }
            ],
            "record_ids": [],
            "unresolved": [
                "Cited works were not externally inspected and were not treated as experiments in this document."
            ],
        },
    ],
}

(RUN / "submitted.json").write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")
(RUN / "raw_model" / "pages_0025_0028_raw_extraction.json").write_text(
    json.dumps(raw_extraction, ensure_ascii=False, indent=2) + "\n"
)
(RUN / "coverage" / "pages_0025_0028.json").write_text(
    json.dumps(coverage, ensure_ascii=False, indent=2) + "\n"
)

print(
    json.dumps(
        {
            "updated_reactions": len(REACTION_IDS),
            "new_reactions": 0,
            "new_molecules": 0,
            "new_source_context": 7,
            "coverage_pages": [25, 26, 27, 28],
        },
        ensure_ascii=False,
    )
)
