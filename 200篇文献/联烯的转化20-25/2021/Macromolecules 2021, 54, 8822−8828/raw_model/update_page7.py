#!/usr/bin/env python3
"""Merge source-derived evidence from original PDF page 7 into the v3.1 bundle.

Page 7 contains bibliography entries and two current-work footnotes. It does
not contain a new reaction arrow, table row, or independent synthesis.
"""
from pathlib import Path
import json


ROOT = Path(__file__).resolve().parents[3]
RUN = ROOT / "runs" / "current"
DOC = "doc_56ca04973542741c"
TEXT = f"{DOC}:page7:text"
IMAGE = f"{DOC}:page7:image"
EVIDENCE = [TEXT, IMAGE]


def append_unique(values, additions):
    for value in additions:
        if value not in values:
            values.append(value)


def extend_fact_evidence(fact):
    append_unique(fact["evidence_ids"], EVIDENCE)


def note(raw, normalized):
    return {
        "raw_value": raw,
        "normalized_value": normalized,
        "acquisition_method": "original_extraction",
        "evidence_ids": EVIDENCE.copy(),
        "data_status": "pending_review",
        "verification": "unreviewed",
        "derivation": None,
        "alternatives": [],
    }


bundle_path = RUN / "submitted.json"
bundle = json.loads(bundle_path.read_text())
# The page-7 image evidence was staged after the earlier partial-page bundle.
# Reuse the untouched document block produced by `standard init` so merge can
# verify it exactly against dataset.json; chemistry records below remain intact.
init_bundle = json.loads((RUN / "raw_model" / "page7_init.json").read_text())
bundle["documents"] = init_bundle["documents"]
reactions = {reaction["reaction_id"]: reaction for reaction in bundle["reactions"]}

# Footnote 45 is explicitly about Table 1 entry 8 and confirms the isolated
# yield plus measured/target Mn comparison. It is duplicate evidence, not a
# new experiment.
entry8 = reactions[f"{DOC}:table1:entry8"]
append_unique(entry8["evidence_ids"], EVIDENCE)
extend_fact_evidence(entry8["location"])
extend_fact_evidence(entry8["experiment_categories"])
entry8["location"]["raw_value"] = "PDF pages 2, 7"
entry8["location"]["normalized_value"] = {"evidence_ids": entry8["evidence_ids"].copy()}
for outcome in entry8["outcomes"]:
    if outcome["outcome_id"] in {"yield", "target_mn", "measured_mn"}:
        extend_fact_evidence(outcome["value"])
        if outcome["outcome_id"] in {"yield", "measured_mn"}:
            extend_fact_evidence(outcome["measurement_basis"])
        if outcome["outcome_id"] == "yield":
            extend_fact_evidence(outcome["yield_type"])

# Footnote 43 reports shared characterization context for the block studies
# in which monomer 1 was polymerized first. The source does not identify the
# individual block experiment, so retain the observation on both linked
# records with an explicit unresolved limitation.
footnote43_raw = (
    "This sharp resonance was also observed during our block copolymerization "
    "studies (where 1 was polymerized first)."
)
footnote43_normalized = (
    "The authors attribute the sharp 11B NMR shoulder signal to residual "
    "monomer 1 and state that protodeboration is unlikely to be a major pathway."
)
shared_unresolved = (
    "Footnote 43 does not assign the shared NMR observation to block "
    "copolymer 6 versus 7 individually; it is retained as shared contextual evidence."
)
for reaction_id in (f"{DOC}:block:6", f"{DOC}:block:7"):
    reaction = reactions[reaction_id]
    append_unique(reaction["evidence_ids"], EVIDENCE)
    extend_fact_evidence(reaction["location"])
    reaction["location"]["raw_value"] = "PDF pages 4, 5, 7"
    reaction["location"]["normalized_value"] = {"evidence_ids": reaction["evidence_ids"].copy()}
    if not any(item["raw_value"] == footnote43_raw for item in reaction["notes"]):
        reaction["notes"].append(note(footnote43_raw, footnote43_normalized))
    if shared_unresolved not in reaction["review"]["unresolved"]:
        reaction["review"]["unresolved"].append(shared_unresolved)

context_id = f"{DOC}:page7"
page7_context = {
    "context_id": context_id,
    "page": 7,
    "kind": "references_and_current_work_footnotes",
    "evidence_ids": EVIDENCE.copy(),
    "finding": (
        "References 33-42, 44, and 46-52 are cited prior literature and were not "
        "imported as experiments of this article. Footnote 43 is shared NMR/"
        "protodeboration interpretation for polymer 2 and block studies, not a "
        "new synthesis. Footnote 45 duplicates and clarifies Table 1 entry 8."
    ),
}
bundle["source_context"] = [
    context for context in bundle["source_context"] if context.get("context_id") != context_id
] + [page7_context]
bundle_path.write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")

raw_page = {
    "document_id": DOC,
    "claimed_page_scope": [7],
    "semantic_accuracy_verified": False,
    "page_image_opened_at_original_detail": True,
    "title": "Vinyl-Addition Polymerizations of Borylated Allenes",
    "regions": [
        {
            "region": "left column, references 33-42",
            "finding": "Cited prior literature only; no current-paper experiment created.",
            "evidence_ids": EVIDENCE.copy(),
            "record_ids": [],
        },
        {
            "region": "left column, footnote 43",
            "finding": (
                "Characterization/author interpretation: a sharp 11B NMR signal is "
                "assigned to residual 1; protodeboration is considered unlikely."
            ),
            "evidence_ids": EVIDENCE.copy(),
            "record_ids": [f"{DOC}:block:6", f"{DOC}:block:7"],
        },
        {
            "region": "left column, footnote 45",
            "finding": (
                "Duplicate source for Table 1 entry 8: multigram-scale polymerization, "
                "96% isolated yield, measured Mn 39,970 versus target 30,000 Da."
            ),
            "evidence_ids": EVIDENCE.copy(),
            "record_ids": [f"{DOC}:table1:entry8"],
        },
        {
            "region": "remaining left column and upper right column, references 44 and 46-52",
            "finding": "Cited prior literature only; no current-paper experiment created.",
            "evidence_ids": EVIDENCE.copy(),
            "record_ids": [],
        },
    ],
    "new_reaction_ids": [],
    "updated_reaction_ids": [
        f"{DOC}:table1:entry8",
        f"{DOC}:block:6",
        f"{DOC}:block:7",
    ],
    "new_molecule_ids": [],
    "notes": [
        "No reaction arrows, tables, or standalone synthesis procedures occur on page 7.",
        "Footnote 43 is not sufficiently experiment-specific to create a new reaction record.",
        "All machine records remain needs_review; no SI was associated.",
    ],
}
(RUN / "raw_model" / "extracted_page_0007.json").write_text(
    json.dumps(raw_page, ensure_ascii=False, indent=2) + "\n"
)

coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {
            "page": 7,
            "state": "attempted",
            "regions": [
                {
                    "region": "left column, references 33-42",
                    "finding": "Cited prior work only; no current-paper reaction record was created.",
                    "evidence_ids": EVIDENCE.copy(),
                },
                {
                    "region": "left column, footnote 43",
                    "finding": (
                        "Characterization/author interpretation rather than an independent "
                        "synthesis: the sharp 11B NMR signal is assigned to residual 1, and "
                        "protodeboration is considered unlikely. It is shared across unspecified "
                        "block copolymerization studies in which 1 was polymerized first."
                    ),
                    "evidence_ids": EVIDENCE.copy(),
                },
                {
                    "region": "left column, footnote 45",
                    "finding": (
                        "Explicit duplicate evidence for existing Table 1 entry 8: multigram-scale "
                        "polymerization of 1, 96% isolated polymer 2, and GPC Mn 39,970 versus "
                        "30,000 Da target."
                    ),
                    "evidence_ids": EVIDENCE.copy(),
                },
                {
                    "region": "remaining left column and upper right column, references 44 and 46-52",
                    "finding": "Cited prior work only; no current-paper reaction record was created.",
                    "evidence_ids": EVIDENCE.copy(),
                },
            ],
            "record_ids": [
                f"{DOC}:table1:entry8",
                f"{DOC}:block:6",
                f"{DOC}:block:7",
            ],
            "unresolved": [
                "Footnote 43 does not identify whether its shared sharp-signal observation applies to block copolymer 6, 7, or both individually.",
                "Figures S4 and S21 and the supporting experimental details are in unavailable, unassociated SI.",
            ],
        }
    ],
}
(RUN / "coverage" / "pages_0007_0007.json").write_text(
    json.dumps(coverage, ensure_ascii=False, indent=2) + "\n"
)

print(json.dumps({"new_reactions": 0, "updated_reactions": 3, "page": 7}))
