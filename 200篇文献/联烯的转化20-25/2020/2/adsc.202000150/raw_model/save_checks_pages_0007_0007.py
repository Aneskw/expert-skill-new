#!/usr/bin/env python3
"""Save validation, audit, and page-7 evidence checks for the bounded pass."""
from __future__ import annotations

import json
import sys
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))

import reaction_tools  # noqa: E402
import standard_records  # noqa: E402


RUN = ROOT / "runs/current"
DOC = "doc_65d1cdf1cb09f398"
PAGE7 = {f"{DOC}:page7:image", f"{DOC}:page7:text"}
TARGETS = {f"{DOC}:table2_3h", f"{DOC}:table3_5c"}


def write(name: str, payload) -> None:
    (RUN / "raw_model" / name).write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n"
    )


dataset = reaction_tools.read(RUN / "dataset.json")
standard = standard_records.read(RUN / "standard_records.json")
coverage = reaction_tools.read(RUN / "coverage/pages_0007_0007.json")

standard_errors = standard_records.validate(standard, RUN)
strict_errors = reaction_tools.validate(
    dataset,
    RUN,
    strict_record_ids={row["record_id"] for row in dataset["experiments"]},
)
audit = reaction_tools.dataset_audit(dataset)

reactions = {row["reaction_id"]: row for row in standard["reactions"]}
molecules = {row["molecule_id"]: row for row in standard["molecules"]}
target_reactions = {rid: reactions[rid] for rid in TARGETS}

evidence_checks = {
    "reaction_top_level": {
        rid: PAGE7.issubset(set(row["evidence_ids"]))
        for rid, row in target_reactions.items()
    },
    "reaction_ccdc_note": {
        rid: any(
            PAGE7.issubset(set(note["evidence_ids"]))
            and isinstance(note.get("normalized_value"), dict)
            and note["normalized_value"].get("CCDC") in {"1937754", "1968335"}
            for note in row["notes"]
        )
        for rid, row in target_reactions.items()
    },
    "molecule_assignment": {
        mid: PAGE7.issubset(
            set(molecules[mid]["stereochemistry"]["assignment"]["evidence_ids"])
        )
        for mid in ("m_3h", "m_5c")
    },
    "markdown": {},
}
for rid in TARGETS:
    markdown_path = RUN / "records" / (rid.replace(":", "__") + ".md")
    markdown = markdown_path.read_text()
    evidence_checks["markdown"][rid] = (
        "page7:image" in markdown
        and "page7:text" in markdown
        and "PDF 第 7 页" in markdown
    )

warning_types = Counter()
for row in audit["warnings"]:
    for warning in row["warnings"]:
        warning_types[warning.split(":", 1)[0]] += 1

all_review_statuses = Counter(row["review"]["status"] for row in standard["reactions"])
coverage_shape_valid = (
    coverage["document_id"] == DOC
    and coverage["semantic_accuracy_verified"] is False
    and len(coverage["pages"]) == 1
    and coverage["pages"][0]["page"] == 7
    and bool(coverage["pages"][0]["regions"])
)

write(
    "standard_validation_pages_0007_0007.json",
    {
        "valid": not standard_errors,
        "errors": standard_errors,
        "chemical_accuracy_verified": False,
    },
)
write(
    "strict_validation_pages_0007_0007.json",
    {
        "valid": not strict_errors,
        "errors": strict_errors,
        "standard_v3_present": "standard" in dataset,
        "chemical_accuracy_verified": False,
    },
)
write("audit_pages_0007_0007.json", audit)
write(
    "page_0007_consistency_checks.json",
    {
        "page_range": [7, 7],
        "coverage_shape_valid": coverage_shape_valid,
        "record_count": len(standard["reactions"]),
        "molecule_count": len(standard["molecules"]),
        "source_context_count": len(standard["source_context"]),
        "new_reaction_count": 0,
        "updated_reaction_ids": sorted(TARGETS),
        "evidence_checks": evidence_checks,
        "all_evidence_checks_pass": all(
            value
            for group in evidence_checks.values()
            for value in group.values()
        ),
        "review_status_counts": dict(all_review_statuses),
        "audit_warning_type_counts": dict(warning_types),
        "reaction_identity_warnings": warning_types["reaction_identity"],
        "outcome_target_warnings": warning_types["outcome_target"],
        "measurement_basis_warnings": warning_types["measurement_basis"],
        "structure_gap_warnings": warning_types["structure_gap"],
        "chemical_accuracy_verified": False,
    },
)

print(
    json.dumps(
        {
            "standard_valid": not standard_errors,
            "strict_valid": not strict_errors,
            "coverage_shape_valid": coverage_shape_valid,
            "all_evidence_checks_pass": all(
                value
                for group in evidence_checks.values()
                for value in group.values()
            ),
            "records": len(standard["reactions"]),
            "molecules": len(standard["molecules"]),
            "source_context": len(standard["source_context"]),
            "review_statuses": dict(all_review_statuses),
            "warning_types": dict(warning_types),
        },
        ensure_ascii=False,
    )
)
