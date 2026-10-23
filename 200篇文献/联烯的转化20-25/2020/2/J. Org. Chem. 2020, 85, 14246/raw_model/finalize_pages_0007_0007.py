#!/usr/bin/env python3
"""Persist strict validation, audit, and handoff checks for PDF page 7."""
from __future__ import annotations

import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[3]
RUN = ROOT / "runs/current"
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))

import reaction_tools
import standard_records as st


DOC = "doc_7ca2cfc676d4ba33"
CONTEXT_ID = f"{DOC}:context_page7_notes_acknowledgments_references"


def write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n")


def main() -> None:
    dataset = json.loads((RUN / "dataset.json").read_text())
    bundle = json.loads((RUN / "standard_records.json").read_text())
    standard_errors = st.validate(bundle, RUN)
    legacy_errors = reaction_tools.validate(
        dataset,
        RUN,
        strict_record_ids={record["record_id"] for record in dataset["experiments"]},
    )
    strict = {
        "valid": not standard_errors and not legacy_errors,
        "standard_v3_errors": standard_errors,
        "strict_compatibility_errors": legacy_errors,
        "standard_v3_present": "standard" in dataset,
        "semantic_accuracy_verified": False,
        "chemical_accuracy_verified": False,
    }
    write_json(RUN / "raw_model/pages_0007_0007_strict_validate.json", strict)

    audit = reaction_tools.dataset_audit(dataset)
    write_json(RUN / "raw_model/pages_0007_0007_audit.json", audit)
    write_json(RUN / "quality_checks.json", audit)

    coverage = json.loads((RUN / "coverage/pages_0007_0007.json").read_text())
    page = coverage.get("pages", [{}])[0]
    context_ids = {item.get("context_id") for item in bundle.get("source_context", [])}
    all_needs_review = all(
        item.get("review", {}).get("status") == "needs_review"
        for item in bundle.get("reactions", [])
    )
    checks = {
        "document_id_matches": coverage.get("document_id") == DOC,
        "exact_page_7_only": [item.get("page") for item in coverage.get("pages", [])] == [7],
        "page_state_attempted": page.get("state") == "attempted",
        "regions_nonempty": bool(page.get("regions")),
        "record_ids_empty": page.get("record_ids") == [],
        "page7_context_present": CONTEXT_ID in context_ids,
        "all_machine_reactions_remain_needs_review": all_needs_review,
        "semantic_accuracy_verified": coverage.get("semantic_accuracy_verified"),
    }
    relevant_warnings = [
        warning
        for warning in audit.get("warnings", [])
        if warning.get("kind") in {
            "reaction_identity",
            "outcome_target",
            "measurement_basis",
            "structure_gap",
        }
    ]
    handoff = {
        "batch_id": "allene-stereo-v31-200-20261017",
        "document_id": DOC,
        "page_range": [7, 7],
        "increment_counts": {"molecules": 0, "reactions": 0, "source_context_items": 1},
        "bundle_counts_after_merge": {
            "molecules": len(bundle.get("molecules", [])),
            "reactions": len(bundle.get("reactions", [])),
            "source_context_items": len(bundle.get("source_context", [])),
        },
        "coverage_checks": checks,
        "audit_warning_count": len(audit.get("warnings", [])),
        "requested_audit_warning_count": len(relevant_warnings),
        "requested_audit_warnings": relevant_warnings,
        "si_association": "unestablished; no SI document imported",
        "unresolved": page.get("unresolved", []),
    }
    write_json(RUN / "raw_model/pages_0007_0007_handoff_checks.json", handoff)

    if not strict["valid"]:
        raise SystemExit("strict validation failed")
    required_checks = [
        "document_id_matches",
        "exact_page_7_only",
        "page_state_attempted",
        "regions_nonempty",
        "record_ids_empty",
        "page7_context_present",
        "all_machine_reactions_remain_needs_review",
    ]
    if not all(checks[name] for name in required_checks):
        raise SystemExit(f"page-7 handoff check failed: {checks}")

    workflow = ROOT / "runs/workflow.md"
    marker = "## Page-7 validation and audit result"
    text = workflow.read_text()
    if marker not in text:
        text += f"""

{marker}

- `standard merge`, post-merge v3 validation, and legacy `validate --strict` passed with no errors; result: `runs/current/raw_model/pages_0007_0007_strict_validate.json`.
- Page-7 coverage contract passed: exactly page 7, state `attempted`, two nonempty located regions, and an intentionally empty `record_ids` list because no experiment occurs on the page.
- Audit produced **{len(audit.get('warnings', []))} total warnings** and **{len(relevant_warnings)} warnings** in the requested `reaction_identity`, `outcome_target`, `measurement_basis`, and `structure_gap` groups. The audit is advisory and does not establish chemistry correctness or recall.
- Formal bundle totals remain **{len(bundle.get('molecules', []))} molecules** and **{len(bundle.get('reactions', []))} reactions**; the only page-7 addition is one bibliography/back-matter source-context item. No prior experiment was re-counted as a new page-7 result.
"""
        workflow.write_text(text)

    print(json.dumps(handoff, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
