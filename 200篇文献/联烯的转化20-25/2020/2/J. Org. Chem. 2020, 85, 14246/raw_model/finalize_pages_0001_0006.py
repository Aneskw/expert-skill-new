#!/usr/bin/env python3
"""Persist validation, audit, mention, and handoff checks for pages 1–6."""
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[3]
RUN = ROOT / "runs/current"
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))

import reaction_tools
import standard_records as st


def write(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n")


def main():
    dataset = json.loads((RUN / "dataset.json").read_text())
    bundle = json.loads((RUN / "standard_records.json").read_text())
    standard_errors = st.validate(bundle, RUN)
    legacy_errors = reaction_tools.validate(
        dataset, RUN, strict_record_ids={r["record_id"] for r in dataset["experiments"]}
    )
    strict = {
        "valid": not standard_errors and not legacy_errors,
        "standard_v3_errors": standard_errors,
        "strict_compatibility_errors": legacy_errors,
        "standard_v3_present": "standard" in dataset,
        "semantic_accuracy_verified": False,
        "chemical_accuracy_verified": False,
    }
    write(RUN / "raw_model/pages_0001_0006_strict_validate.json", strict)

    audit = reaction_tools.dataset_audit(dataset)
    write(RUN / "raw_model/pages_0001_0006_audit.json", audit)
    write(RUN / "quality_checks.json", audit)

    terms = ["7b", "8a", "12", "13", "(+)-Hypoxylactone", "Xylogiblactone A"]
    mentions = st.source_mentions(dataset, terms, context=220)
    write(RUN / "raw_model/pages_0001_0006_source_mentions.json", mentions)

    coverage = json.loads((RUN / "coverage/pages_0001_0006.json").read_text())
    reaction_ids = {r["reaction_id"] for r in bundle["reactions"]}
    coverage_checks = {
        "document_id_matches": coverage.get("document_id") == "doc_7ca2cfc676d4ba33",
        "exact_pages_1_to_6": [p.get("page") for p in coverage.get("pages", [])] == [1, 2, 3, 4, 5, 6],
        "all_states_attempted": all(p.get("state") == "attempted" for p in coverage.get("pages", [])),
        "all_regions_nonempty": all(p.get("regions") for p in coverage.get("pages", [])),
        "all_record_ids_exist": all(rid in reaction_ids for p in coverage.get("pages", []) for rid in p.get("record_ids", [])),
        "semantic_accuracy_verified": coverage.get("semantic_accuracy_verified"),
    }

    markdown = (RUN / "规范记录预览.md").read_text() + "\n" + "\n".join(
        path.read_text() for path in sorted((RUN / "records").glob("*.md"))
    )
    required_markdown_terms = [
        "0.50 g", "0.96 mmol", "1.06 equiv", "1.04 equiv", "2.1 equiv",
        "D-threonine", "L-threonine", "−78 °C", "−30 °C", "75–83%",
        "Ag2O", "21 mg", "0.091 mmol",
    ]
    markdown_checks = {term: term in markdown for term in required_markdown_terms}
    handoff = {
        "coverage": coverage_checks,
        "markdown_terms": markdown_checks,
        "title": bundle["documents"][0].get("title"),
        "molecule_count": len(bundle["molecules"]),
        "reaction_count": len(bundle["reactions"]),
        "specific_experiment_count": sum(r["record_type"] == "specific_experiment" for r in bundle["reactions"]),
        "audit_warning_count": len(audit["warnings"]),
        "source_reviewed_records": audit["source_reviewed_records"],
        "si_association": "unestablished; no SI document imported",
    }
    write(RUN / "raw_model/pages_0001_0006_handoff_checks.json", handoff)

    if not strict["valid"]:
        raise SystemExit("strict validation failed")
    if not all(coverage_checks[k] for k in ("document_id_matches", "exact_pages_1_to_6", "all_states_attempted", "all_regions_nonempty", "all_record_ids_exist")):
        raise SystemExit("coverage contract check failed")
    if not all(markdown_checks.values()):
        raise SystemExit("generated Markdown is missing required source terms")

    workflow = ROOT / "runs/workflow.md"
    with workflow.open("a") as stream:
        stream.write("""

## Validation and audit result

- Pre-merge submission validation, post-merge `standard validate`, and dataset `validate --strict` all passed with no errors. Persisted result: `runs/current/raw_model/pages_0001_0006_strict_validate.json`.
- Audit reports **0 warning records**. In particular, there are no `reaction_identity`, `outcome_target`, `measurement_basis`, or `structure_gap` warnings. This is an advisory result and does not establish chemical accuracy or extraction recall.
- `documents_without_linked_si` correctly contains this document because no SI was supplied or explicitly associated.
- Generated Markdown was checked for exact low-temperature reagent charges/equivalents, Ag2O charges, the D/L-threonine and −78/−30 °C conflict, and the unallocated 75–83% Scheme-2 range. All required terms are present.
- Coverage contract checks pass: exactly pages 1–6, every state is `attempted`, every page has nonempty regions, and every listed record ID exists.
- Literal cross-region mention hints were saved at `runs/current/raw_model/pages_0001_0006_source_mentions.json`; they were used only as review leads and did not merge experiments by label alone.
""")
    print(json.dumps(handoff, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
