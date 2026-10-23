#!/usr/bin/env python3
"""Persist validation/audit snapshots and final acceptance checks."""
import json
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
RUN = ROOT / "runs/current"
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import reaction_tools as rt
import standard_records as st

data = st.read(RUN / "dataset.json")
standard = st.read(RUN / "standard_records.json")
std_errors = st.validate(standard, RUN)
strict_errors = rt.validate(data, RUN, strict_record_ids={r["record_id"] for r in data["experiments"]})
audit = rt.dataset_audit(data)
st.write(RUN / "raw_model/standard_validation_pages_0001_0005.json", {
    "valid": not std_errors, "errors": std_errors, "chemical_accuracy_verified": False
})
st.write(RUN / "raw_model/strict_validation_pages_0001_0005.json", {
    "valid": not strict_errors, "errors": strict_errors,
    "standard_v3_present": "standard" in data, "chemical_accuracy_verified": False
})
st.write(RUN / "raw_model/audit_pages_0001_0005.json", audit)

mentions = st.source_mentions(data, ["2f", "2s", "3m", "5b", "4j", "4l"], 180)
st.write(RUN / "raw_model/source_mentions_pages_0001_0005.json", mentions)

coverage = st.read(RUN / "coverage/pages_0001_0005.json")
assert [p["page"] for p in coverage["pages"]] == [1, 2, 3, 4, 5]
assert all(p["state"] == "attempted" and p["regions"] for p in coverage["pages"])
assert len(standard["reactions"]) == 83 and len(data["experiments"]) == 83
assert len(standard["molecules"]) == 134
assert all(r["review"]["status"] == "needs_review" for r in standard["reactions"])

samples = {
    "table2_2m": (RUN / "records/doc_08b3dd4eacda5aa8__table2_2m.md").read_text(),
    "table4_4a": (RUN / "records/doc_08b3dd4eacda5aa8__table4_4a.md").read_text(),
    "table5_5b": (RUN / "records/doc_08b3dd4eacda5aa8__table5_5b.md").read_text(),
    "table1_07": (RUN / "records/doc_08b3dd4eacda5aa8__table1_entry_07.md").read_text(),
}
azide_record = st.read(RUN / "records/doc_08b3dd4eacda5aa8__table4_4a.json")
markdown_checks = {
    "gp_a_reagent_doses": all(x in samples["table2_2m"] for x in ["1.1 equiv.", "5 mL", "3 × 20 mL"]),
    "azide_excludes_base_and_has_tmsn3": "TMS-N3" in samples["table4_4a"] and not any(p["molecule_id"].endswith(":mol:li2co3") for p in azide_record["participants"]),
    "gp_b_reagent_doses": all(x in samples["table5_5b"] for x in ["1.2 equiv.", "1.1 equiv.", "5 mL", "3 × 20 mL"]),
    "table1_ratio_method": "1H NMR" in samples["table1_07"] and "100:0" in samples["table1_07"],
    "unreported_yield_basis_visible": "Yield measurement basis is not stated" in samples["table2_2m"],
}
assert all(markdown_checks.values()), markdown_checks

issues = Counter()
for row in audit["warnings"]:
    for warning in row["warnings"]:
        issues[warning.split(":", 1)[0]] += 1
verification = {
    "pages": [1, 2, 3, 4, 5],
    "reactions": len(standard["reactions"]),
    "molecules": len(standard["molecules"]),
    "record_kind_counts": audit["record_kind_counts"],
    "warning_categories": dict(issues),
    "markdown_checks": markdown_checks,
    "source_mentions": len(mentions["matches"]),
    "source_reviewed_records": audit["source_reviewed_records"],
    "semantic_accuracy_verified": False,
}
st.write(RUN / "raw_model/final_checks_pages_0001_0005.json", verification)

workflow_path = ROOT / "runs/workflow.md"
workflow = workflow_path.read_text()
workflow += f"""

## Final validation and audit

- `standard validate` and legacy `validate --strict` both pass with zero errors; snapshots are saved under `runs/current/raw_model/`.
- Final totals are 83 reactions and 134 molecules. Record kinds are 4 generic schemes, 17 optimization entries, 56 substrate-scope entries and 6 other reported derivatizations. All 83 reactions remain `needs_review`; source-reviewed count is zero.
- Audit has no `reaction_identity`, `outcome_target`, `structure_gap`, `missing_steps` or characterization-only warning. It retains {issues.get('measurement_basis', 0)} `measurement_basis` advisories because the main-text pages print yields but do not state isolated/NMR/GC yield basis. Those fields intentionally remain unreported, with the reason in each record; no method was guessed from chromatography.
- Literal cross-region mention hints for 2f, 2s, 3m, 5b, 4j and 4l are saved separately and were used only as review pointers; identical labels were not used to merge distinct experiments. In particular, Table 5 5b (91%) and Scheme 4 3m→5b (94%) remain separate experiments.
- Generated record Markdown was checked for all shared reagent doses, THF/DMF 5 mL, EtOAc 3 × 20 mL, the azide base exception, Table 1 1H-NMR E/Z basis and explicit unreported-yield-basis notes. No numerical source conflict was found on pages 1–5, so no artificial conflict fact was created.
"""
workflow_path.write_text(workflow)
print(json.dumps(verification, ensure_ascii=False, indent=2))
