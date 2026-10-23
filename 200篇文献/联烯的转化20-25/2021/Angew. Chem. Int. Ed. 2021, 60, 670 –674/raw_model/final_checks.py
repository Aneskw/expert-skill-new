#!/usr/bin/env python3
import json
from pathlib import Path

RUN = Path(__file__).resolve().parents[1]
data = json.loads((RUN / "dataset.json").read_text())
standard = json.loads((RUN / "standard_records.json").read_text())
coverage = json.loads((RUN / "coverage/pages_0001_0005.json").read_text())

assert data["standard"] == standard
assert [p["page"] for p in coverage["pages"]] == [1, 2, 3, 4, 5]
assert all(p["state"] == "attempted" and p["regions"] for p in coverage["pages"])
assert len(standard["reactions"]) == 76
assert len(standard["molecules"]) == 84
assert len(standard["source_context"]) == 5
assert all(r["review"]["status"] == "needs_review" for r in standard["reactions"])
assert all(
    m["stereochemistry"][key]["normalized_value"] is not None
    for m in standard["molecules"]
    for key in ("chirality_status", "geometry_status", "representation_status")
)

checks = {
    "table1_entry12": {
        "file": "records/doc_151fad1677f63120__table1_entry12.md",
        "tokens": ["0.2 mmol", "0.25 mmol", "1 mol%", "1.1 equiv", "0.1 equiv", "1.0 mL", "room temperature", "8 h", "83%", "NMR spectroscopy using anisole", ">20:1 d.r.", "cross-page duplicate merged"],
    },
    "table1_entry10": {
        "file": "records/doc_151fad1677f63120__table1_entry10.md",
        "tokens": ["Pd-AmP-MCF", "5 mol%", "AgOTf", "DCE", "12 h", "94%", "NMR spectroscopy using anisole"],
    },
    "scheme3_5a_mcf": {
        "file": "records/doc_151fad1677f63120__scheme3_5a_PdAmPMCF_1mol.md",
        "tokens": ["1.25 equiv", "1.1 equiv", "1 mol%", "CHCl3", "8 h", "59%", ">90% conversion", "59:40"],
    },
    "scheme4b": {
        "file": "records/doc_151fad1677f63120__scheme4b_borylation_6.md",
        "tokens": ["B2pin2", "1.25 equiv", "1 mol%", "1.1 equiv", "0.1 equiv", "DCE", "room temperature", "8 h", "89%", ">20:1 d.r."],
    },
    "footnote38": {
        "file": "records/doc_151fad1677f63120__footnote38_aromatic_R_1i_failure.md",
        "tokens": ["messy mixture", "TLC and 1H NMR spectroscopy", "terminal alkyne identity", "exact conditions"],
    },
}
for item in checks.values():
    text = (RUN / item["file"]).read_text()
    item["missing_tokens"] = [token for token in item["tokens"] if token not in text]
    assert not item["missing_tokens"], (item["file"], item["missing_tokens"])

summary = {
    "documents_equal": True,
    "pages": [1, 2, 3, 4, 5],
    "reaction_count": len(standard["reactions"]),
    "molecule_count": len(standard["molecules"]),
    "source_context_count": len(standard["source_context"]),
    "all_machine_records_need_review": True,
    "all_molecules_have_explicit_stereo_statuses": True,
    "markdown_checks": checks,
    "semantic_accuracy_verified": False,
    "chemical_accuracy_verified": False,
}
(RUN / "raw_model/verification_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n")
print(json.dumps(summary, ensure_ascii=False, indent=2))
