#!/usr/bin/env python3
"""Write page-7 coverage, audit snapshot, task state and workflow decisions."""

from collections import Counter
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[3]
RUN = ROOT / "runs/current"
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import reaction_tools as rt


DOC = "doc_68d6021490685b69"
RID = f"{DOC}:footnote19_racemic_0C_2plus2_detection"
E_TEXT = f"{DOC}:page7:text"
E_IMAGE = f"{DOC}:page7:image"

data = json.loads((RUN / "dataset.json").read_text())
standard = data["standard"]
reactions = standard["reactions"]
assert any(r["reaction_id"] == RID for r in reactions)

audit = rt.dataset_audit(data)
rt.write(RUN / "audit.json", audit)
rt.write(RUN / "quality_checks.json", audit)

coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [{
        "page": 7,
        "state": "attempted",
        "regions": [
            {
                "region": "left column and upper right column, references 12–18",
                "finding": (
                    "Prior-literature/review citations only; no current-paper "
                    "reaction arrow, table row or characterization peak list."
                ),
                "evidence_ids": [E_IMAGE, E_TEXT],
            },
            {
                "region": "right column, footnote 19",
                "finding": (
                    "Distinct mechanistic control: the racemic Table 1 entry-1 "
                    "reaction was run at 0 °C and an unnumbered [2+2] "
                    "cycloaddition product was qualitatively detected.  One "
                    "partial needs_review experiment record was added."
                ),
                "evidence_ids": [E_IMAGE, E_TEXT],
            },
            {
                "region": "right column, footnote 20",
                "finding": (
                    "CCDC 2042864 (3aa), 2052253 (rac-4ka) and 2079315 "
                    "(rac-5aa) are crystallographic characterization context, "
                    "not new synthesis experiments."
                ),
                "evidence_ids": [E_IMAGE, E_TEXT],
            },
            {
                "region": "right column, footnote 21",
                "finding": (
                    "Qualitative statement that a one-pot procedure gave "
                    "product family 4 as major product.  It lacks a substrate "
                    "instance, yield, ratio and conditions on this page; it is "
                    "retained in source_context rather than fabricated as a "
                    "fully distinct reaction row."
                ),
                "evidence_ids": [E_IMAGE, E_TEXT],
            },
            {
                "region": "right column, references 22–23 and publication footer",
                "finding": (
                    "Prior-literature citations and administrative footer only; "
                    "no additional current-paper experiment."
                ),
                "evidence_ids": [E_IMAGE, E_TEXT],
            },
        ],
        "record_ids": [RID],
        "unresolved": [
            "Footnote 19 does not print the detected product number, structure, amount/yield, detection method or product stereochemical composition.",
            "Footnote 21 is not uniquely distinguishable as a concrete experiment without the cited ESI Table S5.",
            "The article cites ESI Tables S1 and S5, but no SI file is explicitly associated in source_mapping.json; no SI facts were inferred.",
        ],
    }],
}
rt.write(RUN / "coverage/pages_0007_0007.json", coverage)

tasks_path = RUN / "tasks.json"
tasks = json.loads(tasks_path.read_text())
page7_tasks = [t for t in tasks if t.get("evidence_id") == E_TEXT]
assert len(page7_tasks) == 1
page7_tasks[0]["state"] = "attempted"
page7_tasks[0]["reason"] = (
    "Original page-7 image inspected; references and footnotes are inventoried "
    "in coverage/pages_0007_0007.json."
)
rt.write(tasks_path, tasks)

new_warnings = [w for w in audit["warnings"] if w["record_id"] == RID]
issue_counts = Counter(
    warning.split(":", 1)[0]
    for item in new_warnings
    for warning in item["warnings"]
)

marker = "## Page 7 extension (current session)"
workflow_path = ROOT / "runs/workflow.md"
existing = workflow_path.read_text() if workflow_path.exists() else "# Extraction workflow\n"
if marker in existing:
    existing = existing.split(marker, 1)[0].rstrip() + "\n\n"

section = f"""{marker}

- Batch: `allene-stereo-v31-200-20261017`
- Document: `{DOC}`
- Actual visual scope: original PDF page 7 only (`7–7`, one page).
- Source mapping: main article only; `supporting_information_association` remains `unestablished`.
- Result: one new concrete experiment record, `{RID}`. The merged bundle now contains {len(standard['molecules'])} molecules/agents and {len(reactions)} reactions; all remain `needs_review`.

### Page inventory and decisions

1. Read the complete extraction skill and the required extraction, stereochemistry, structured-fields, quality, schema and bond-taxonomy references before editing records.
2. Inspected `runs/current/assets/{DOC}/page7.jpg` at original resolution and a lineage-tracked crop `runs/current/raw_model/page7_notes19_23.png`.
3. References 12–18 and 22–23 are prior literature only. They remain bibliographic source context and are not counted as this paper's experiments.
4. Footnote 19 is a distinct low-temperature mechanistic control: the racemic Table 1 entry-1 experiment was run at 0 °C and an unnumbered [2+2] product was detected. The explicit Table 1 link supports inheritance of entry-1 scale, reactant ratio, Ni(OTf)2 loading, CH2Cl2 concentration, 16 h and N2; 0 °C and qualitative detection come from page 7. No ESI detail is imported.
5. The detected product has a separate molecule ID rather than being merged into enantioenriched `(1R,2R,4S)-3aa`. The author's phrase “racemic reaction” is retained, but no product ee/er, 1:1 composition or absolute configuration is invented. No product SMILES is entered; representation remains `not_encoded`, chirality/geometry remain `unresolved`, and `standard stereo` is not applicable to that unprinted structure.
6. The reused concrete structures 1a and 2a were rechecked with the standard stereo diagnostic; the derived output is saved in `runs/current/raw_model/page7_reused_stereo_diagnostics.json`. This graph check does not certify source identity.
7. Footnote 20 contains only CCDC characterization identifiers. Footnote 21 says an under-specified one-pot procedure gave product family 4 as major product, but gives no instance, yield, ratio or conditions and points only to unavailable ESI Table S5. Both are retained in `source_context`; footnote 21 is not promoted to a fabricated distinct reaction.

### Record and audit decisions

- Footnote 19 is classified partially as two new C–C connections characteristic of the explicitly named [2+2] product. No atom-mapped reaction SMILES or complete classification is claimed.
- Its outcome is qualitative `detected`, targeted to the separate detected-product molecule. Yield, amount, detection method and measurement basis remain null/unreported.
- The page contains no independent mass-spectrum peak table, molecular-weight table, reaction scheme or optimization/scope table.
- The page-7 record contributes {len(new_warnings)} audit warning item(s): {dict(issue_counts)}. Any `structure_gap` is intentional because the page does not print a concrete product structure; guessing a SMILES would violate the source rules. No missing measurement basis is repaired because no quantitative yield is reported.
- `standard validate`, full `validate --strict`, and the deterministic audit are run after merge. Schema validation is not treated as chemistry verification; `semantic_accuracy_verified` and `chemical_accuracy_verified` remain false.

### Current-session outputs

- Submission extension: `runs/current/submitted.json`
- Merged v3.1 records: `runs/current/standard_records.json`, `dataset.json`, `molecules.json`, `records/`, `classification_index.json`
- Required page coverage: `runs/current/coverage/pages_0007_0007.json`
- Raw visual extraction and source mentions: `runs/current/raw_model/page7_visual_extraction.json`, `page7_source_mentions.json`
- Image crop and lineage receipt: `runs/current/raw_model/page7_notes19_23.png` and `.evidence.json`
- Audit: `runs/current/audit.json`

Unresolved items are limited to source omissions and unavailable, unassociated ESI: the note-19 product identity/stereochemistry and detection basis, plus the concrete identity/conditions behind note 21. No CAS number, product structure or SI association was guessed.
"""
workflow_path.write_text(existing + section)

inventory = {
    "document_id": DOC,
    "pages_reviewed": [7],
    "new_reaction_ids": [RID],
    "source_context_ids": [
        f"{DOC}:page7_references",
        f"{DOC}:page7_footnote19_interpretation",
        f"{DOC}:page7_footnote20_crystallography",
        f"{DOC}:page7_footnote21_onepot",
    ],
    "merged_reaction_count": len(reactions),
    "merged_molecule_count": len(standard["molecules"]),
    "semantic_accuracy_verified": False,
}
rt.write(RUN / "raw_model/page7_final_inventory.json", inventory)

print(json.dumps({
    "coverage_pages": [7],
    "new_reactions": 1,
    "merged_reactions": len(reactions),
    "merged_molecules": len(standard["molecules"]),
    "new_warning_items": len(new_warnings),
    "new_issue_counts": dict(issue_counts),
}, ensure_ascii=False))
