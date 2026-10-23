#!/usr/bin/env python3
"""Create bounded coverage, saved validation/audit reports, and workflow log."""
import collections
import json
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
RUN = ROOT / "runs/current"
DOC = "doc_9fe9c91a39ef967d"
PY = "/EvoLab/student-workspaces/xllin/OpenEvo-chem-integration/.venv/bin/python"
TOOL = ROOT / "skills/chemical-reaction-extraction/scripts/reaction_tools.py"


def run_json(args):
    p = subprocess.run([PY, str(TOOL), *args], cwd=ROOT, text=True,
                       capture_output=True, check=True)
    return json.loads(p.stdout)


standard_validation = run_json(["standard", "validate", "--bundle", str(RUN / "standard_records.json")])
strict_validation = run_json(["validate", "--dataset", str(RUN / "dataset.json"), "--strict"])
audit = run_json(["audit", "--dataset", str(RUN / "dataset.json")])
(RUN / "standard_validation.json").write_text(json.dumps(standard_validation, ensure_ascii=False, indent=2) + "\n")
(RUN / "strict_validation.json").write_text(json.dumps(strict_validation, ensure_ascii=False, indent=2) + "\n")
(RUN / "audit.json").write_text(json.dumps(audit, ensure_ascii=False, indent=2) + "\n")

bundle = json.loads((RUN / "standard_records.json").read_text())
by_page = {i: [] for i in range(1, 7)}
for reaction in bundle["reactions"]:
    for page in range(1, 7):
        eid = f"{DOC}:page{page}:image"
        if eid in reaction["evidence_ids"]:
            by_page[page].append(reaction["reaction_id"])

regions = {
    1: [
        {"region": "title, abstract, and bibliographic header", "finding": "Present-paper title/DOI and high-level claims; no concrete experimental row on this page.", "evidence_ids": [f"{DOC}:page1:image"]},
        {"region": "Figure 1a-b", "finding": "Prior literature reaction classes and cited examples; not imported as present-paper experiments.", "evidence_ids": [f"{DOC}:page1:image"]},
        {"region": "Figure 1c, This work", "finding": "Two generic [8+2] reaction families (allenamides and ynamides), mapped to generic scheme records only.", "evidence_ids": [f"{DOC}:page1:image"]},
    ],
    2: [
        {"region": "Table 1, entries 1-17 including footnotes", "finding": "Seventeen separate optimization/scale-up experiments; conversion and yield bases retained row by row.", "evidence_ids": [f"{DOC}:page2:image"]},
        {"region": "Results text around Table 1", "finding": "Solvent/catalyst controls, decomposition statements, isolated-yield notes, and scope general conditions.", "evidence_ids": [f"{DOC}:page2:image"]},
    ],
    3: [
        {"region": "Scheme 1", "finding": "Sixteen displayed product rows/routes; 4a is merged with Table 1 entry 16, while the other fifteen are separate scope experiments.", "evidence_ids": [f"{DOC}:page3:image"]},
        {"region": "Scheme 2 and adjacent prose", "finding": "Six crude-mixture time-course experiments plus one isolated-6p thermal rearrangement; 6p isolation footnote is cross-linked to page 6.", "evidence_ids": [f"{DOC}:page3:image"]},
        {"region": "Scheme 3", "finding": "Author-proposed mechanism, stored separately from mechanistic experiments.", "evidence_ids": [f"{DOC}:page3:image"]},
    ],
    4: [
        {"region": "Table 2, entries 1-6 including footnotes", "finding": "Six independent ynamide optimization/control experiments; entry 3 is merged with the 7a scope/procedure evidence.", "evidence_ids": [f"{DOC}:page4:image"]},
        {"region": "Scheme 4", "finding": "Thirteen displayed scope routes total; 7a duplicates Table 2 entry 3, leaving twelve additional records including two distinct 7f routes and two complex-mixture/0% cases.", "evidence_ids": [f"{DOC}:page4:image"]},
    ],
    5: [
        {"region": "Scheme 5", "finding": "Author-proposed ynamide mechanism, not treated as an experimentally proven elementary sequence.", "evidence_ids": [f"{DOC}:page5:image"]},
        {"region": "Typical procedure for 4a", "finding": "Exact charges, 60 °C/1.5 h, PPh3 workup, alumina purification, 79.5 mg and 95% are merged into Table 1 entry 16; spectra remain source context.", "evidence_ids": [f"{DOC}:page5:image"]},
        {"region": "Typical procedure for 7a", "finding": "Exact charges, 100 °C/4.5 h, PPh3 workup, silica purification, 86.4 mg and 90% are merged into Table 2 entry 3; 100 °C conflicts with 110 °C on page 4.", "evidence_ids": [f"{DOC}:page5:image"]},
    ],
    6: [
        {"region": "references 2-23", "finding": "Cited literature only; not imported as present-paper reaction experiments.", "evidence_ids": [f"{DOC}:page6:image"]},
        {"region": "footnotes 19-21", "finding": "Crystallographic note for 4b, instability of 6c, and isolated 6p (52%, 40 °C, 3 h); the latter is linked to its Scheme 2 experiment.", "evidence_ids": [f"{DOC}:page6:image"]},
    ],
}
unresolved = {
    1: ["Figure 1a-b structures are literature context, not present-paper experiments.", "No SI was provided despite the availability statement."],
    2: ["Table 1 entry 10 reports 37% conversion and 40% NMR yield; retained without correction.", "Entry 17 scale-up denominator is ambiguous between footnote and prose.", "The textual '8-phenyl' identity of 1a conflicts with the displayed para-methyl aryl structure/mass evidence."],
    3: ["Scheme 1 does not state the measurement basis for scope yields.", "Complete atom mapping, hydrogen shifts, and bond-order changes are unresolved.", "6c cannot be isolated according to footnote 20."],
    4: ["Scheme 4 does not state the measurement basis for scope yields.", "7d/7e are complex mixtures with no reported component identities or reaction times.", "The 7f-from-3m desilylation reagents/operations are absent from the supplied main text."],
    5: ["7a temperature conflict: 110 °C in Table 2/Scheme 4 versus 100 °C in the Typical Procedure.", "Complex product SMILES and stereochemical identities remain unresolved; spectra do not independently establish the machine structures."],
    6: ["Supporting Information mentioned in footnote 21 was not supplied or associated.", "Page 6 is mostly references; no cited reaction was recast as a present-paper experiment."],
}
coverage = {"document_id": DOC, "semantic_accuracy_verified": False, "pages": []}
for page in range(1, 7):
    coverage["pages"].append({"page": page, "state": "attempted", "regions": regions[page],
                              "record_ids": sorted(set(by_page[page])), "unresolved": unresolved[page]})
(RUN / "coverage").mkdir(parents=True, exist_ok=True)
(RUN / "coverage/pages_0001_0006.json").write_text(json.dumps(coverage, ensure_ascii=False, indent=2) + "\n")

# Mark the bounded page tasks attempted; this is not a semantic-completeness claim.
tasks = json.loads((RUN / "tasks.json").read_text())
for task in tasks:
    if task["document_id"] == DOC and task["task_id"].endswith(":text"):
        task["state"] = "attempted"
        task["reason"] = "Original page image and extracted text inspected; see bounded coverage file."
(RUN / "tasks.json").write_text(json.dumps(tasks, ensure_ascii=False, indent=2) + "\n")

issues = collections.Counter()
for row in audit["warnings"]:
    for warning in row["warnings"]:
        issues[warning.split(":", 1)[0]] += 1

workflow = f"""

---

# Session allene-stereo-v31-200-20261017 — {DOC} — pages 1–6 (2026-10-18)

## Scope and source handling

- Read `AGENTS.md`, the complete extraction skill, and the extraction-standard, literature-schema, stereochemistry, structured-field, quality, and bond-taxonomy references before extraction.
- Visually inspected the supplied original images for PDF pages 1–6. The pre-existing dataset contained zero candidates/experiments; no old experiment record was used as source data.
- Added the printed title: *Gold(I)-Catalyzed [8+2]-Cycloaddition of 8-Aryl-8-azaheptafulvenes with Allenamides and Ynamides: Regioselective Synthesis of Dihydrocycloheptapyrrole Derivatives* (DOI 10.1002/chem.202005348).
- The paper points to Supporting Information, but no SI file was supplied or explicitly associated; SI-only details remain unresolved.

## Page inventory and record decisions

- Page 1: Figure 1a-b is cited prior work; Figure 1c is the present-work overview and maps to two generic-scheme records, not extra concrete experiments.
- Page 2: Table 1 contributes 17 independent rows. Entry 16 is the same 4a experiment repeated in Scheme 1/page 3 and the Typical Procedure/page 5, so those sources are merged. Entry 17 is a distinct scale-up.
- Page 3: Scheme 1 contributes 15 additional scope experiments after deduplicating 4a. Scheme 2 contributes six time-course experiments plus the prose-only isolated-6p rearrangement. Scheme 3 is a proposed mechanism.
- Page 4: Table 2 contributes six independent rows. Entry 3 is the same 7a experiment as Scheme 4/page 4 and the Typical Procedure/page 5. Scheme 4 contributes 12 additional routes after deduplicating 7a, including separate 7f routes from 3f and 3m and the 0%/complex-mixture 7d and 7e cases.
- Page 5: Scheme 5 is a proposed mechanism. The two Typical Procedures enrich the matching 4a/7a records with exact charges, PPh3 workup, purification agents/ratios, isolated mass, and characterization context; NMR/HRMS blocks are not separate syntheses.
- Page 6: references are literature context. Footnote 21's 6p isolation (52%, 40 °C, 3 h) is linked to the matching Scheme 2 record; no reference-paper experiment is imported.

## Conflicts, structures, stereochemistry, and classification

- Retained the 7a temperature conflict exactly: 110 °C in Table 2/Scheme 4 versus 100 °C in the Typical Procedure. The temperature fact is `conflict`, its chosen normalized value is null, both alternatives retain their own page evidence, and the review note remains unresolved. The 1,4-dihydro versus 2-amido-1,6-dihydro naming discrepancy for 7a is also retained without choosing corrected locants.
- Retained Table 1 entry 10's 37% conversion/40% NMR yield without correction, and entry 17's ambiguous 1 mmol scale basis without guessing the denominator.
- Retained the 1a identity wording conflict: prose says 8-phenyl, while the displayed structure and mass evidence indicate a para-methyl aryl group. The drawing-derived candidate SMILES is explicitly unverified.
- Complex scope substrates/products and mechanistic intermediates remain unresolved rather than receiving guessed SMILES. CAS numbers were not supplied and were not guessed. Eleven concrete candidate SMILES were passed through the standard CLI stereo diagnostic; outputs are in `runs/current/raw_model/stereo_checks.json`.
- Every molecule carries explicit chirality, geometry, and representation statuses. Unassigned allene/double-bond or complex product stereochemistry remains unresolved/partial; no absent label is converted to achiral or racemic. No ee/er/dr/E/Z result is reported in the covered pages; Scheme 1 and 2 isomer ratios are stored as product-targeted selectivity outcomes.
- Partial inferred cycloaddition classification contains only the two visually explicit new sigma connections (`C-N`, `C-C`). Complete atom mapping, bond-order relocation, and the proposed 1,5-H shift remain unresolved. Negative/no-reaction rows receive no intended-product labels, and `model_votes` is empty.

## Accounting, validation, Markdown, and audit disposition

- Final v3.1 bundle: {len(bundle['molecules'])} molecules, {len(bundle['reactions'])} reaction/scheme records, and {len(bundle['source_context'])} source-context entries. It contains {sum(r['record_type']=='specific_experiment' for r in bundle['reactions'])} specific experiments, {sum(r['record_type']=='generic_scheme' for r in bundle['reactions'])} generic schemes, and {sum(r['record_type']=='proposed_mechanism' for r in bundle['reactions'])} proposed mechanisms.
- All machine records remain `needs_review` and `partial`; semantic and chemical accuracy flags remain false. `runs/current/coverage/pages_0001_0006.json` contains exactly pages 1–6, each attempted with nonempty regions.
- `standard validate`, `standard merge`, and `validate --strict` passed. Generated Markdown was checked for the 4a/7a exact reagent doses and purification operations, and for the explicit 100/110 °C conflict note.
- Final advisory audit is saved as `runs/current/audit.json`. It has no `reaction_identity`, `outcome_target`, or `missing_steps` warnings. It retains {issues.get('measurement_basis',0)} `measurement_basis` warnings because scope figures do not state a method, and {issues.get('structure_gap',0)} `structure_gap` warning for the isolated 6p→4p rearrangement because no source-safe concrete product/reactant SMILES was asserted. These are documented source limitations, not silently filled assumptions.

## Outputs

- Raw model work: `runs/current/raw_model/build_submission.py`, `visual_inventory.json`, page crops, `source_mentions.json`, `stereo_checks.json`, and the initial audit.
- Formal records: `runs/current/submitted.json`, `standard_records.json`, `dataset.json`, `molecules.json`, and `records/`.
- Rendered/indexed review: `runs/current/规范记录预览.md`, `runs/current/文献实验记录.md`, and `runs/current/classification_index.json`.
- Coverage/QA: `runs/current/coverage/pages_0001_0006.json`, `standard_validation.json`, `strict_validation.json`, and `audit.json`.
"""
workflow_path = ROOT / "runs/workflow.md"
existing = workflow_path.read_text()
session_header = f"# Session allene-stereo-v31-200-20261017 — {DOC} — pages 1–6 (2026-10-18)"
if session_header in existing:
    start = existing.rfind("\n---\n\n", 0, existing.index(session_header))
    existing = existing[:start if start >= 0 else existing.index(session_header)]
workflow_path.write_text(existing.rstrip() + workflow)

print(json.dumps({"coverage_pages": len(coverage["pages"]), "records": len(bundle["reactions"]),
                  "warnings": dict(issues), "standard_valid": standard_validation["valid"],
                  "strict_valid": strict_validation["valid"]}, ensure_ascii=False))
