#!/usr/bin/env python3
"""Write coverage, audit snapshot, task states, and workflow decisions."""
import json
from collections import Counter
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[3]
RUN = ROOT / "runs/current"
sys.path.insert(0, str(ROOT / "skills/chemical-reaction-extraction/scripts"))
import reaction_tools as rt
import standard_records as sr

DOC = "doc_68d6021490685b69"
data = json.loads((RUN / "dataset.json").read_text())
audit = rt.dataset_audit(data)
rt.write(RUN / "audit.json", audit)
rt.write(RUN / "raw_model/source_mentions_1a_2a_3aa_4aa_5aa.json",
         sr.source_mentions(data, ["1a", "2a", "3aa", "4aa", "5aa"], context=180))

reactions = data["standard"]["reactions"]
page_records = {}
for page in range(1, 7):
    eid = f"{DOC}:page{page}:text"
    page_records[page] = [r["reaction_id"] for r in reactions if eid in r["evidence_ids"]]

coverage = {
 "document_id": DOC,
 "semantic_accuracy_verified": False,
 "pages": [
  {"page": 1, "state": "attempted", "regions": [
    {"region": "title, abstract and bibliographic block", "finding": "Article title/DOI and abstract verified visually; abstract is context rather than a concrete experiment.", "evidence_ids": [f"{DOC}:page1:image", f"{DOC}:page1:text"]},
    {"region": "Scheme 1a", "finding": "Two generic prior-literature [2+2] routes; retained as literature_example records, not current-work experiments.", "evidence_ids": [f"{DOC}:page1:image"]},
    {"region": "Scheme 1b", "finding": "Two generic current-work transformations: asymmetric [2+2] product family 3 and sequential product family 4.", "evidence_ids": [f"{DOC}:page1:image"]}],
   "record_ids": page_records[1], "unresolved": ["Generic R groups and individual conditions are not concrete on this page."]},
  {"page": 2, "state": "attempted", "regions": [
    {"region": "Table 1 rows 1–11 and footnotes a–e", "finding": "Eleven independent optimization experiments transcribed, including combined yield, 3aa:4aa ratio, signed −11% ee, entry-3 1:1.5 ligand/Ni ratio, and entry-9 condition changes.", "evidence_ids": [f"{DOC}:page2:image", f"{DOC}:page2:text"]},
    {"region": "Table 2 rows 1–6 and footnotes a–e", "finding": "Six independent Lewis-acid/sequential-transformation experiments transcribed; entry 6 is isolated-3aa two-step control rather than a one-pot 1a+2a experiment.", "evidence_ids": [f"{DOC}:page2:image", f"{DOC}:page2:text"]},
    {"region": "results prose", "finding": "Prose repetitions of Table 1 outcomes were linked to the table records and not double counted.", "evidence_ids": [f"{DOC}:page2:text"]}],
   "record_ids": page_records[2], "unresolved": ["Optimization workups/purifications and additional ESI trials are unavailable because SI was not provided."]},
  {"page": 3, "state": "attempted", "regions": [
    {"region": "Table 3 products 3aa–3ag", "finding": "Twenty-three distinct [2+2] scope arrows transcribed with isolated yields and HPLC ee values.", "evidence_ids": [f"{DOC}:page3:image", f"{DOC}:page3:text"]},
    {"region": "Table 3 products 4aa–4ag", "finding": "Twenty-three distinct one-pot sequential-product arrows transcribed, including Cu-loading footnotes and 4fa:3fa / 4la:3la ratios.", "evidence_ids": [f"{DOC}:page3:image", f"{DOC}:page3:text"]},
    {"region": "Table 3 cond. 1 and footnote a", "finding": "Source conflict retained: scheme prints 35 °C while footnote inherits Table 1 entry 8 at 30 °C; no value chosen.", "evidence_ids": [f"{DOC}:page3:image", f"{DOC}:page3:text"]},
    {"region": "left-column control prose", "finding": "Forward 3aa→4aa (37%, tiny unidentified by-product) and reverse no-transformation controls recorded separately.", "evidence_ids": [f"{DOC}:page3:text"]},
    {"region": "bottom-right 5aa prose", "finding": "Repeated evidence for Scheme 2a transformations and 5aa dr/ee linked to the page-4 records, not duplicated.", "evidence_ids": [f"{DOC}:page3:text"]}],
   "record_ids": page_records[3], "unresolved": ["Compact structures are retained by source labels when concrete stereochemical SMILES cannot be defended.", "Individual procedures and measurement methods outside the table footnote require unavailable SI."]},
  {"page": 4, "state": "attempted", "regions": [
    {"region": "Scheme 2a", "finding": "TFA/H2O 4aa→5aa and Cu(OTf)2/air 3aa→5aa arrows inspected; the latter is the same 3aa→5aa experiment repeated in Scheme 2b and counted once.", "evidence_ids": [f"{DOC}:page4:image", f"{DOC}:page4:text"]},
    {"region": "Scheme 2b", "finding": "Five substrate-to-gamma-lactone arrows transcribed with yield, dr and paired ee values.", "evidence_ids": [f"{DOC}:page4:image"]},
    {"region": "Scheme 2c", "finding": "Three alkoxyallene scope arrows to 7aa–7ac transcribed.", "evidence_ids": [f"{DOC}:page4:image", f"{DOC}:page4:text"]},
    {"region": "Scheme 3a–e", "finding": "Two scale-ups, one gamma-lactone scale conversion, and two hydrogenations transcribed; the 58% value for 4qa remains an overall two-step yield.", "evidence_ids": [f"{DOC}:page4:image", f"{DOC}:page4:text"]},
    {"region": "mechanistic discussion", "finding": "Thermal 3aa→E-4aa control at 80 °C for 24 h recorded qualitatively; DFT discussion retained as context.", "evidence_ids": [f"{DOC}:page4:text"]}],
   "record_ids": page_records[4], "unresolved": ["Many Scheme 2/3 yield and selectivity measurement methods are not printed.", "Reaction times and workups absent from the scheme remain unreported unless explicitly inherited from a named table procedure."]},
  {"page": 5, "state": "attempted", "regions": [
    {"region": "Figure 1", "finding": "Computational energy profile for E-/Z-4aa formation retained as non-synthetic source_context.", "evidence_ids": [f"{DOC}:page5:image", f"{DOC}:page5:text"]},
    {"region": "Figure 2a", "finding": "Proposed stereochemical working model and transition state; retained within one proposed_mechanism record, not an observed experiment.", "evidence_ids": [f"{DOC}:page5:image"]},
    {"region": "Figure 2b and accompanying prose", "finding": "Proposed cycles A/B for 3aa/4aa→5aa; arrow-pushing intermediates are not registered as isolated compounds.", "evidence_ids": [f"{DOC}:page5:image", f"{DOC}:page5:text"]}],
   "record_ids": page_records[5], "unresolved": ["Detailed computations and supporting experiments cited to ESI are not available."]},
  {"page": 6, "state": "attempted", "regions": [
    {"region": "conclusion and data availability", "finding": "No new reaction arrow or experiment; ESI availability is stated but no SI file is associated in this workspace.", "evidence_ids": [f"{DOC}:page6:image", f"{DOC}:page6:text"]},
    {"region": "author contributions, conflicts, acknowledgements and references", "finding": "No current-work synthetic experiment or independent characterization table to extract.", "evidence_ids": [f"{DOC}:page6:image"]}],
   "record_ids": page_records[6], "unresolved": ["PDF page 7 and the ESI were outside the requested visual coverage and are not claimed as reviewed."]}
 ]}
rt.write(RUN / "coverage/pages_0001_0006.json", coverage)

tasks_path = RUN / "tasks.json"
tasks = json.loads(tasks_path.read_text())
for task in tasks:
    page = int(task["evidence_id"].split(":page", 1)[1].split(":", 1)[0])
    if 1 <= page <= 6:
        task["state"] = "attempted"
        task["reason"] = "Original page image inspected; experiment/context inventory is recorded in coverage/pages_0001_0006.json."
rt.write(tasks_path, tasks)

issue_counts = Counter()
for item in audit["warnings"]:
    for warning in item["warnings"]:
        issue_counts[warning.split(":", 1)[0]] += 1

workflow = f"""# Extraction workflow: Chem. Sci. 2021, 12, 9991–9997, PDF pages 1–6

- Batch: `allene-stereo-v31-200-20261017`
- Document: `{DOC}`
- Title: *{data['documents'][0]['title']}*
- DOI: `10.1039/d1sc02681j`
- Requested visual scope: original PDF pages 1–6 exactly. PDF page 7 was available as text context but is not claimed as visually covered.
- Review state: all {len(reactions)} v3.1 reactions remain `needs_review`; `semantic_accuracy_verified` and `chemical_accuracy_verified` remain false.

## Source handling and workflow

1. Read `AGENTS.md`, the complete extraction skill, extraction standard, stereochemistry rules, structured-fields rules, quality rules, v3.1 schema, and bond taxonomy before extraction.
2. Confirmed that the prepared dataset initially contained zero experiments. The pre-existing `runs/workflow.md` described a different paper and was treated as stale output, never as source data.
3. Visually inspected each supplied original page image for pages 1–6. Saved lineage-tracked crops of page 2 Tables 1–2, page 3 Table 3, page 4 Schemes 2–3, and page 5 Figures 1–2 under `runs/current/raw_model/` with `.evidence.json` receipts.
4. Initialized `runs/current/submitted.json` with `standard init`, authored only v3.1 molecule/reaction records, and imported with `standard merge`. The legacy experiment view was generated automatically.
5. `runs/current/source_mapping.json` identifies only the main article. The article says ESI exists, but no SI file was supplied or explicitly associated; no SI procedure, spectrum, yield basis, or condition was inferred.
6. Added the source title to document metadata. Source PDF/evidence, tools, schemas, tests, and scientific rules were not modified.

## Page and experiment inventory

- Page 1: title/abstract; two Scheme 1a prior-literature generic routes; two Scheme 1b current-work generic transformations.
- Page 2: Table 1 rows 1–11 and Table 2 rows 1–6, for 17 independent optimization records.
- Page 3: Table 3 has 23 distinct product-3 arrows and 23 distinct product-4 arrows. Two real forward/reverse controls are separate. Repeated 5aa prose is linked to page-4 records rather than double counted.
- Page 4: Scheme 2 has nine unique experiments (the 3aa→5aa arrow is repeated in panels a/b), Scheme 3 has five experiments, and the thermal 3aa→E-4aa control is one mechanistic-study experiment.
- Page 5: Figure 1 is computational context; Figure 2 is one proposed-mechanism record, not a performed experiment.
- Page 6: conclusion/data availability/authorship/references; no new current-work reaction.

Final bundle: {len(data['standard']['molecules'])} molecule/agent records and {len(reactions)} reaction records. Compatibility counts are {json.dumps(audit['record_kind_counts'], ensure_ascii=False)}.

## Identity and stereochemistry decisions

- Every molecule has explicit `chirality_status`, `geometry_status`, and `representation_status`. Missing source assignments remain `unspecified`/`unresolved`; none were converted to racemic or achiral merely by omission.
- Candidate structures were entered only for fully legible model substrates 1a and 2a. `standard stereo` diagnostics are saved in `raw_model/stereo_diagnostics.json`. 1a retains the source E geometry; 2a remains `partial` because toolkit cumulene diagnostics do not certify axial representation.
- Product 3aa retains the source `(1R,2R,4S)` X-ray assignment. Other product-3 drawings retain a single drawn diastereomer and the authors' CD-comparison statement without inventing explicit R/S labels.
- Product 4ka retains `(R,E)`; 4aa retains the explicit mechanistic designation `E-4aa`. Other product-4 absolute configurations remain unnamed. Gamma-lactones 5 retain experiment-specific dr and separate ee values for each diastereomer.
- Products 7–9 retain ee/dr outcomes but no unprinted absolute configurations. Compact source structures without defensible stereochemical SMILES remain label-linked and `not_encoded`.
- CAS numbers were not supplied and remain unreported. No external database enrichment was performed.

## Conditions and outcomes

- Table 1 retains 1a 0.05 mmol, 2a 1.0 or 1.2 equiv, catalyst/ligand loadings, CH2Cl2 concentration, 30/35 °C, 16 h, N2, combined isolated yield, 3aa:4aa ratio, and separate HPLC ee values. Entry 11's signed `−11% ee` is preserved as a signed textual selectivity rather than rejected as a negative yield.
- Table 2 retains two ordered stages. Entry 5 includes 4 Å MS (10 mg), CHCl3 (0.05 M), 30 °C, 3 h and air; entry 6 correctly starts from isolated 3aa and uses 4 h under N2.
- Table 3 product-3 records preserve the unresolved source conflict: `cond. 1` prints 35 °C, while footnote a inherits Table 1 entry 8 at 30 °C. Both alternatives are evidence-linked; no chosen value is asserted.
- Table 3 product-4 records retain footnote-specific Cu loadings (10, 7.5 or 12.5 mol%) and the reported 4fa:3fa / 4la:3la ratios.
- Scheme 2 retains TFA/H2O, Cu(OTf)2 loadings, air, temperatures, yields, dr and ee exactly where printed. Missing times and measurement methods remain null.
- Scheme 3 retains all printed quantities: 2.0/2.4 mmol and 1.12 g; 1.60/1.96 mmol, Cu(OTf)2 7.5 mol%, 4 Å MS 400 mg, CHCl3 0.07 M and 0.62 g; 3ab 2.5 mmol and 0.72 g; Pd/C 10 or 100 wt%, H2 5 MPa, EtOAc 0.05 M. The 58% 4qa value remains an overall two-step yield.
- Every quantitative yield/selectivity/product amount targets molecule IDs. Measurement basis is a separate fact; when absent in the supplied pages it remains missing rather than being called isolated or HPLC by convention.

## Bond-formation classification

- Product-3 [2+2] records are partially classified as two new C–C connections. Sequential/formal [4+2] products retain partial C–C/C–O formation labels.
- 3→4 rearrangements retain partial C–O formation; gamma-lactone formation retains partial C–O and N–H formation when supported by deprotection/lactonization drawings.
- Hydrogenations are labelled C–H only. Reduction of C=C is retained conceptually as bond-order decrease and is not mislabelled as new C–C sigma formation.
- Negative control and proposed mechanism classifications are `not_applicable`; no expected product or model-vote set was fabricated. All classification statements are inferred/structure-derived and no complete atom mapping is claimed.

## Audit dispositions

- `standard validate` and `validate --strict` pass with zero errors and v3 present.
- Audit flags {len(audit['warnings'])} records: {issue_counts.get('structure_gap', 0)} `structure_gap` and {issue_counts.get('measurement_basis', 0)} `measurement_basis` warnings; no `reaction_identity`, `outcome_target`, `product_identity`, or `missing_steps` warning remains.
- `structure_gap` is intentionally retained for complex label-linked structures that cannot be safely encoded from the compact pages. Guessing SMILES solely to silence the warning would violate the source rules.
- `measurement_basis` is intentionally retained where the main-text prose/scheme says “yield” but does not say isolated, NMR, or another method. No method was borrowed from a different table or experiment.
- The main document remains in `documents_without_linked_si`, matching the explicit unestablished SI association.

## Outputs

- Submission/merged data: `runs/current/submitted.json`, `dataset.json`, `standard_records.json`, `molecules.json`, `records/`, and `classification_index.json`
- Human-readable review: `runs/current/规范记录预览.md` and compatibility `runs/current/文献实验记录.md`
- Required coverage: `runs/current/coverage/pages_0001_0006.json`
- Raw extraction artifacts: `runs/current/raw_model/page1_6_visual_inventory.json`, crops and crop-lineage receipts, `stereo_diagnostics.json`, and the session builder scripts
- Deterministic audit: `runs/current/audit.json`

All chemistry-tool commands used `/EvoLab/student-workspaces/xllin/OpenEvo-chem-integration/.venv/bin/python`. All records remain machine-generated `needs_review`; passing validation does not establish chemistry correctness, completeness, or reproducibility.
"""
(ROOT / "runs/workflow.md").write_text(workflow)
print(json.dumps({"coverage_pages": 6, "warnings": len(audit["warnings"]),
                  "issue_counts": issue_counts, "reactions": len(reactions)},
                 ensure_ascii=False, default=dict))
