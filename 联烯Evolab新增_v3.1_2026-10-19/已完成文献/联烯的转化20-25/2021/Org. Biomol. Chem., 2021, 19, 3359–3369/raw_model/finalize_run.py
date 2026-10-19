#!/usr/bin/env python3
"""Write bounded page coverage, audit snapshot, task states, and workflow log."""
from pathlib import Path
import collections
import json
import sys

ROOT = Path(__file__).resolve().parents[3]
RUN = ROOT / "runs/current"
sys.path.insert(0, str(ROOT / "skills/chemical-reaction-extraction/scripts"))
import reaction_tools as rt

DOC = "doc_ce3af83f9dcdab7f"
dataset = rt.read(RUN / "dataset.json")
audit = rt.dataset_audit(dataset)
rt.write(RUN / "audit.json", audit)

prefix = DOC + ":"
scheme5_success = [prefix + f"scheme5_6{x}" for x in "abcde"]

coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {
            "page": 7,
            "state": "attempted",
            "regions": [
                {"region": "top continuation of compound 4j characterization", "finding": "remaining 13C NMR/IR/HRMS lines for 4j; linked to the existing Scheme 3 4j experiment, not a new synthesis", "evidence_ids": [f"{DOC}:page7:text", f"{DOC}:page7:image"]},
                {"region": "compound blocks 4k–4r", "finding": "names, reported masses/yields, melting points, eluents, NMR/IR/HRMS; masses/yields augment the corresponding Scheme 3 records and spectral lists remain characterization context", "evidence_ids": [f"{DOC}:page7:text", f"{DOC}:page7:image"]},
            ],
            "record_ids": [prefix + "scheme3_4" + c for c in "jklmnopqr"],
            "unresolved": [
                "Spectral peak lists were inventoried as characterization context but not exhaustively normalized into reaction outcomes.",
                "Sulfoximine sulfur absolute configuration is not assigned; no racemic state is inferred.",
            ],
        },
        {
            "page": 8,
            "state": "attempted",
            "regions": [
                {"region": "compound blocks 4s–4z", "finding": "names, reported masses/yields, melting points, eluents, NMR/IR/HRMS; linked to existing Scheme 3 records, not duplicated as new reactions", "evidence_ids": [f"{DOC}:page8:text", f"{DOC}:page8:image"]},
                {"region": "General procedure for the synthesis of products 6 at page bottom", "finding": "beginning of a cross-page procedure continued on page 9; merged into the five successful Scheme 5 product records", "evidence_ids": [f"{DOC}:page8:text", f"{DOC}:page8:image"]},
            ],
            "record_ids": [prefix + "scheme3_4" + c for c in "stuvwxyz"] + scheme5_success,
            "unresolved": [
                "The general-procedure text calls the starting material N-allenyl-2-iodoaniline 5, while Scheme 5 and its prose call it N-allenylaniline 7; both source statements are retained.",
                "The procedure continues on page 9; neither page alone contains the full operation sequence.",
            ],
        },
        {
            "page": 9,
            "state": "attempted",
            "regions": [
                {"region": "continuation of General procedure for products 6", "finding": "charges, CO purge, room-temperature hold, silica-pad workup and column purification; merged into Scheme 5 records 6a–6e", "evidence_ids": [f"{DOC}:page9:text", f"{DOC}:page9:image"]},
                {"region": "compound blocks 6a–6e", "finding": "reported masses/yields and characterization; product amounts and isolated basis linked to Scheme 5 records", "evidence_ids": [f"{DOC}:page9:text", f"{DOC}:page9:image"]},
                {"region": "General procedure for gram-scale synthesis of 4a", "finding": "one gram-scale current-paper experiment with charges, CO purge, 6 h hold, silica workup, 801 mg and 86% isolated result", "evidence_ids": [f"{DOC}:page9:text", f"{DOC}:page9:image"]},
                {"region": "General procedure for preparation of compound 8 and beginning of compound 8 block", "finding": "4a reduction with ordered 0 °C addition/hold, warming overnight, quench/extraction/drying/purification, followed by 117.5 mg/65% product heading", "evidence_ids": [f"{DOC}:page9:text", f"{DOC}:page9:image"]},
            ],
            "record_ids": scheme5_success + [prefix + "scheme6_gram_scale_4a", prefix + "scheme6_reduction_4a_to_8"],
            "unresolved": [
                "Scheme 5 gives PPh3 10 mol% and 2–6 h; its caption/pages 8–9 procedure give PPh3 20 mol% and the procedure gives 2–3 h.",
                "Scheme 6 gram-scale arrow gives PPh3 10 mol% and toluene 2.0 mL, whereas the page-9 procedure gives 20 mol% and 10 mL.",
                "The reduction procedure says pure product 7, while Scheme 6 and the characterization heading identify product 8.",
            ],
        },
        {
            "page": 10,
            "state": "attempted",
            "regions": [
                {"region": "top continuation of compound 8 characterization", "finding": "remaining NMR/IR/HRMS data for product 8; linked to the existing Scheme 6 reduction record", "evidence_ids": [f"{DOC}:page10:text", f"{DOC}:page10:image"]},
                {"region": "Conflicts of interest, Acknowledgements, Notes and references 1–16", "finding": "non-experimental back matter and cited literature; no current-paper synthesis arrow or table row", "evidence_ids": [f"{DOC}:page10:text", f"{DOC}:page10:image"]},
            ],
            "record_ids": [prefix + "scheme6_reduction_4a_to_8"],
            "unresolved": ["Characterization peak lists are retained in source_context and were not treated as synthetic experiments."],
        },
        {
            "page": 11,
            "state": "attempted",
            "regions": [
                {"region": "Notes and references 16–28", "finding": "reference list only, including CCDC 2048086 and a cited precursor-preparation paper; no current-paper experimental record", "evidence_ids": [f"{DOC}:page11:text", f"{DOC}:page11:image"]},
            ],
            "record_ids": [],
            "unresolved": ["The reference/CCDC/ESI notices do not establish an explicitly supplied SI association in this workspace."],
        },
    ],
}
rt.write(RUN / "coverage/pages_0007_0011.json", coverage)

tasks = rt.read(RUN / "tasks.json")
for task in tasks:
    evidence_id = task.get("evidence_id", "")
    if task.get("document_id") != DOC or ":page" not in evidence_id:
        continue
    page_text = evidence_id.split(":page", 1)[1].split(":", 1)[0]
    if page_text.isdigit() and 7 <= int(page_text) <= 11:
        task["state"] = "attempted"
        task["reason"] = "Page image and extracted text inspected; detailed regions and unresolved items are in coverage/pages_0007_0011.json."
rt.write(RUN / "tasks.json", tasks)

issue_counts = collections.Counter()
for row in audit["warnings"]:
    for warning in row["warnings"]:
        issue_counts[warning.split(":", 1)[0]] += 1

workflow = f"""# Extraction workflow — {DOC}

## Session scope

- Batch: `allene-stereo-v31-200-20261017`.
- Source: *One-pot multi-step cascade protocols toward β-indolyl sulfoximidoyl amides via intermolecular trapping of an α-indolylpalladium complex by CO* (Org. Biomol. Chem. 2021, 19, 3359–3369; DOI 10.1039/D1OB00128K).
- Bounded visual review in this session: original PDF pages 7–11, viewed from `runs/current/assets/{DOC}/page7.jpg` through `page11.jpg` at original resolution. Coverage is `runs/current/coverage/pages_0007_0011.json` and contains exactly those five pages.
- Source mapping/hash came from `runs/current/source_mapping.json`. The fresh submission skeleton came from `standard init` at `runs/current/raw_model/submitted_template_pages_0007_0011.json`.
- No SI file was explicitly supplied or associated. Reference, CCDC and ESI notices on page 11 were not treated as an SI association.
- Existing records were used only for cross-page deduplication. Page-7–11 evidence was merged into the same experiment where it continued a reaction already shown on pages 4–6; no old record was presented as a newly extracted result.

## Page inventory and record accounting

- Page 7: continuation of 4j and complete characterization blocks for 4k–4r. Reported masses/yields were linked to the nine existing Scheme 3 records; spectra remain source context.
- Page 8: characterization blocks for 4s–4z, plus the beginning of the cross-page general procedure for products 6. The latter supports Scheme 5 products 6a–6e.
- Page 9: completion of the products-6 procedure; characterization of 6a–6e; detailed gram-scale 4a procedure; detailed 4a→8 reduction procedure and beginning of product-8 characterization.
- Page 10: completion of product-8 characterization, followed by non-experimental back matter/references.
- Page 11: references only; no current-paper reaction record.
- No new reaction IDs were needed for pages 7–11: all source-derived facts belong to 24 existing records (17 Scheme 3 products 4j–4z, five Scheme 5 products 6a–6e, and two Scheme 6 experiments; Scheme 5/6 overlap in page accounting but reaction IDs remain unique). The cumulative merged dataset remains at {len(dataset['standard']['reactions'])} reactions and {len(dataset['standard']['molecules'])} molecules.

## Source and chemistry decisions

- Characterization headings supplied product masses and yields; independent NMR/IR/HRMS peaks, molecular-ion values and fragments were kept in `source_context`, not turned into synthesis experiments.
- The products-6 general procedure spans pages 8–9. Its Pd(OAc)2 (5 mol%), PPh3 (20 mol%), DBU (0.80 mmol), toluene (2.0 mL), CO purge, 2–3 h hold, silica-pad workup and column purification were associated only with successful products 6a–6e. The Scheme 5 arrow conflicts (PPh3 10 mol%, 2–6 h) remain explicit alternatives.
- The pages 8–9 general procedure calls the starting material “N-allenyl-2-iodoaniline 5”, while Scheme 5 and the surrounding prose identify it as N-allenylaniline 7. The record follows the Scheme 5 reactant relationship but preserves the page-8 wording as an unresolved source-number conflict.
- The gram-scale Scheme 6 arrow/procedure conflict remains unresolved: PPh3 10 versus 20 mol%, and toluene 2.0 versus 10 mL. Both values and evidence are retained.
- Reduction 4a→8 preserves the actual ordered stages: ice-cold dropwise addition, 0 °C for 3 h, warm to room temperature overnight, MeOH/H2O quench, CH2Cl2 extraction, MgSO4 drying and silica purification. The procedure's “product 7” wording conflicts with the Scheme 6/heading label 8 and is retained in notes.
- Every product-specific yield/mass targets its molecule. Isolation basis is attached only where the cross-page general procedure explicitly supports silica purification.
- Sulfoximine sulfur configuration remains unspecified/unencoded. No absent wedge or descriptor is converted into achiral/racemic status. No ee/er/dr/E/Z selectivity is reported on pages 7–11.
- Concrete source-legible candidates were checked with `standard stereo`; diagnostics are in `runs/current/raw_model/stereo_diagnostics.json`. Parsing is treated only as a graph diagnostic, not proof of structure.
- Existing inferred bond labels remain partial and unmapped. No unchanged double bond or bond-order decrease was classified as a newly formed bond, and no model votes were fabricated.

## Validation and audit

- Submission source: fresh `standard init` template, followed by source-derived edits in `runs/current/raw_model/build_submission.py` and import through `standard merge`.
- All machine records remain `needs_review`; semantic and chemical accuracy flags remain false.
- `standard validate`, `validate --strict`, the audit command and the repository tests were run after merge.
- Final audit warning records: {len(audit['warnings'])}; issue counts: {dict(issue_counts)}. Remaining warnings are not auto-guessed; their missing bases or structure uncertainties are preserved for review.

## Outstanding expert review

- Resolve the author numbering conflict (starting material 5 versus 7) and the Scheme 5 time/ligand conflicts if explicitly associated SI becomes available.
- Resolve the gram-scale ligand/solvent conflict and the reduction product-number typo from an authoritative source.
- Confirm drawing-derived connectivity and sulfur stereochemistry before promoting any record beyond `needs_review`.
- Spectral peak lists were inventoried but not exhaustively normalized into molecular characterization fields.
"""
(ROOT / "runs/workflow.md").write_text(workflow)
print(json.dumps({"coverage_pages": len(coverage["pages"]), "audit_warning_records": len(audit["warnings"]), "issue_counts": dict(issue_counts)}, ensure_ascii=False))
