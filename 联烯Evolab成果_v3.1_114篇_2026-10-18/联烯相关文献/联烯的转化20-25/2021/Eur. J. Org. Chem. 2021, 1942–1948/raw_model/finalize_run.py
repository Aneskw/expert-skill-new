#!/usr/bin/env python3
"""Write deterministic audit, requested page coverage, and task states."""
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "skills/chemical-reaction-extraction/scripts"))
import reaction_tools


DOC = "doc_a10334c5f49b8d92"
base = ROOT / "runs/current"
dataset = json.loads((base / "dataset.json").read_text())
standard = dataset["standard"]
audit = reaction_tools.dataset_audit(dataset)
(base / "audit.json").write_text(json.dumps(audit, ensure_ascii=False, indent=2) + "\n")


def ids(page):
    needles = {f"{DOC}:page{page}:image", f"{DOC}:page{page}:text"}
    return [r["reaction_id"] for r in standard["reactions"] if needles.intersection(r["evidence_ids"])]


coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {
            "page": 1,
            "state": "attempted",
            "regions": [
                {"region": "title, abstract, and author block", "finding": "Article title/DOI and abstract identify the present DMAP-mediated one-pot chromene method; no reaction arrow or table is present here.", "evidence_ids": [f"{DOC}:page1:image", f"{DOC}:page1:text"]},
                {"region": "lower-right Results and Discussion opening", "finding": "Prose reports three initial catalyst trials with 1a/2a and no desired product; the detailed condition drawing is on page 3 and evidence was merged into the three control records.", "evidence_ids": [f"{DOC}:page1:image", f"{DOC}:page1:text", f"{DOC}:page3:image"]},
                {"region": "page footnote", "finding": "Supporting Information is explicitly said to be available online, but no SI file was supplied or associated; status remains mentioned_not_supplied.", "evidence_ids": [f"{DOC}:page1:image", f"{DOC}:page1:text"]},
            ],
            "record_ids": ids(1),
            "unresolved": ["No SI procedure or catalyst-screening details beyond the supplied main article were imported."]
        },
        {
            "page": 2,
            "state": "attempted",
            "regions": [
                {"region": "Scheme 1a-b", "finding": "Two cited prior-work generic transformations (refs. 11 and 12), classified as literature_example rather than current-paper experiments.", "evidence_ids": [f"{DOC}:page2:image"]},
                {"region": "Scheme 1c", "finding": "Present-work generic Rauhut-Currier/aza- or oxa-Michael cascade design, recorded as a generic scheme with R-group templates.", "evidence_ids": [f"{DOC}:page2:image"]},
                {"region": "Figure 1", "finding": "Natural products/pharmaceuticals and a chroman core are contextual structures, not current-paper synthesis experiments; retained in source_context.", "evidence_ids": [f"{DOC}:page2:image"]},
                {"region": "bottom discussion columns", "finding": "Narrative reports 3 at 39%, conversion to 4a, one-pot 4a at 68%, and optimization trends; detailed records cite the corresponding page 3 figure/table.", "evidence_ids": [f"{DOC}:page2:image", f"{DOC}:page2:text", f"{DOC}:page3:image"]},
            ],
            "record_ids": ids(2),
            "unresolved": ["Figure 1 compound structures were not treated as experimental products."]
        },
        {
            "page": 3,
            "state": "attempted",
            "regions": [
                {"region": "Scheme 2 eq. 1", "finding": "Three distinct failed controls: PPh3 (10 mol%, rt-to-60 °C), DABCO (10 mol%, rt), and DMAP (10 mol%, rt), all in toluene with desired 4a not observed.", "evidence_ids": [f"{DOC}:page3:image"]},
                {"region": "Scheme 2 eq. 2", "finding": "Two concrete transformations: 1b+2a to isolated 3 (39%, partial decomposition during purification), then TBAF/THF at rt to 4a (87%).", "evidence_ids": [f"{DOC}:page3:image", f"{DOC}:page2:text"]},
                {"region": "Scheme 2 eq. 3", "finding": "One-pot 1b+2a sequence: DMAP/toluene/rt/12 h, then TBAF in THF/toluene at rt, giving 4a in 68% over two steps.", "evidence_ids": [f"{DOC}:page3:image"]},
                {"region": "Table 1 rows 1-15 and footnotes", "finding": "All 15 optimization rows were recorded independently; entries 1-2 preserve trace product plus decomposition, and entries 3-15 preserve isolated yields. Entry 15 is merged with the identical Table 2 4b presentation.", "evidence_ids": [f"{DOC}:page3:image"]},
            ],
            "record_ids": ids(3),
            "unresolved": ["Scheme 2 eq. 1/2 reaction times or reagent loadings are partly unreported; PPh3 ramp duration is unreported; Scheme 2 eq. 2/3 yield measurement bases remain missing where not explicitly stated."]
        },
        {
            "page": 4,
            "state": "attempted",
            "regions": [
                {"region": "Table 2 common scheme and footnote", "finding": "Common scope conditions: p-QM 1 (0.1 mmol), allenoate 2 (0.12 mmol), DMAP (0.01 mmol), Et2O (1.0 mL), rt 2 h, then TBAF (1 M in THF, 0.1 mmol) 1 h.", "evidence_ids": [f"{DOC}:page4:image"]},
                {"region": "Table 2 products 4a-4o", "finding": "Fifteen scope products were inspected. Fourteen new scope records were created; 4b/75% is the same experiment as Table 1 entry 15 and was not duplicated.", "evidence_ids": [f"{DOC}:page4:image"]},
                {"region": "4e X-ray inset", "finding": "X-ray/CCDC 2044635 is structure characterization, retained in source_context and the 4e note rather than counted as a separate synthesis.", "evidence_ids": [f"{DOC}:page4:image", f"{DOC}:page5:text"]},
            ],
            "record_ids": ids(4),
            "unresolved": ["Absolute/relative configuration at the potential C4 stereocenter is not assigned for 4a-4o; complex structures remain not encoded rather than guessed.", "4a table yield (82%) conflicts with the page 6 experimental paragraph (68%, 33 mg); both are retained with no chosen normalized yield."]
        },
        {
            "page": 5,
            "state": "attempted",
            "regions": [
                {"region": "upper-left scope discussion", "finding": "Narrative confirms yield ranges and identifies allenoates 2b-2d and products 4l-4o; author-proposed explanation for low 4o yield is retained as a note.", "evidence_ids": [f"{DOC}:page5:image", f"{DOC}:page5:text"]},
                {"region": "Scheme 3 enlarged scale", "finding": "1c (2 mmol, 920 mg)+2a (2.4 mmol, 360 mg) under standard one-pot conditions gives 4b (730 mg, 71%).", "evidence_ids": [f"{DOC}:page5:image", f"{DOC}:page5:text"]},
                {"region": "Scheme 3 eqs. 1-3", "finding": "Three derivatizations of 4b: LAH/THF 0 °C-to-rt to 5 (>99%); MeI/NaH/DMF 0 °C-to-rt to 6 (87%); Pd/C/H2/EtOH/rt to 7 (96%).", "evidence_ids": [f"{DOC}:page5:image", f"{DOC}:page5:text"]},
                {"region": "Scheme 3 eq. 4", "finding": "Amino-p-QM 8 control decomposed/messy; intended product was not detected. No product yield or intended bond class was invented.", "evidence_ids": [f"{DOC}:page5:image", f"{DOC}:page5:text"]},
                {"region": "right-column Scheme 4/5 discussion", "finding": "Prose supports the β-ICD reaction to 10 (34%, 70% ee), failed β-ICD trial, and author-proposed mechanism shown on page 6.", "evidence_ids": [f"{DOC}:page5:image", f"{DOC}:page5:text", f"{DOC}:page6:image"]},
            ],
            "record_ids": ids(5),
            "unresolved": ["Derivatization reagent quantities, reaction times, workups, purification and yield measurement bases are unreported.", "The scale-up yield basis is unreported.", "The failed amino-p-QM control lacks a detection method/limit."]
        },
        {
            "page": 6,
            "state": "attempted",
            "regions": [
                {"region": "Scheme 4 eq. 1", "finding": "β-ICD (10 mol%) in 2-Me-THF at rt converts 9+2a to 10 in 34% yield and 70% ee; the starred center lacks an R/S assignment and the ee method is unreported.", "evidence_ids": [f"{DOC}:page6:image", f"{DOC}:page5:text"]},
                {"region": "Scheme 4 eq. 2", "finding": "β-ICD trial with an OTBS-depicted substrate labelled 1a gives intended 3 not detected. This label conflicts with page 3's OH-depicted 1a, so separate molecule IDs are retained.", "evidence_ids": [f"{DOC}:page6:image", f"{DOC}:page3:image"]},
                {"region": "Scheme 5", "finding": "Author-proposed Rauhut-Currier/oxa-Michael cycle with intermediates A-E; recorded as proposed_mechanism, not a mechanistic experiment.", "evidence_ids": [f"{DOC}:page6:image", f"{DOC}:page5:text"]},
                {"region": "Experimental Section General Procedure", "finding": "General scope operation, all quantities, staged times/temperature, water/ether workup, Na2SO4 drying and PE:EA silica chromatography were structured and inherited only by records explicitly following the general procedure.", "evidence_ids": [f"{DOC}:page6:image", f"{DOC}:page6:text"]},
                {"region": "specific 4a paragraph", "finding": "4a is reported as 33 mg, 68% from 1b (0.1 mmol, 43 mg) and 2a (0.12 mmol, 21 mg), conflicting with Table 2's 82%; characterization begins but is not a separate synthesis record.", "evidence_ids": [f"{DOC}:page6:image", f"{DOC}:page6:text", f"{DOC}:page4:image"]},
            ],
            "record_ids": ids(6),
            "unresolved": ["The main article points to unavailable SI for catalyst screening and detailed compound data.", "No absolute configuration or ee measurement method is supplied for 10.", "Page 7 continuation of 4a characterization was outside the requested visual coverage and is not claimed as inspected."]
        },
    ]
}
coverage_dir = base / "coverage"
coverage_dir.mkdir(parents=True, exist_ok=True)
(coverage_dir / "pages_0001_0006.json").write_text(json.dumps(coverage, ensure_ascii=False, indent=2) + "\n")

tasks_path = base / "tasks.json"
tasks = json.loads(tasks_path.read_text())
for task in tasks:
    if task["document_id"] == DOC and any(task["evidence_id"] == f"{DOC}:page{p}:text" for p in range(1, 7)):
        task["state"] = "attempted"
        task["reason"] = "PDF page visually inspected and represented in coverage/pages_0001_0006.json; attempted does not certify semantic accuracy."
tasks_path.write_text(json.dumps(tasks, ensure_ascii=False, indent=2) + "\n")

print(json.dumps({"coverage_pages": len(coverage["pages"]), "audit_warnings": len(audit["warnings"])}))
