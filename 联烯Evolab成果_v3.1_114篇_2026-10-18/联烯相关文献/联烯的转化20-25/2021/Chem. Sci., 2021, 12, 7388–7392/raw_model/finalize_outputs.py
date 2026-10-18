#!/usr/bin/env python3
"""Capture validation/audit and page coverage after the standard merge."""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
RUN = ROOT / "runs/current"
sys.path.insert(0, str(ROOT / "skills/chemical-reaction-extraction/scripts"))
import reaction_tools as rt

data = json.loads((RUN / "dataset.json").read_text())
strict_ids = {r["record_id"] for r in data["experiments"]}
errors = rt.validate(data, RUN, strict_record_ids=strict_ids)
validation = {
    "valid": not errors,
    "errors": errors,
    "standard_v3_present": "standard" in data,
    "chemical_accuracy_verified": False,
}
audit = rt.dataset_audit(data)
(RUN / "validate_strict.json").write_text(json.dumps(validation, ensure_ascii=False, indent=2) + "\n")
(RUN / "audit.json").write_text(json.dumps(audit, ensure_ascii=False, indent=2) + "\n")

doc = "doc_8a0a9af017c80023"
reactions = data["standard"]["reactions"]
by_page = {}
for page in range(1, 6):
    ids = {f"{doc}:page{page}:text", f"{doc}:page{page}:image"}
    by_page[page] = [r["reaction_id"] for r in reactions if ids.intersection(r["evidence_ids"])]

coverage = {
    "document_id": doc,
    "semantic_accuracy_verified": False,
    "pages": [
        {
            "page": 1,
            "state": "attempted",
            "regions": [
                {"region": "title, abstract, DOI and article metadata", "finding": "Article identity confirmed; abstract describes geometry-, enantio- and diastereoselective hydroalkylation and acidic hydrolysis.", "evidence_ids": [f"{doc}:page1:image", f"{doc}:page1:text"]},
                {"region": "Scheme 1(a)", "finding": "Clayden et al. intramolecular N→C migration is cited prior work; stored as source_context, not this-paper experiment.", "evidence_ids": [f"{doc}:page1:image", f"{doc}:page1:text"]},
                {"region": "Scheme 1(b)", "finding": "This-work generic intermolecular hydroalkylation of alkynes or allenes; one generic-scheme record with R/Ar templates.", "evidence_ids": [f"{doc}:page1:image", f"{doc}:page1:text"]},
            ],
            "record_ids": by_page[1],
            "unresolved": ["Generic R and Ar groups do not define concrete molecules.", "ESI availability is stated but no SI file was explicitly supplied or associated."],
        },
        {
            "page": 2,
            "state": "attempted",
            "regions": [
                {"region": "Table 1, entries 1–8 including footnotes", "finding": "Eight distinct optimization experiments transcribed with base loading, temperature, solvent exception, 48 h exception, isolated/NMR yields, dr and er.", "evidence_ids": [f"{doc}:page2:image", f"{doc}:page2:text"]},
                {"region": "Scheme 2, products 3a–3ai", "finding": "Thirty-five substrate-scope experiments transcribed individually with isolated yield, dr and reported er/not-determined state.", "evidence_ids": [f"{doc}:page2:image", f"{doc}:page2:text"]},
                {"region": "scope discussion below Scheme 2", "finding": "Halogen, naphthyl, thiophenyl, benzodioxole, polyfluoro, enyne and bis-alkyne observations linked as interpretation; no duplicate experiments created.", "evidence_ids": [f"{doc}:page2:image", f"{doc}:page2:text"]},
            ],
            "record_ids": by_page[2],
            "unresolved": ["Compact stereochemical drawings were not expanded into unverified complete SMILES.", "Scheme 2 does not state an er measurement method; footnote-a entries explicitly have er not determined."],
        },
        {
            "page": 3,
            "state": "attempted",
            "regions": [
                {"region": "Scheme 3(a)", "finding": "Seven arrows recorded: esterification, amidation, diastereomer-forming cyclization with two products, two methylations and two stereodivergent hydroalkylations.", "evidence_ids": [f"{doc}:page3:image", f"{doc}:page3:text"]},
                {"region": "Scheme 3(b)", "finding": "Large-scale 3a preparation and microwave acidic hydrolysis to amino-acid 4a recorded as two experiments; 74% is explicitly an isolated yield in prose.", "evidence_ids": [f"{doc}:page3:image", f"{doc}:page3:text"]},
                {"region": "Scheme 3(c)", "finding": "Hydrolysis products 4a–4e inventoried; 4a merged with the immediately preceding hydrolysis arrow, and 4b–4e retained as separate derivatizations.", "evidence_ids": [f"{doc}:page3:image", f"{doc}:page3:text"]},
                {"region": "Scheme 4 and X-ray inset", "finding": "Seven allene-scope experiments recorded; X-ray of 3ac stored as characterization context rather than a synthesis.", "evidence_ids": [f"{doc}:page3:image", f"{doc}:page3:text"]},
            ],
            "record_ids": by_page[3],
            "unresolved": ["Scheme 3(c) does not report full hydrolysis conditions or yield measurement basis for 4b–4e.", "The source reuses/overlaps label 4a for the amino-acid product and the allene caption; these are separate molecule IDs.", "Several Scheme 3 yields lack an explicit isolated/NMR basis, and Scheme 4 footnote-a er values were not determined."],
        },
        {
            "page": 4,
            "state": "attempted",
            "regions": [
                {"region": "Figure 1(a)", "finding": "Author-proposed pathway stored as proposed_mechanism; the LiHMDS/CD3OD deuteration control is a separate mechanistic experiment.", "evidence_ids": [f"{doc}:page4:image", f"{doc}:page4:text"]},
                {"region": "Figure 1(b) and accompanying prose", "finding": "In situ IR experiment and n-BuLi control recorded separately; 1926, 1666 and 1712 cm−1 retained as source strings because cm−1 is not a supported structured unit.", "evidence_ids": [f"{doc}:page4:image", f"{doc}:page4:text"]},
                {"region": "Figure 1(c), entries 1–4", "finding": "Four deuterium-labeling experiments transcribed row by row with x/y site-specific D incorporation; entry 4 includes the D2O quench.", "evidence_ids": [f"{doc}:page4:image", f"{doc}:page4:text"]},
                {"region": "methyl iodide quench prose", "finding": "No-alkylated-product control recorded without inventing a product identity or bond classification.", "evidence_ids": [f"{doc}:page4:image", f"{doc}:page4:text"]},
                {"region": "Figure 1(d)", "finding": "DFT transition states and relative energies stored as computational source_context, not synthetic experiments.", "evidence_ids": [f"{doc}:page4:image", f"{doc}:page4:text"]},
            ],
            "record_ids": by_page[4],
            "unresolved": ["The prose does not unambiguously assign the 65% D quench result to the n-BuLi versus LiHMDS IR mixture; it is recorded once with the Figure 1(a) LiHMDS sequence and flagged.", "Main text does not state the measurement method for Figure 1(c) x/y percentages.", "Transient intermediates and DFT structures are not isolated molecular identities."],
        },
        {
            "page": 5,
            "state": "attempted",
            "regions": [
                {"region": "conclusion, author contributions, conflicts and acknowledgements", "finding": "No additional reaction arrow; conclusion repeats the method scope without creating duplicate records.", "evidence_ids": [f"{doc}:page5:image", f"{doc}:page5:text"]},
                {"region": "footnote 13", "finding": "HMPA and [12]crown-4 lithium-chelation controls recorded separately with reported ee values of 73% and 51%.", "evidence_ids": [f"{doc}:page5:image", f"{doc}:page5:text"]},
                {"region": "references 1–22", "finding": "Bibliographic citations are context only; no cited-paper synthesis was imported as this work.", "evidence_ids": [f"{doc}:page5:image", f"{doc}:page5:text"]},
            ],
            "record_ids": by_page[5],
            "unresolved": ["Footnote 13 does not report additive loading, reaction yield, time/temperature or ee measurement method.", "Reference 20 points to ESI details, but SI was not supplied."],
        },
    ],
}
(RUN / "coverage/pages_0001_0005.json").write_text(json.dumps(coverage, ensure_ascii=False, indent=2) + "\n")

tasks_path = RUN / "tasks.json"
tasks = json.loads(tasks_path.read_text())
for task in tasks:
    if task.get("document_id") == doc and task.get("evidence_id") in {f"{doc}:page{p}:text" for p in range(1, 6)}:
        task["state"] = "attempted"
        task["reason"] = "Original page image and extracted text inspected; page disposition is recorded in coverage/pages_0001_0005.json."
tasks_path.write_text(json.dumps(tasks, ensure_ascii=False, indent=2) + "\n")

print(json.dumps({"validation": validation, "flagged_records": len(audit["warnings"]), "coverage_pages": len(coverage["pages"])}, ensure_ascii=False))
