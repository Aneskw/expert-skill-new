#!/usr/bin/env python3
"""Write coverage and advisory audit artifacts for the current extraction."""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
RUN = ROOT / "runs/current"
sys.path.insert(0, str(ROOT / "skills/chemical-reaction-extraction/scripts"))
from reaction_tools import dataset_audit

data = json.loads((RUN / "dataset.json").read_text())
doc = "doc_e1c2e1b1afa31a18"
reactions = data["standard"]["reactions"]


def ids_for_page(page):
    token = f"{doc}:page{page}:"
    return [r["reaction_id"] for r in reactions
            if any(e.startswith(token) for e in r.get("evidence_ids", []))]


pages = [
    {"page": 1, "state": "attempted", "regions": [
        {"region": "title, graphical abstract, and abstract", "finding": "Article metadata and a duplicate graphical summary of antipodal JiaPhos-catalyzed annulations; not counted as extra experiments.", "evidence_ids": [f"{doc}:page1:image", f"{doc}:page1:text"]},
        {"region": "right-column catalyst-synthesis prose", "finding": "Prose evidence for Scheme 2 cycloaddition/sulfurization, deprotection, condensation, Grignard addition, Raney-Ni reduction, yields, masses, and scale-up totals.", "evidence_ids": [f"{doc}:page1:image", f"{doc}:page1:text"]},
        {"region": "introduction and prior-work discussion", "finding": "Background/cited reactions only; not current-paper experiments.", "evidence_ids": [f"{doc}:page1:image", f"{doc}:page1:text"]},
    ], "record_ids": ids_for_page(1), "unresolved": ["Graphical abstract structures are generic R-group summaries.", "Detailed catalyst-route operations reside in unavailable SI."]},
    {"page": 2, "state": "attempted", "regions": [
        {"region": "Scheme 1A/B", "finding": "Prior catalyst examples and current design rationale; no new preparative row fabricated.", "evidence_ids": [f"{doc}:page2:image"]},
        {"region": "Scheme 2", "finding": "Six displayed reaction branches/arrows recorded: cycloaddition+sulfurization, deprotection, condensation, Grignard separation, and two Raney-Ni reductions.", "evidence_ids": [f"{doc}:page2:image", f"{doc}:page2:text"]},
        {"region": "Scheme 2 X-ray panels", "finding": "4d/4d′ crystallography retained as characterization/source context, not extra syntheses.", "evidence_ids": [f"{doc}:page2:image", f"{doc}:page2:text"]},
        {"region": "Table 1 header, reaction drawing, rows 1–21, and footnotes a–f", "finding": "One generic optimization scheme plus 21 distinct optimization entries; isolated yield, 8aa/8aa′:9aa ratio by 1H NMR, signed ee by chiral HPLC, additives, substrate charges, catalyst loadings, solvent, time, and temperature retained.", "evidence_ids": [f"{doc}:page2:image", f"{doc}:page2:text"]},
        {"region": "bottom prose", "finding": "Scale-up and prior-art discussion linked to current records where applicable; cited prior reactions remain context.", "evidence_ids": [f"{doc}:page2:image", f"{doc}:page2:text"]},
    ], "record_ids": ids_for_page(2), "unresolved": ["Scheme 2 lacks most times, charges, atmospheres, workups, and yield bases.", "Scheme 2 drawing versus prose conflicts on which 5/5′ label bears the Rp/Sp descriptor string."]},
    {"page": 3, "state": "attempted", "regions": [
        {"region": "Scheme 3 general arrow, 29 product rows, footnotes a–f", "finding": "One generic scheme, 29 standard-scale product rows, and a separate parenthesized 8aa 1 mmol scale-up record; isolated yields, dr by 1H NMR, ee by HPLC, and 80 °C exceptions retained.", "evidence_ids": [f"{doc}:page3:image", f"{doc}:page3:text"]},
        {"region": "Scheme 4 general arrow, 28 product rows, footnotes a–e", "finding": "One generic scheme and 28 5a′-catalyzed substrate-scope experiments; THF exception for 8ae′ and signed negative ee values retained as reported.", "evidence_ids": [f"{doc}:page3:image", f"{doc}:page3:text"]},
        {"region": "lower optimization discussion", "finding": "Duplicate prose for selected Table 1 entries (1, 2, 3, 8, 11, 17) linked to existing optimization rows rather than counted again.", "evidence_ids": [f"{doc}:page3:image", f"{doc}:page3:text"]},
        {"region": "scope discussion", "finding": "Narrative grouping checked against the individual Scheme 3/4 rows.", "evidence_ids": [f"{doc}:page3:image", f"{doc}:page3:text"]},
    ], "record_ids": ids_for_page(3), "unresolved": ["Scheme 3 caption says substrates 1/2 whereas the drawing and prose identify 6/7; conflict retained.", "Complete structures, absolute R/S assignments, HPLC stationary phases, workups, and purification details are unavailable."]},
    {"page": 4, "state": "attempted", "regions": [
        {"region": "continuation of Scheme 4 scope discussion", "finding": "Narrative support for primed scope records, including N-protecting groups, allenoate aryl/ester variation, and high yield/ee statements.", "evidence_ids": [f"{doc}:page4:image", f"{doc}:page4:text"]},
        {"region": "alkyl-substituted methyleneindolinone statement", "finding": "A genuine no-product scope limitation is reported, but no individual substrate labels or row-level charges are given; retained in source_context without fabricating a control product/row.", "evidence_ids": [f"{doc}:page4:image", f"{doc}:page4:text"]},
        {"region": "conclusion", "finding": "Summary only; no additional experiment.", "evidence_ids": [f"{doc}:page4:image", f"{doc}:page4:text"]},
        {"region": "Supporting Information and CCDC notices", "finding": "SI availability and crystallographic accessions are stated, but no SI/CIF was supplied or associated.", "evidence_ids": [f"{doc}:page4:image", f"{doc}:page4:text"]},
        {"region": "author information, acknowledgments, and start of references", "finding": "No current-paper synthetic experiment.", "evidence_ids": [f"{doc}:page4:image", f"{doc}:page4:text"]},
    ], "record_ids": ids_for_page(4), "unresolved": ["No concrete negative-control record can be built from the grouped alkyl-substrate statement.", "SI and CIF details remain unavailable and unassociated."]},
    {"page": 5, "state": "attempted", "regions": [
        {"region": "references 2–15", "finding": "Bibliography/prior work only; no experiment of the current paper.", "evidence_ids": [f"{doc}:page5:image", f"{doc}:page5:text"]},
    ], "record_ids": ids_for_page(5), "unresolved": []},
    {"page": 6, "state": "attempted", "regions": [
        {"region": "continuation of references 15–17", "finding": "Bibliography/prior work only; no experiment of the current paper.", "evidence_ids": [f"{doc}:page6:image", f"{doc}:page6:text"]},
    ], "record_ids": ids_for_page(6), "unresolved": []},
]
coverage = {"document_id": doc, "semantic_accuracy_verified": False, "pages": pages}
(RUN / "coverage").mkdir(exist_ok=True)
(RUN / "coverage/pages_0001_0006.json").write_text(json.dumps(coverage, ensure_ascii=False, indent=2) + "\n")

tasks_path = RUN / "tasks.json"
tasks = json.loads(tasks_path.read_text())
for task in tasks:
    task["state"] = "attempted"
    task["reason"] = "Original page image and text inspected; see coverage/pages_0001_0006.json. Attempted does not certify semantic accuracy."
tasks_path.write_text(json.dumps(tasks, ensure_ascii=False, indent=2) + "\n")

audit = dataset_audit(data)
(RUN / "raw_model/audit.json").write_text(json.dumps(audit, ensure_ascii=False, indent=2) + "\n")
print(json.dumps({"pages": len(pages), "records_by_page": {str(p["page"]): len(p["record_ids"]) for p in pages},
                  "audit_flagged_records": len(audit["warnings"])}, ensure_ascii=False))
