#!/usr/bin/env python3
"""Save final validation/audit receipts and exact page-coverage manifest."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import reaction_tools as rt
import standard_records as st

RUN = ROOT / "runs/current"
data = rt.read(RUN / "dataset.json")
standard = st.read(RUN / "standard_records.json")
base = RUN
strict_ids = {r["record_id"] for r in data["experiments"]}
standard_errors = st.validate(standard, base)
strict_errors = rt.validate(data, base, strict_record_ids=strict_ids)
validation = {
    "standard_bundle": {"valid": not standard_errors, "errors": standard_errors},
    "strict_dataset": {"valid": not strict_errors, "errors": strict_errors, "standard_v3_present": "standard" in data},
    "chemical_accuracy_verified": False,
}
audit = rt.dataset_audit(data)
st.write(RUN / "raw_model/validation.json", validation)
st.write(RUN / "raw_model/audit.json", audit)

ev_page = {e["evidence_id"]: e["page"] for d in data["documents"] for e in d["evidence"]}
records_by_page = {
    page: [r["reaction_id"] for r in standard["reactions"]
           if any(ev_page.get(eid) == page for eid in r["evidence_ids"])]
    for page in range(1, 7)
}

regions = {
    1: [
        {"region": "title, authors, abstract and graphical abstract", "finding": "Current-paper summary of formal haloalkynylation, 27 examples, seven applications, and NMR/DFT studies; graphical summary is not an additional experiment.", "evidence_ids": ["doc_26447b9b5c1dd197:page1:image", "doc_26447b9b5c1dd197:page1:text"]},
        {"region": "two-column introduction", "finding": "Background and cited prior work only; the prose introduces current-work intermediate A but gives no independent procedure on this page.", "evidence_ids": ["doc_26447b9b5c1dd197:page1:image", "doc_26447b9b5c1dd197:page1:text"]},
    ],
    2: [
        {"region": "Scheme 1A-B", "finding": "General/cited iminium chemistry and authors' previous alkynylations; retained as prior-work context, not current-paper experiments.", "evidence_ids": ["doc_26447b9b5c1dd197:page2:image"]},
        {"region": "Scheme 1C", "finding": "One generic current-work formal haloalkynylation scheme; no scheme-level yield assigned.", "evidence_ids": ["doc_26447b9b5c1dd197:page2:image"]},
        {"region": "Scheme 2A", "finding": "Seven first-stage N,O-aminal preparations (2a and 2a-i to 2a-vi), each recorded separately with 0 °C/5 min conditions and its printed yield.", "evidence_ids": ["doc_26447b9b5c1dd197:page2:image"]},
        {"region": "Scheme 2B optimization table and footnotes", "finding": "Every numbered row inventoried; rows with individually assigned substrates were split, two omission controls were separated, and source-aggregated screens 5-10 retain collective ranges rather than invented per-run values.", "evidence_ids": ["doc_26447b9b5c1dd197:page2:image", "doc_26447b9b5c1dd197:page2:text"]},
        {"region": "Scheme 3", "finding": "Nineteen 0.2 mmol alkyne-scope reactions 4aa-4as plus a distinct 5 mmol/96% 4aa scale-up.", "evidence_ids": ["doc_26447b9b5c1dd197:page2:image", "doc_26447b9b5c1dd197:page2:text"]},
        {"region": "prose below Scheme 2", "finding": "Separate 33% sequential one-pot experiment retained; 73% isolated prose outcome for entry 2 retained alongside the table's NMR basis.", "evidence_ids": ["doc_26447b9b5c1dd197:page2:image", "doc_26447b9b5c1dd197:page2:text"]},
    ],
    3: [
        {"region": "Scheme 4A-E and footnotes", "finding": "Four allene-derivative scope reactions, four regioisomer-mixture reactions, four substrate-limit controls, two other-halogen reactions and one elimination product were mapped to records.", "evidence_ids": ["doc_26447b9b5c1dd197:page3:image", "doc_26447b9b5c1dd197:page3:text"]},
        {"region": "Figure 1", "finding": "Computed chlorinated-intermediate energies; non-synthetic source context.", "evidence_ids": ["doc_26447b9b5c1dd197:page3:image", "doc_26447b9b5c1dd197:page3:text"]},
        {"region": "Figure 2", "finding": "One 2a-to-A NMR treatment and one analogous 2a-v complex-spectrum treatment recorded as mechanistic experiments; spectral peaks/HRMS values remain characterization context.", "evidence_ids": ["doc_26447b9b5c1dd197:page3:image", "doc_26447b9b5c1dd197:page3:text"]},
        {"region": "Figure 3", "finding": "Computed transition states only; not an observed reaction experiment.", "evidence_ids": ["doc_26447b9b5c1dd197:page3:image", "doc_26447b9b5c1dd197:page3:text"]},
        {"region": "lower-column prose", "finding": "Confirms the 5 mmol/96% 4aa scale-up and discusses scope; linked to existing records rather than duplicated.", "evidence_ids": ["doc_26447b9b5c1dd197:page3:image", "doc_26447b9b5c1dd197:page3:text"]},
    ],
    4: [
        {"region": "Scheme 5A", "finding": "Two natural-product-bound bromoalkynylations, 4at and 4au.", "evidence_ids": ["doc_26447b9b5c1dd197:page4:image", "doc_26447b9b5c1dd197:page4:text"]},
        {"region": "Scheme 5B and footnotes", "finding": "Seven independent arrows 4aa to 7-13; all printed reagent quantities, solvent ratios, temperatures, times, H2 pressure text and yields were structured.", "evidence_ids": ["doc_26447b9b5c1dd197:page4:image", "doc_26447b9b5c1dd197:page4:text"]},
        {"region": "Scheme 4 continuation prose", "finding": "Links substrate limitations, halogen scope and elimination to page-3 Scheme 4 records; no duplicate experiments created.", "evidence_ids": ["doc_26447b9b5c1dd197:page4:image", "doc_26447b9b5c1dd197:page4:text"]},
        {"region": "computational and NMR discussion", "finding": "Links Figures 1-3 and the two NMR treatments to existing records/context; calculated pathways were not promoted to observed experiments.", "evidence_ids": ["doc_26447b9b5c1dd197:page4:image", "doc_26447b9b5c1dd197:page4:text"]},
    ],
    5: [
        {"region": "left-column Scheme 5 discussion", "finding": "Textual restatement of 4at/4au and products 7-13; linked as supporting evidence to the existing page-4 arrow records.", "evidence_ids": ["doc_26447b9b5c1dd197:page5:image", "doc_26447b9b5c1dd197:page5:text"]},
        {"region": "conclusion and Associated Content", "finding": "No new experiment. SI availability is stated, but no SI was supplied or associated.", "evidence_ids": ["doc_26447b9b5c1dd197:page5:image", "doc_26447b9b5c1dd197:page5:text"]},
        {"region": "CCDC, author information, acknowledgements and beginning of references", "finding": "CCDC 2019531 notice retained as context; no CIF supplied; bibliography is not current-paper experimental evidence.", "evidence_ids": ["doc_26447b9b5c1dd197:page5:image", "doc_26447b9b5c1dd197:page5:text"]},
    ],
    6: [
        {"region": "both columns", "finding": "References only (continuation of citations 2-9); reaction terms occur in cited titles and were not treated as current-paper experiments.", "evidence_ids": ["doc_26447b9b5c1dd197:page6:image", "doc_26447b9b5c1dd197:page6:text"]},
    ],
}
unresolved = {
    1: ["No independent procedure is printed on this page; graphical structures were used only as summary context."],
    2: ["SI is unavailable; individual values behind aggregate optimization rows 5-10 cannot be recovered.", "First-stage yield measurement bases are not stated.", "Many scope structures remain drawing-linked unresolved identities rather than definitive SMILES."],
    3: ["Definitive structures for complex scope products and failed substrates 1i-1k were not safely transcribed.", "Scheme 4 yield bases are generally not reported; ratios are explicitly 1H NMR.", "Computational structures do not establish observed intermediates."],
    4: ["Scheme 5 yield measurement bases and workup/purification are not reported.", "For 4aa to 8, the 2 h applies to the overall 0 °C-to-rt program; stage allocation and ramp time are absent.", "The analogous 2a-v NMR experiment does not separately restate all charges/duration."],
    5: ["SI and CIF are mentioned but not supplied/associated; no details were inferred from them."],
    6: ["Page 6 is bibliography-only; no current-paper chemistry record is expected."],
}
coverage = {
    "document_id": "doc_26447b9b5c1dd197",
    "semantic_accuracy_verified": False,
    "pages": [
        {"page": page, "state": "attempted", "regions": regions[page],
         "record_ids": records_by_page[page], "unresolved": unresolved[page]}
        for page in range(1, 7)
    ],
}
st.write(RUN / "coverage/pages_0001_0006.json", coverage)
print({"validation": validation, "warning_records": len(audit["warnings"]),
       "coverage_record_counts": {p: len(records_by_page[p]) for p in records_by_page}})
