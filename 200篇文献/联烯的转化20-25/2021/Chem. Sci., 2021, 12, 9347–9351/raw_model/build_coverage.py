#!/usr/bin/env python3
import json
from pathlib import Path

RUN = Path(__file__).resolve().parents[1]
DOC = "doc_3252f1d5f252557b"
prefix = DOC + ":"


def ids(*suffixes):
    return [prefix + s for s in suffixes]


table1 = [f"table1_entry_{i}" for i in range(1, 12)] + ["table1_entry_9_1mmol_isolated"]
table2 = ["table2_" + x for x in ["3ba", "3ca", "3da", "3ea", "3fa", "3ga", "3ha", "3ia", "3ja", "3ka", "3la", "3ma", "3na", "3oa", "3pa", "3qa", "3ra", "3sa", "3ta", "3ua", "3va"]]
table3 = ["table3_" + x for x in ["3bb", "3pc", "3pd", "3oe", "3bf", "3pg", "3ph", "3pi", "3bj", "3pk", "3pl", "3bm", "3pn", "3bo", "3pp"]]
scheme2 = ["scheme2a_gram_scale_3ga", "scheme2a_3ga_to_4", "scheme2a_3ga_to_5", "scheme2a_3ga_to_6", "scheme2a_3ga_to_7", "scheme2b_eucalyptene_A", "scheme2c_glucose_to_8", "scheme2c_8_to_9", "scheme2c_make_10", "scheme2c_9_plus_10_to_terricollene"]

coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {
            "page": 1,
            "state": "attempted",
            "regions": [
                {"region": "title, abstract, and bibliographic header", "finding": "Actual article title and DOI identified; abstract summarizes the current-work coupling and total syntheses but is not a separate experiment.", "evidence_ids": [f"{DOC}:page1:image", f"{DOC}:page1:text"]},
                {"region": "Scheme 1a", "finding": "Natural-product structures are background examples, not experiments reported in this work; retained in source_context.", "evidence_ids": [f"{DOC}:page1:image"]},
                {"region": "Scheme 1b", "finding": "Generic Pd/Cu aryl-halide ether synthesis is cited prior work; retained as literature context, not counted as a this-paper experiment.", "evidence_ids": [f"{DOC}:page1:image"]},
                {"region": "Scheme 1c", "finding": "One current-work generic phenol/allenylic-carbonate coupling scheme and one summarized unsuccessful aryl-halide/allenol survey. The latter cannot be decomposed into individual trials without the unprovided ESI.", "evidence_ids": [f"{DOC}:page1:image", f"{DOC}:page2:text"]},
            ],
            "record_ids": ids("scheme1c_general", "preliminary_aryl_halide_allenol_survey"),
            "unresolved": ["Natural-product structures and compact generic R groups were not promoted to concrete reaction SMILES.", "The preliminary survey lacks concrete allenol identities, catalyst-by-catalyst conditions, and measurement basis in the supplied pages."],
        },
        {
            "page": 2,
            "state": "attempted",
            "regions": [
                {"region": "Table 1 header, rows 1–11, and footnotes a–b", "finding": "Eleven optimization experiments transcribed separately. Entry 9 footnote b is a distinct 1 mmol repeat with 95% isolated yield, separate from the 0.2 mmol/99% NMR row. NMR yields and starting-material recoveries retain CH3NO2/1H-NMR basis.", "evidence_ids": [f"{DOC}:page2:image", f"{DOC}:page2:text"]},
                {"region": "Table 2 product grid 3ba–3va", "finding": "Twenty-one substrate-scope arrows transcribed individually with isolated yields and footnote-specific catalyst, solvent, temperature, 2a loading, or sealed-tube changes.", "evidence_ids": [f"{DOC}:page2:image", f"{DOC}:page2:text"]},
                {"region": "Table 2 footnotes a–f", "finding": "Default 1.0 mmol phenol, 1.5 mmol 2a, Pd2dba3 2.5 mol%, Xantphos 10 mol%, Et2O 5 mL, rt preserved; b–f exceptions mapped only to marked products.", "evidence_ids": [f"{DOC}:page2:image"]},
                {"region": "discussion prose below Table 1 and beside Table 2", "finding": "Narrative mentions are duplicate evidence for Table 1 entries and Table 2 products; no additional concrete experiment counted.", "evidence_ids": [f"{DOC}:page2:image", f"{DOC}:page2:text"]},
            ],
            "record_ids": ids("preliminary_aryl_halide_allenol_survey", *table1, *table2),
            "unresolved": ["Table 2 does not print reaction times, workups, or purification details.", "Entry 9 footnote b does not restate solvent volume for the 1 mmol repeat; no volume was guessed.", "Most compact scope structures remain label-linked unresolved identities rather than speculative SMILES."],
        },
        {
            "page": 3,
            "state": "attempted",
            "regions": [
                {"region": "Table 3 product grid 3bb–3pp and footnote a", "finding": "Fifteen 2,3-allenylic-carbonate scope arrows transcribed separately with common charges, conditions, and isolated yields. Potential allene axial stereochemistry is unassigned in the source and marked unspecified/not_encoded.", "evidence_ids": [f"{DOC}:page3:image", f"{DOC}:page3:text"]},
                {"region": "Scheme 2a central 3ga and gram-scale arrow", "finding": "One 5 mmol scale-up record with 1.0465 g/93% product, separate from Table 2 3ga.", "evidence_ids": [f"{DOC}:page3:image", f"{DOC}:page3:text"]},
                {"region": "Scheme 2a arrows from 3ga to 4, 5, 6, and (Z)-7", "finding": "Four derivatization experiments recorded individually; yields and Z/E or E/Z order preserved, with unreported analytical bases left missing.", "evidence_ids": [f"{DOC}:page3:image", f"{DOC}:page3:text"]},
                {"region": "Scheme 2b", "finding": "One eucalyptene A synthesis arrow recorded with 1 mmol substrate, 2a 1.5 equiv, Pd(PPh3)4 5 mol%, Et2O, rt, and 90% yield.", "evidence_ids": [f"{DOC}:page3:image", f"{DOC}:page3:text"]},
                {"region": "Scheme 2c", "finding": "Four reaction arrows recorded: D-glucose to 8, 8 to 9, 70 mmol preparation of 10, and 9 plus 10 through methanol/resin treatment to terricollene A. The 58% value remains an overall yield based on glucose, not a last-step yield.", "evidence_ids": [f"{DOC}:page3:image", f"{DOC}:page4:text"]},
            ],
            "record_ids": ids(*table3, *scheme2),
            "unresolved": ["Table 3 and several Scheme 2 arrows omit reaction time, workup, purification, or yield measurement basis.", "Compact structures 4–10 and terricollene A were not converted to complete stereochemical SMILES.", "The source does not report methods for the 94/6 and 98/2 geometric ratios or the sole-anomer statement."],
        },
        {
            "page": 4,
            "state": "attempted",
            "regions": [
                {"region": "Scheme 3", "finding": "Author-proposed Pd(0)/alpha-methylene-pi-allylic-Pd catalytic cycle; recorded as proposed_mechanism, not an observed experiment, and excluded from observed bond-formation classification.", "evidence_ids": [f"{DOC}:page4:image", f"{DOC}:page4:text"]},
                {"region": "left-column Scheme 2c prose continuation", "finding": "Duplicate/expanded evidence for preparation of 10 and the 9 + 10 glycosylation/desilylation sequence. It supplies 4 h for the DCM stage, 4.5 g terricollene A, 58% overall based on glucose, and sole-anomer wording; merged with page-3 records rather than duplicated.", "evidence_ids": [f"{DOC}:page4:image", f"{DOC}:page4:text"]},
                {"region": "conclusions and data availability", "finding": "No new experiment arrow. Data-availability text says ESI exists, but no SI PDF is explicitly associated; retained as a source limitation.", "evidence_ids": [f"{DOC}:page4:image", f"{DOC}:page4:text"]},
            ],
            "record_ids": ids("scheme2c_make_10", "scheme2c_9_plus_10_to_terricollene", "scheme3_proposed_mechanism"),
            "unresolved": ["Quantities for 9, base, DCM, methanol, and resin are absent from the supplied pages.", "Duration of the methanol/resin stage and measurement bases for overall yield/anomeric selectivity are unreported.", "Mechanistic atom mapping and proton transfers are proposed, not established."],
        },
        {
            "page": 5,
            "state": "attempted",
            "regions": [
                {"region": "two-column references 7–20", "finding": "References only; no new table row, reaction arrow, characterization-only spectrum, or current-work experiment appears on this page.", "evidence_ids": [f"{DOC}:page5:image", f"{DOC}:page5:text"]},
            ],
            "record_ids": [],
            "unresolved": ["None beyond normal citation context; cited papers were not treated as experiments in this article."],
        },
    ],
}

known = {r["reaction_id"] for r in json.loads((RUN / "standard_records.json").read_text())["reactions"]}
listed = {rid for page in coverage["pages"] for rid in page["record_ids"]}
assert listed <= known, sorted(listed - known)
assert [p["page"] for p in coverage["pages"]] == [1, 2, 3, 4, 5]
assert all(p["regions"] for p in coverage["pages"])
(RUN / "coverage").mkdir(exist_ok=True)
(RUN / "coverage/pages_0001_0005.json").write_text(json.dumps(coverage, ensure_ascii=False, indent=2) + "\n")
