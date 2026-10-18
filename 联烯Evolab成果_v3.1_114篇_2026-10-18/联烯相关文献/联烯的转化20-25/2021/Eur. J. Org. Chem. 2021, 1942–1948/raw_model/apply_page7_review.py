#!/usr/bin/env python3
"""Apply the visually reviewed page-7 characterization/context update.

Page 7 contains no new synthesis experiment.  This script deliberately keeps
the existing molecule and reaction counts stable, adds characterization and SI
context, writes the one-page coverage artifact, and marks only the page-7 task
as attempted.
"""
from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
RUN = ROOT / "runs/current"
DOC = "doc_a10334c5f49b8d92"
PAGE7_IMAGE = f"{DOC}:page7:image"
PAGE7_TEXT = f"{DOC}:page7:text"
REACTION_4A = f"{DOC}:table2_4a"


def write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n")


submitted_path = RUN / "submitted.json"
bundle = json.loads(submitted_path.read_text())
init_template = json.loads((RUN / "raw_model/page7_init_template.json").read_text())

assert bundle["schema_version"] == "3.1"
assert any(m["molecule_id"] == "m_4a" for m in bundle["molecules"])
assert any(r["reaction_id"] == REACTION_4A for r in bundle["reactions"])
original_counts = (len(bundle["molecules"]), len(bundle["reactions"]))

# The merge contract requires this read-only source envelope to exactly match
# the freshly initialized template.  In particular, the prior page-1--6 bundle
# did not yet contain the page-7 rendered-image evidence object.
bundle["documents"] = init_template["documents"]

contexts = {item["context_id"]: item for item in bundle["source_context"]}
contexts[f"{DOC}:page7_4a_characterization"] = {
    "context_id": f"{DOC}:page7_4a_characterization",
    "document_id": DOC,
    "evidence_ids": [PAGE7_IMAGE, PAGE7_TEXT],
    "kind": "characterization",
    "target_molecule_ids": ["m_4a"],
    "related_reaction_ids": [REACTION_4A],
    "raw_data": {
        "1H_NMR_continuation": (
            "4H), 7.02-7.00 (m, 2H), 6.98 (s, 2H), 5.12 (d, J = 1.2 Hz, 2H), "
            "5.03 (s, 1H), 4.98 (s, 1H), 2.51 (s, 3H), 1.34 (s, 18H)"
        ),
        "13C_NMR": (
            "100 MHz, CDCl3: 167.1, 161.0, 152.3, 149.5, 137.2, 136.4, "
            "135.5, 129.0, 128.4, 127.7, 127.4, 127.2, 125.6, 124.5, "
            "124.0, 116.0, 106.3, 65.8, 41.0, 34.2, 30.2, 19.5 ppm"
        ),
        "IR": "3635, 2958, 1712, 1643, 1596, 1493, 1209, 1062, 821, 696 cm-1",
        "melting_point": "148-149 degC",
        "HRMS": "C32H36O4Na+ [M+Na]+: calculated 507.2511; found 507.2506",
    },
    "content": (
        "Page 7 continues the compound 4a characterization begun on page 6. "
        "It contains the end of the 1H NMR list; 13C NMR (100 MHz, CDCl3) "
        "signals at 167.1, 161.0, 152.3, 149.5, 137.2, 136.4, 135.5, "
        "129.0, 128.4, 127.7, 127.4, 127.2, 125.6, 124.5, 124.0, 116.0, "
        "106.3, 65.8, 41.0, 34.2, 30.2 and 19.5 ppm; IR bands at 3635, "
        "2958, 1712, 1643, 1596, 1493, 1209, 1062, 821 and 696 cm-1; m.p. "
        "148-149 degC; and HRMS [M+Na]+ calculated 507.2511, found 507.2506. "
        "These are characterization data, not a separate synthesis experiment."
    ),
}

si_id = f"{DOC}:si_statement"
assert si_id in contexts
contexts[si_id]["evidence_ids"] = list(dict.fromkeys(
    contexts[si_id]["evidence_ids"] + [PAGE7_IMAGE, PAGE7_TEXT]
))
contexts[si_id]["content"] = (
    "The article explicitly states that Supporting Information contains "
    "spectroscopic data, detailed experimental procedures, and the crystal "
    "structure of 4e (CCDC 2044635). No SI file was supplied or explicitly "
    "associated in this workspace, so the SI remains unavailable."
)

xray_id = f"{DOC}:xray_4e"
assert xray_id in contexts
contexts[xray_id]["evidence_ids"] = list(dict.fromkeys(
    contexts[xray_id]["evidence_ids"] + [PAGE7_IMAGE, PAGE7_TEXT]
))

bundle["source_context"] = list(contexts.values())
assert original_counts == (len(bundle["molecules"]), len(bundle["reactions"]))
write_json(submitted_path, bundle)

tasks_path = RUN / "tasks.json"
tasks = json.loads(tasks_path.read_text())
matched = 0
for task in tasks:
    if task["document_id"] == DOC and task["evidence_id"] == PAGE7_TEXT:
        task["visual_evidence_id"] = PAGE7_IMAGE
        task["state"] = "attempted"
        task["reason"] = (
            "Original PDF page 7 image was visually inspected. It contains "
            "4a characterization continuation and non-experimental back matter; "
            "see coverage/pages_0007_0007.json. Attempted does not certify "
            "semantic or chemical accuracy."
        )
        matched += 1
assert matched == 1
write_json(tasks_path, tasks)

coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {
            "page": 7,
            "state": "attempted",
            "regions": [
                {
                    "region": "upper-left continuation of the 4a analytical paragraph",
                    "finding": (
                        "Characterization only: the page continues the 1H NMR "
                        "list for 4a and reports 13C NMR, IR, melting point "
                        "148-149 degC, and HRMS calculated/found values. It is "
                        "linked to the existing 4a preparation and is not a new "
                        "synthetic experiment."
                    ),
                    "evidence_ids": [PAGE7_IMAGE, PAGE7_TEXT],
                },
                {
                    "region": "upper-left Supporting Information paragraph",
                    "finding": (
                        "The article says SI contains spectra, detailed procedures, "
                        "and the 4e crystal structure (CCDC 2044635). No SI was "
                        "supplied or explicitly associated, so no SI facts were imported."
                    ),
                    "evidence_ids": [PAGE7_IMAGE, PAGE7_TEXT],
                },
                {
                    "region": "remaining two columns below/alongside the SI paragraph",
                    "finding": (
                        "Acknowledgements, conflict-of-interest statement, keywords, "
                        "references 1-17, and manuscript dates. These are back matter "
                        "and cited literature, with no current-paper reaction arrow, "
                        "scope/optimization row, mechanistic experiment, or new synthesis."
                    ),
                    "evidence_ids": [PAGE7_IMAGE, PAGE7_TEXT],
                },
            ],
            "record_ids": [REACTION_4A],
            "unresolved": [
                "The 1H NMR list starts on page 6, so page 7 alone is only a continuation; it was not reconstructed as an independent complete spectrum.",
                "The explicitly mentioned Supporting Information was not supplied or associated.",
                "The analytical data do not assign absolute or relative stereochemistry for 4a; its potential C4 stereocenter remains unspecified and its complex structure remains not encoded as a full stereochemical SMILES.",
            ],
        }
    ],
}
coverage_dir = RUN / "coverage"
coverage_dir.mkdir(parents=True, exist_ok=True)
write_json(coverage_dir / "pages_0007_0007.json", coverage)

inventory = {
    "document_id": DOC,
    "scope": {"first_page": 7, "last_page": 7, "claimed_pages": 1},
    "pages": [
        {
            "page": 7,
            "visual_findings": [
                "4a characterization continuation: 1H NMR tail, 13C NMR, IR, melting point, and HRMS",
                "Supporting Information content statement; SI itself not supplied",
                "acknowledgements, conflict of interest, keywords, references, and manuscript dates",
                "no reaction arrows, reaction tables, or new synthesis experiments",
            ],
            "evidence_ids": [PAGE7_IMAGE, PAGE7_TEXT],
            "linked_reaction_ids": [REACTION_4A],
        }
    ],
    "decision": (
        "No new molecule or reaction record was created. Characterization was "
        "stored in source_context and linked to the existing 4a preparation."
    ),
    "semantic_accuracy_verified": False,
}
write_json(RUN / "raw_model/page7_visual_inventory.json", inventory)

print(json.dumps({
    "molecules": len(bundle["molecules"]),
    "reactions": len(bundle["reactions"]),
    "new_reactions": 0,
    "page7_context_added": True,
    "page7_task_attempted": True,
}, ensure_ascii=False))
