#!/usr/bin/env python3
"""Build the bounded pages 7-9 v3.1 update for doc_403ed550996f2a54."""

from __future__ import annotations

import copy
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))

import standard_records as st


DOC = "doc_403ed550996f2a54"
RUN = ROOT / "runs/current"
DATASET = st.read(RUN / "dataset.json")
CURRENT = st.read(RUN / "standard_records.json")
EV = {page: f"{DOC}:page{page}:text" for page in range(1, 10)}
IMG = {page: f"{DOC}:page{page}:image" for page in range(1, 10)}


def unique(items):
    return list(dict.fromkeys(items))


def add_evidence(fact, *evidence_ids):
    fact["evidence_ids"] = unique(fact["evidence_ids"] + list(evidence_ids))


def find_molecule(molecule_id):
    return copy.deepcopy(next(m for m in CURRENT["molecules"] if m["molecule_id"] == molecule_id))


def find_reaction(reaction_id):
    return copy.deepcopy(next(r for r in CURRENT["reactions"] if r["reaction_id"] == reaction_id))


def f(raw=None, normalized=None, evidence=(), method="original_extraction", derivation=None):
    return st.fact(raw, normalized, evidence, method=method, derivation=derivation)


# Page 7 completes the analytical paragraph for 5d. It adds characterization
# evidence to the already-existing experiment but no new synthetic operation or
# outcome.
r5d = find_reaction(f"{DOC}:diyne_4")
r5d["evidence_ids"] = unique(r5d["evidence_ids"] + [EV[7]])
location_5d = copy.deepcopy(r5d["location"]["normalized_value"])
location_5d["pages"] = sorted(set(location_5d.get("pages", []) + [7]))
r5d["location"] = f(
    location_5d,
    location_5d,
    unique(r5d["location"]["evidence_ids"] + [EV[7]]),
)
if not any(n["normalized_value"] == "5d characterization continuation; no new experiment" for n in r5d["notes"]):
    r5d["notes"].append(
        f(
            "The 13C NMR list and HRMS characterization for 5d continue at the top of PDF page 7; these analytical data are not a separate synthesis.",
            "5d characterization continuation; no new experiment",
            [EV[7]],
        )
    )


# Page 7 supplies the full source name and characterization paragraph for 5e.
m5e = find_molecule("m_5e")
full_name_5e = (
    "Cholesteryl (1S*,6R*)-7-cyclopropyl-8-(cyclopropylethynyl)-"
    "9-azabicyclo[4.2.1]nona-2,4,7-triene-9-carboxylate equivalent with "
    "cholesteryl (1R*,6S*)-7-cyclopropyl-8-(cyclopropylethynyl)-"
    "9-azabicyclo[4.2.1]nona-2,4,7-triene-9-carboxylate"
)
m5e["names"]["original"] = f(full_name_5e, full_name_5e, [EV[7]])
m5e["names"]["english"] = f(full_name_5e, full_name_5e, [EV[7]])
occ = m5e["occurrences"][0]
occ["evidence_ids"] = unique(occ["evidence_ids"] + [EV[7]])
add_evidence(occ["compound_label"], EV[7])

# Retain the source's relative bridgehead descriptors and the cholesterol
# stereochemistry without claiming a complete isomeric graph.
atom_centers = {
    "bridgeheads": "relative descriptors (1S*,6R*) equivalent with (1R*,6S*) are reported"
}
m5e["stereochemistry"]["atom_centers"] = f(atom_centers, atom_centers, [EV[7]])
# Page 3, not page 7, explicitly supplies the 1:1 ratio. Page 7 only states
# that two rotamers exist, so it supports the assignment but not the ratio.
m5e["stereochemistry"]["mixture"]["evidence_ids"] = [EV[3]]
m5e["stereochemistry"]["assignment"]["evidence_ids"] = [EV[3], EV[7]]
for key in ("chirality_status", "geometry_status", "representation_status"):
    add_evidence(m5e["stereochemistry"][key], EV[7])
m5e["identity_note"] = f(
    "The source gives relative bridgehead descriptors and identifies two N-(CO)O-cholesteryl rotamers, but no complete checked isomeric SMILES was safely encoded from the bounded pages.",
    "Relative bridgehead stereochemistry retained; complete stereochemical graph not encoded.",
    [EV[7]],
)

r5e = find_reaction(f"{DOC}:diyne_5")
r5e["evidence_ids"] = unique(r5e["evidence_ids"] + [EV[7]])
location_5e = copy.deepcopy(r5e["location"]["normalized_value"])
location_5e["pages"] = sorted(set(location_5e.get("pages", []) + [7]))
r5e["location"] = f(
    location_5e,
    location_5e,
    unique(r5e["location"]["evidence_ids"] + [EV[7]]),
)

yield_outcome = next(o for o in r5e["outcomes"] if o["outcome_id"] == "o_yield")
for key in ("product_role", "value", "scope"):
    add_evidence(yield_outcome[key], EV[7])

mass_outcome = st.new_outcome("o_product_mass", "product_amount", ["m_5e"], [EV[7]])
mass_outcome["product_role"] = f("target product", "target product", [EV[7]])
mass_outcome["value"] = f("0.490 g", {"number": 0.490, "unit": "g"}, [EV[7]])
mass_outcome["measurement_basis"] = f(
    "Yield 77% (0.490 g), yellowish solid",
    "isolated product mass reported in the product characterization paragraph",
    [EV[7]],
)
mass_outcome["scope"] = f(
    "compound 5e product characterization",
    {"compound_label": "5e", "single_isolated_product": True},
    [EV[7]],
)
r5e["outcomes"] = [o for o in r5e["outcomes"] if o["outcome_id"] != "o_product_mass"]
r5e["outcomes"].append(mass_outcome)
if not any(n["normalized_value"] == "5e characterization linked to the existing Table 2 entry" for n in r5e["notes"]):
    r5e["notes"].append(
        f(
            "Page 7 reports the full name, isolated yield/mass, physical properties, NMR, and HRMS for 5e; analytical peak lists are characterization, not another synthesis.",
            "5e characterization linked to the existing Table 2 entry",
            [EV[7]],
        )
    )
r5e["review"]["unresolved"] = [
    item for item in r5e["review"]["unresolved"]
    if item != "Product mass is not reported within the bounded pages used for this record."
]


source_context = [
    {
        "context_id": f"{DOC}:characterization_page_7",
        "kind": "characterization",
        "finding": (
            "PDF page 7 contains the end of the 5d 13C NMR/HRMS paragraph and the full 5e name, "
            "yield/mass, physical properties, 1H/13C NMR, and HRMS. The spectra/peak lists characterize "
            "the existing Table 2 experiments and do not create standalone synthesis records."
        ),
        "evidence_ids": [EV[7], IMG[7]],
        "related_reaction_ids": [r5d["reaction_id"], r5e["reaction_id"]],
    },
    {
        "context_id": f"{DOC}:associated_content_page_7",
        "kind": "source_limitation",
        "finding": (
            "Page 7 advertises Supporting Information containing NMR spectra, but source_mapping.json "
            "marks the SI association unestablished and no SI document was supplied; no SI facts were imported."
        ),
        "evidence_ids": [EV[7], IMG[7]],
    },
    {
        "context_id": f"{DOC}:bibliography_pages_7_9",
        "kind": "bibliography",
        "finding": (
            "The remainder of page 7 and all visible content on pages 8-9 are author/funding/acknowledgment "
            "matter and bibliography. Reference titles are citations, not current-paper reaction experiments."
        ),
        "evidence_ids": [EV[7], IMG[7], EV[8], IMG[8], EV[9], IMG[9]],
    },
]

dependency_molecule_ids = {
    "m_azepine_1",
    "m_4d",
    "m_5d",
    "m_4e",
    "m_5e",
    "m_co_acac2_dppe",
    "m_zinc",
    "m_zni2",
    "m_dce",
    "m_petroleum_ether",
    "m_silica_gel",
    "m_ethyl_acetate",
    "m_diyne_template",
    "m_product5_template",
}
dependency_molecules = [
    copy.deepcopy(m) for m in CURRENT["molecules"]
    if m["molecule_id"] in dependency_molecule_ids and m["molecule_id"] != "m_5e"
]

bundle = st.empty(DATASET["documents"])
bundle["molecules"] = dependency_molecules + [m5e]
bundle["reactions"] = [
    copy.deepcopy(next(r for r in CURRENT["reactions"] if r["reaction_id"] == f"{DOC}:current_diyne_generic")),
    r5d,
    r5e,
]
bundle["source_context"] = source_context
st.write(RUN / "submitted.json", bundle)

raw_inventory = {
    "document_id": DOC,
    "bounded_pages": [7, 8, 9],
    "semantic_accuracy_verified": False,
    "inventory": {
        "new_reaction_records": 0,
        "updated_existing_reaction_records": 2,
        "new_molecule_records": 0,
        "updated_existing_molecule_records": 1,
        "standalone_characterization_records": 0,
    },
    "page_findings": {
        "7": [
            "5d 13C NMR/HRMS continuation linked to the existing diyne_4 experiment",
            "5e full characterization linked to the existing diyne_5 experiment",
            "5e isolated product mass 0.490 g added as a product_amount outcome",
            "SI is advertised but not associated or imported",
            "No new reaction arrow, table row, mechanism, or control experiment",
        ],
        "8": ["Bibliography only; no current-paper reaction record"],
        "9": ["Bibliography continuation only; no current-paper reaction record"],
    },
    "limitations": [
        "No explicitly associated SI document was supplied.",
        "Complete stereochemical SMILES for the cholesterol-containing products remain unencoded.",
        "Page 7 analytical peak lists were retained as characterization context rather than synthesis records.",
        "References on pages 7-9 were not converted into experiments without an explicit reaction depiction or procedure on those pages.",
    ],
}
st.write(RUN / "raw_model/pages_0007_0009_extraction.json", raw_inventory)

coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {
            "page": 7,
            "state": "attempted",
            "regions": [
                {
                    "region": "left column top, continuation of compound 5d characterization",
                    "finding": "The 13C NMR list and HRMS for 5d continue from page 6; this is characterization of the existing Table 2 entry, not a new synthesis.",
                    "evidence_ids": [EV[7], IMG[7]],
                },
                {
                    "region": "left column, compound 5e characterization paragraph",
                    "finding": "Full source name for 5e, 77% yield (0.490 g), physical properties, 1H/13C NMR, and HRMS; linked to the existing diyne_5 experiment.",
                    "evidence_ids": [EV[7], IMG[7]],
                },
                {
                    "region": "left column lower half, Associated Content and author information",
                    "finding": "The article advertises NMR Supporting Information, but SI association remains unestablished and no SI facts were imported; author/funding metadata are not reactions.",
                    "evidence_ids": [EV[7], IMG[7]],
                },
                {
                    "region": "right column, acknowledgments and references 1-3",
                    "finding": "Acknowledgments and cited bibliography only; cited titles are not current-paper experimental records on this page.",
                    "evidence_ids": [EV[7], IMG[7]],
                },
            ],
            "record_ids": [r5d["reaction_id"], r5e["reaction_id"]],
            "unresolved": [
                "Complete stereochemical graphs for 5d and 5e were not safely encoded from the bounded source.",
                "The advertised SI was not supplied or explicitly associated.",
            ],
        },
        {
            "page": 8,
            "state": "attempted",
            "regions": [
                {
                    "region": "both columns, references 4-12",
                    "finding": "Bibliographic references only; there are no reaction arrows, tables, procedures, controls, or characterization blocks to record as this paper's experiments.",
                    "evidence_ids": [EV[8], IMG[8]],
                }
            ],
            "record_ids": [],
            "unresolved": [
                "The cited papers were not supplied as source documents, so their experimental details were not extracted from citation titles."
            ],
        },
        {
            "page": 9,
            "state": "attempted",
            "regions": [
                {
                    "region": "upper left, end of reference 12 and reference 13; remainder blank",
                    "finding": "Bibliography continuation only; no reaction arrow, table, procedure, mechanistic experiment, or characterization record is present.",
                    "evidence_ids": [EV[9], IMG[9]],
                }
            ],
            "record_ids": [],
            "unresolved": [
                "The cited book and catalyst paper were not supplied as source documents and were not treated as experiments in the current article."
            ],
        },
    ],
}
st.write(RUN / "coverage/pages_0007_0009.json", coverage)

print(json.dumps({
    "pages": [7, 8, 9],
    "new_reactions": 0,
    "updated_reactions": [r5d["reaction_id"], r5e["reaction_id"]],
    "updated_molecules": [m5e["molecule_id"]],
    "source_context_added": len(source_context),
}, ensure_ascii=False, indent=2))
