#!/usr/bin/env python3
"""Build the page-7--12 v3.1 extraction increment after visual inspection."""
from __future__ import annotations

import copy
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

RUN = ROOT / "runs/current"
DOC = "doc_2a371ced54a752e8"
TITLE = (
    "Bi(cyclopentyl)diol-Derived Boronates in Highly Enantioselective "
    "Chiral Phosphoric Acid-Catalyzed Allylation, Propargylation, and Crotylation of Aldehydes"
)


def ev(page: int, kind: str = "image") -> str:
    return f"{DOC}:page{page}:{kind}"


def evid(*pages: int) -> list[str]:
    out = []
    for page in pages:
        out.extend([ev(page), ev(page, "text")])
    return out


def src(raw, normalized=None, evidence=()):
    return st.fact(raw, raw if normalized is None else normalized, evidence)


def inf(normalized, evidence, why):
    return st.fact(None, normalized, evidence, method="inference", derivation=why)


def derived(normalized, evidence, why):
    return st.fact(None, normalized, evidence, method="structure_derivation", derivation=why)


def missing(evidence=(), status="unreported"):
    return st.fact(evidence=evidence, status=status)


standard = st.read(RUN / "standard_records.json")
dataset = st.read(RUN / "dataset.json")
dataset["documents"][0]["title"] = TITLE
st.write(RUN / "dataset.json", dataset)
existing_molecules = {m["molecule_id"]: m for m in standard["molecules"]}
existing_reactions = {r["reaction_id"]: r for r in standard["reactions"]}
bundle = st.empty(dataset["documents"])
molecules: dict[str, dict] = {}
reactions: list[dict] = []
stereo_diagnostics = {}
new_ids: list[str] = []
updated_ids: list[str] = []


def append_unique(seq: list, values) -> None:
    for value in values:
        if value not in seq:
            seq.append(value)


def add_evidence_to_fact(fact: dict, evidence) -> None:
    append_unique(fact.setdefault("evidence_ids", []), evidence)


def existing_molecule(mid: str) -> dict:
    if mid not in molecules:
        molecules[mid] = copy.deepcopy(existing_molecules[mid])
    return molecules[mid]


def add_occurrence(mid: str, label: str, evidence) -> None:
    m = existing_molecule(mid)
    if any(o["document_id"] == DOC and o["compound_label"].get("normalized_value") == label for o in m["occurrences"]):
        return
    m["occurrences"].append({
        "document_id": DOC,
        "compound_label": src(label, label, evidence),
        "evidence_ids": list(evidence),
    })


def set_candidate_structure(m: dict, smiles: str, evidence, *, chirality=None,
                            geometry=None, representation="encoded", assignment=None,
                            mixture=None, note=None) -> None:
    diag = st.stereo_inspect(smiles)
    stereo_diagnostics[m["molecule_id"]] = diag
    m["identity_status"] = "source_identified"
    m["representation_type"] = "smiles"
    m["smiles"] = derived(
        smiles, evidence,
        "Candidate SMILES transcribed from the inspected source structure/name; chemistry remains needs_review.",
    )
    stdata = m["stereochemistry"]
    stdata["isomeric_smiles"] = derived(
        diag["canonical_isomeric_smiles"], evidence,
        "RDKit canonicalization of the source-transcribed candidate graph.",
    )
    stdata["identity_key"] = derived(
        diag["stereo_graph_key"], evidence,
        "SHA-256 key of the canonical isomeric graph.",
    )
    stdata["atom_centers"] = (
        derived(diag["atom_centers"], evidence, "RDKit diagnostics for the candidate graph.")
        if diag["atom_centers"] else missing(evidence, "not_applicable")
    )
    db = [d for d in diag["double_bonds"] if d["stereo"] != "STEREONONE"]
    stdata["double_bonds"] = (
        derived(db, evidence, "Double-bond geometry encoded in the candidate SMILES.")
        if db else missing(evidence, "not_applicable")
    )
    potential_atoms = [p for p in diag["potential_stereo"] if str(p["type"]).startswith("Atom_")]
    potential_bonds = [p for p in diag["potential_stereo"] if str(p["type"]).startswith("Bond_")]
    if chirality is None:
        chirality = "specified" if diag["atom_centers"] else ("unspecified" if potential_atoms else "achiral")
    if geometry is None:
        geometry = "specified" if potential_bonds and all(p["specified"] == "Specified" for p in potential_bonds) else ("unspecified" if potential_bonds else "not_applicable")
    reason = "Status follows the inspected source identity and standard stereo diagnostics; it is not source-reviewed."
    stdata["chirality_status"] = inf(chirality, evidence, reason)
    stdata["geometry_status"] = inf(geometry, evidence, reason)
    stdata["representation_status"] = inf(representation, evidence, reason)
    if assignment is not None:
        stdata["assignment"] = src(assignment, assignment, evidence)
    if mixture is not None:
        stdata["mixture"] = src(mixture[0], mixture[1], evidence)
    if note:
        m["identity_note"] = src(note, note, evidence)


def add_molecule(mid: str, label: str, name: str, evidence, *, smiles=None,
                 chirality=None, geometry=None, representation=None,
                 assignment=None, mixture=None, note=None, identity_status="source_identified") -> dict:
    if mid in existing_molecules:
        m = existing_molecule(mid)
        for field in (m["occurrences"][0]["compound_label"], m["names"]["original"], m["names"]["english"]):
            add_evidence_to_fact(field, evidence)
    else:
        m = st.new_molecule(mid, DOC, evidence)
        m["identity_status"] = identity_status
        m["occurrences"][0]["compound_label"] = src(label, label, evidence)
        m["names"]["original"] = src(name, name, evidence)
        m["names"]["english"] = src(name, name, evidence)
        if label != name:
            m["names"]["abbreviations"] = src(label, [label], evidence)
        molecules[mid] = m
    if smiles:
        set_candidate_structure(
            m, smiles, evidence, chirality=chirality, geometry=geometry,
            representation=representation or "encoded", assignment=assignment,
            mixture=mixture, note=note,
        )
    elif mid not in existing_molecules or note:
        m["identity_status"] = identity_status
        m["representation_type"] = "unresolved"
        reason = "Structure or stereochemical encoding is not fully available in the inspected source region."
        m["stereochemistry"]["chirality_status"] = inf(chirality or "unresolved", evidence, reason)
        m["stereochemistry"]["geometry_status"] = inf(geometry or "unresolved", evidence, reason)
        m["stereochemistry"]["representation_status"] = inf(representation or "not_encoded", evidence, reason)
        if assignment:
            m["stereochemistry"]["assignment"] = src(assignment, assignment, evidence)
        if mixture:
            m["stereochemistry"]["mixture"] = src(mixture[0], mixture[1], evidence)
        if note:
            m["identity_note"] = src(note, note, evidence)
    return m


def participant(pid, mid, role, evidence, *, amount=None, equiv=None,
                molpct=None, mass=None, volume=None, ref=None, stage="reaction"):
    p = st.new_participant(pid, mid, role, evidence, stage=stage)
    for key, value in (("amount", amount), ("equivalents", equiv), ("mol_percent", molpct),
                       ("mass", mass), ("volume", volume)):
        if value is not None:
            p["charges"][key] = src(value[0], value[1], evidence)
    if ref:
        p["equivalence_reference"] = src(ref, ref, evidence)
    return p


def step(step_id, evidence, pids, operation, *, temp=None, time=None, stage="reaction"):
    s = st.new_step(step_id, evidence, stage=stage)
    s["participant_ids"] = list(pids)
    s["operation"] = src(operation, operation, evidence)
    if temp is not None:
        s["temperature"] = src(temp[0], temp[1], evidence)
    if time is not None:
        s["time"] = src(time[0], time[1], evidence)
    return s


def outcome(oid, kind, target, raw, normalized, evidence, *, basis=None,
            yield_type=None, scope=None):
    o = st.new_outcome(oid, kind, [target], evidence)
    o["value"] = src(raw, normalized, evidence)
    if basis is not None:
        o["measurement_basis"] = src(basis, basis, evidence)
    if yield_type is not None:
        o["yield_type"] = src(yield_type, yield_type, evidence)
    if scope is not None:
        o["scope"] = src(scope, scope, evidence)
    return o


def reaction(rid, evidence, record_type, categories, location, transformation, unresolved=()):
    full = f"{DOC}:{rid}"
    r = st.new_reaction(full, DOC, evidence)
    r["record_type"] = record_type
    r["experiment_categories"] = src(", ".join(categories), categories, evidence)
    pages = sorted({int(e.split(":page")[1].split(":")[0]) for e in evidence})
    r["location"] = src(location, {"pages": pages}, evidence)
    r["transformation_type"] = src(transformation, transformation, evidence)
    r["review"]["unresolved"] = list(unresolved) or ["Machine-extracted record requires expert source review."]
    r["completeness"] = "partial"
    new_ids.append(full)
    return r


def copied_reaction(short_id: str, new_evidence) -> dict:
    full = f"{DOC}:{short_id}"
    r = copy.deepcopy(existing_reactions[full])
    append_unique(r["evidence_ids"], new_evidence)
    add_evidence_to_fact(r["location"], new_evidence)
    pages = set(r["location"]["normalized_value"].get("pages", []))
    pages.update(int(e.split(":page")[1].split(":")[0]) for e in new_evidence)
    r["location"]["normalized_value"]["pages"] = sorted(pages)
    r["review"]["status"] = "needs_review"
    r["review"]["method"] = "model"
    r["review"]["reviewer"] = None
    updated_ids.append(full)
    return r


def add_label(r, evidence, code, category, atoms, change_type="formed", after=1, coverage="partial"):
    cid = "b1"
    change = {
        "change_id": cid,
        "atoms": [{"atom_map": None, "element": atoms[0]}, {"atom_map": None, "element": atoms[1]}],
        "before_order": 0,
        "after_order": after,
        "change_type": change_type,
    }
    r["classification"]["bond_changes"] = [derived(change, evidence, "Net bond change inferred from source reactant/product identities; no atom map is asserted.")]
    r["classification"]["labels"] = derived(
        [{"category_number": category, "code": code, "supporting_change_ids": [cid]}], evidence,
        "Taxonomy label inferred from the listed net bond change.",
    )
    r["classification"]["coverage"] = inf(coverage, evidence, "Complete atom mapping and omitted partners are unresolved.")
    r["classification"]["limitations"] = ["Inferred net change only; no atom-mapped reaction SMILES."]


def classify_carbonyl_addition(r, evidence):
    changes = [
        {"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
        {"change_id": "b2", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "O"}], "before_order": 2, "after_order": 1, "change_type": "order_decreased"},
        {"change_id": "b3", "atoms": [{"atom_map": None, "element": "O"}, {"atom_map": None, "element": "H"}], "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain", "hydrogen_count": 1},
    ]
    r["classification"]["bond_changes"] = [derived(c, evidence, "Net change inferred from the displayed transformation; no atom map is asserted.") for c in changes]
    r["classification"]["labels"] = derived(
        [{"category_number": 2, "code": "C-C", "supporting_change_ids": ["b1"]},
         {"category_number": 24, "code": "O-H", "supporting_change_ids": ["b3"]}],
        evidence, "Multi-label taxonomy assignment inferred from the displayed net transformation.",
    )
    r["classification"]["coverage"] = inf("partial", evidence, "Hydrolysis/proton source and full atom mapping are not explicit.")
    r["classification"]["limitations"] = ["Hydrolysis and proton-transfer partners are omitted.", "No complete atom mapping."]


# Reuse and extend the pre-existing registry for entities cited on pages 7--12.
for mid in [
    "m_product_2", "m_product_3", "m_product_4", "m_product_5", "m_product_6", "m_product_7",
    "m_product_8", "m_product_9", "m_product_10", "m_product_11", "m_product_12", "m_product_13",
    "m_product_14", "m_product_15", "m_product_16", "m_product_17", "m_product_18", "m_product_19",
    "m_product_20", "m_product_21", "m_product_22", "m_product_23", "m_product_24", "m_product_25",
    "m_product_26", "m_product_27", "m_product_28", "m_product_29", "m_product_30", "m_product_31",
    "m_product_32", "m_product_33", "m_product_34", "m_product_36", "m_aldehyde_1",
    "m_aldehyde_for_8", "m_B2", "m_B3", "m_B10", "m_PA1", "m_MS4A", "m_toluene",
    "m_cyclohexane", "m_grubbs_II",
]:
    existing_molecule(mid)
add_occurrence("m_product_2", "51", evid(7))
add_occurrence("m_product_22", "52", evid(7))
add_occurrence("m_product_36", "53", evid(7))

# Application structures and reagents from Scheme 5.
add_molecule(
    "m_valine_46", "46", "(S)-2-((tert-butoxycarbonyl)amino)-3-hydroxy-3-methylbutanoic acid",
    evid(7), smiles="CC(C)(O)[C@H](NC(=O)OC(C)(C)C)C(=O)O", chirality="specified",
    assignment="(S)", note="The complete name is read from the Scheme 5 structure; no CAS number is inferred.",
)
add_molecule(
    "m_product_47", "47", "(1S,2R)-1-(4-bromophenyl)-2-methylbut-3-en-1-yl (S)-2-((tert-butoxycarbonyl)amino)-3-hydroxy-3-methylbutanoate",
    evid(7), smiles="C=C[C@@H](C)[C@H](OC(=O)[C@@H](NC(=O)OC(C)(C)C)C(O)(C)C)c1ccc(Br)cc1",
    chirality="specified", assignment="(1S,2R; amino-acid S)",
)
add_molecule(
    "m_product_49", "49", "(1R,2R)-1-(4-bromophenyl)-2-methylbut-3-en-1-yl acrylate",
    evid(7), smiles="C=C[C@@H](C)[C@@H](OC(=O)C=C)c1ccc(Br)cc1", chirality="specified", assignment="(1R,2R)",
)
add_molecule(
    "m_product_50", "50", "(5R,6R)-6-(4-bromophenyl)-5-methyl-5,6-dihydro-2H-pyran-2-one",
    evid(7), smiles="C[C@H]1[C@H](c2ccc(Br)cc2)OC(=O)C=C1", chirality="specified", assignment="(5R,6R)",
)
add_molecule("m_DIAD", "DIAD", "diisopropyl azodicarboxylate", evid(7), smiles="CC(C)OC(=O)N=NC(=O)OC(C)C",
             geometry="unspecified", representation="partial",
             note="The source does not specify azo-bond geometry; the candidate connectivity therefore remains partial.")
add_molecule("m_PPh3", "PPh3", "triphenylphosphine", evid(7), smiles="P(c1ccccc1)(c1ccccc1)c1ccccc1")
add_molecule("m_acryloyl_chloride", "acryloyl chloride", "acryloyl chloride", evid(7), smiles="C=CC(=O)Cl")
add_molecule("m_triethylamine", "TEA", "triethylamine", evid(7), smiles="CCN(CC)CC")
add_molecule("m_DCM", "DCM", "dichloromethane", evid(7), smiles="ClCCl")
add_molecule("m_product_48", "48", "nannocystin analogue 48 (literature target)", evid(7),
             chirality="unresolved", geometry="unresolved", representation="not_encoded",
             note="Downstream structure is retained as a literature target; it is not converted into a current-paper experiment.")

# Scheme 6 comparison reagents and the distinct product stereochemical states.
add_molecule("m_B11", "B11", "1,3-propanediol-derived allylboronate B11", evid(8), smiles="C=CCB1OCCCO1")
add_molecule("m_B12", "B12", "ethylene-glycol-derived allylboronate B12", evid(8), smiles="C=CCB1OCCO1")
add_molecule("m_B13", "B13", "diisopropyl allylboronate B13", evid(8), smiles="C=CCB(OC(C)C)OC(C)C")
add_molecule("m_product_2_S", "S-2", "(S)-1-phenylbut-3-en-1-ol", evid(8),
             smiles="C=CC[C@H](O)c1ccccc1", chirality="specified", assignment="(S)")
add_molecule(
    "m_product_2_rac", "rac-2", "racemic 1-phenylbut-3-en-1-ol", evid(8),
    smiles="C=CCC(O)c1ccccc1", chirality="racemic", representation="partial",
    mixture=("racemic product", {"type": "racemate", "components": ["R", "S"]}),
    note="A single nonisomeric connectivity SMILES represents the reported racemate; the two enantiomers are not collapsed into an enantiopure identity.",
)

# Preparative-section registry additions/updates.
add_molecule("m_tetrahydroxydiboron", "B2(OH)4", "tetrahydroxydiboron", evid(9), smiles="OB(O)B(O)O")
add_molecule("m_bicyclopentyl_diol", "[1,1'-bicyclopentyl]-1,1'-diol", "[1,1'-bicyclopentyl]-1,1'-diol", evid(9),
             smiles="OC1(CCCC1)C1(O)CCCC1", chirality="achiral")
add_molecule("m_diboron_bicyclopentyl_ester", "diboron [1,1'-bi(cyclopentane)]-1,1'-diol ester",
             "diboron [1,1'-bi(cyclopentane)]-1,1'-diol ester", evid(9),
             chirality="unresolved", geometry="not_applicable", representation="not_encoded",
             note="The source gives a name/formula/characterization but no unambiguous connection table in the inspected region; candidate SMILES is withheld.")
add_molecule("m_allylboronic_acid", "allylboronic acid", "allylboronic acid", evid(3, 9), smiles="C=CCB(O)O")

acetates = {
    "AC1": ("(Z)-3,7-dimethylocta-2,6-dien-1-yl acetate", "CC(=O)OC/C=C(/C)CCC=C(C)C", "specified"),
    "AC2": ("(E)-3,7-dimethylocta-2,6-dien-1-yl acetate", "CC(=O)OC/C=C(\\C)CCC=C(C)C", "specified"),
    "AC3": ("cyclohex-2-en-1-yl acetate", "CC(=O)OC1C=CCCC1", "unspecified"),
    "AC4": ("(E)-hex-2-en-1-yl acetate", "CCC/C=C/COC(C)=O", "specified"),
}
for label, (name, smi, geom) in acetates.items():
    add_molecule(f"m_{label}", label, name, evid(9), smiles=smi,
                 chirality="unspecified" if label == "AC3" else "achiral", geometry=geom,
                 representation="partial" if label == "AC3" else "encoded",
                 assignment=name.split(")", 1)[0] + ")" if name.startswith("(") else None,
                 note="The allylic ring carbon is potentially stereogenic but unspecified by the source." if label == "AC3" else None)
    add_molecule(f"m_precursor_{label}", f"precursor to {label}",
                 f"corresponding allylic alcohol or allylic chloride precursor to {label}", evid(9),
                 chirality="unresolved", geometry="unresolved", representation="not_encoded",
                 identity_status="unresolved",
                 note="The paper gives two precursor classes but does not identify which member was used for this individual acetate.")

boronate_structures = {
    4: ("12-(3-methylbut-2-en-1-yl)-11,13-dioxa-12-boradispiro[4.0.4.6.3.5]tridecane", "CC(C)=CCB1OC2(CCCC2)C2(CCCC2)O1", "not_applicable"),
    5: ("(E)-12-(3,7-dimethylocta-2,6-dien-1-yl)-11,13-dioxa-12-boradispiro[4.0.4.6.3.5]tridecane", "CC(C)=CCC/C(C)=C/CB1OC2(CCCC2)C2(CCCC2)O1", "specified"),
    6: ("(Z)-12-(3,7-dimethylocta-2,6-dien-1-yl)-11,13-dioxa-12-boradispiro[4.0.4.6.3.5]tridecane", "CC(C)=CCC/C(C)=C\\CB1OC2(CCCC2)C2(CCCC2)O1", "specified"),
    7: ("(E)-12-(hex-2-en-1-yl)-11,13-dioxa-12-boradispiro[4.0.4.6.3.5]tridecane", "CCC/C=C/CB1OC2(CCCC2)C2(CCCC2)O1", "specified"),
    8: ("12-(cyclohex-2-en-1-yl)-11,13-dioxa-12-boradispiro[4.0.4.6.3.5]tridecane", "B1(C2C=CCCC2)OC2(CCCC2)C2(CCCC2)O1", "unspecified"),
}
for n, (name, smi, geom) in boronate_structures.items():
    add_molecule(f"m_B{n}", f"B{n}", name, evid(9), smiles=smi,
                 chirality="unspecified" if n == 8 else "achiral", geometry=geom,
                 representation="partial" if n == 8 else "encoded",
                 assignment=name.split(")", 1)[0] + ")" if name.startswith("(") else None,
                 note="The cyclohexenyl attachment center is potentially stereogenic but unspecified by the source." if n == 8 else None)
for n in range(3, 9):
    add_molecule(f"m_precursor_B{n}", f"allylic substrate for B{n}", f"allylic substrate used to prepare B{n}", evid(9),
                 chirality="unresolved", geometry="unresolved", representation="not_encoded", identity_status="unresolved",
                 note="Exact substrate/reagent charges are deferred to a cited procedure and are not reported on page 9.")
add_molecule("m_allenyl_Bpin", "allenylboronic acid pinacol ester", "allenylboronic acid pinacol ester", evid(9),
             smiles="C=C=CB1OC(C)(C)C(C)(C)O1", chirality="achiral", representation="partial",
             note="Terminal-CH2 allene is not axially chiral; allene encoding support is partial.")


# Page 7, Scheme 5a: current Mitsunobu transformation. The 47 -> 48 arrow is literature only.
e = evid(6, 7)
r = copied_reaction("page6_19_to_47", evid(7))
r["reaction_name"] = src("Mitsunobu preparation of 47", "Mitsunobu preparation of 47", e)
r["participants"] = [
    participant("alcohol19", "m_product_19", "reactant", e),
    participant("valine46", "m_valine_46", "reactant", e),
    participant("DIAD", "m_DIAD", "reagent", e),
    participant("PPh3", "m_PPh3", "reagent", e),
    participant("product", "m_product_47", "product", e),
]
r["steps"] = [
    step("low_temperature_charge", evid(7), ["alcohol19", "valine46", "DIAD", "PPh3"], "combine under the depicted Mitsunobu conditions", temp=("−20 °C", {"number": -20, "unit": "°C"})),
    step("warm_to_room_temperature", evid(7), ["alcohol19", "valine46", "DIAD", "PPh3"], "allow the reaction to warm from −20 °C to room temperature", temp=("room temperature after warming from −20 °C", "room temperature")),
]
r["outcomes"] = [
    outcome("yield", "yield", "m_product_47", "86%", {"number": 86, "unit": "%"}, evid(7), yield_type="product yield"),
    outcome("ee", "selectivity", "m_product_47", "99% ee", {"metric": "ee", "value": {"number": 99, "unit": "%"}, "configuration": "(1S,2R; amino-acid S)"}, evid(7)),
    outcome("dr", "selectivity", "m_product_47", "dr >20:1", {"metric": "dr", "ratio": ">20:1", "order": "major:minor", "major": "depicted product 47"}, evid(7)),
]
r["review"]["unresolved"] = ["Reagent amounts, solvent, duration, and all outcome measurement bases are unreported in Scheme 5a."]
add_label(r, evid(7), "C-O", 3, ("C", "O"))
reactions.append(r)

# Page 7, Scheme 5b: acrylation and RCM are distinct experiments.
r = copied_reaction("page6_19_to_49", evid(7))
r["participants"] = [
    participant("alcohol19", "m_product_19", "reactant", e),
    participant("acryloyl_chloride", "m_acryloyl_chloride", "reagent", e),
    participant("TEA", "m_triethylamine", "reagent", e),
    participant("DCM", "m_DCM", "solvent", e),
    participant("product", "m_product_49", "product", e),
]
r["solvent_conditions"]["composition"] = src("DCM", "dichloromethane", evid(7))
r["steps"] = [
    step("low_temperature_charge", evid(7), ["alcohol19", "acryloyl_chloride", "TEA", "DCM"], "combine alcohol 19, acryloyl chloride and triethylamine in DCM", temp=("0 °C", {"number": 0, "unit": "°C"})),
    step("warm_to_room_temperature", evid(7), ["alcohol19", "acryloyl_chloride", "TEA", "DCM"], "allow the reaction to warm from 0 °C to room temperature", temp=("room temperature after warming from 0 °C", "room temperature")),
]
r["outcomes"] = [
    outcome("yield", "yield", "m_product_49", "92%", {"number": 92, "unit": "%"}, evid(7), yield_type="product yield"),
    outcome("ee", "selectivity", "m_product_49", "87% ee", {"metric": "ee", "value": {"number": 87, "unit": "%"}, "configuration": "(1R,2R)"}, evid(7)),
    outcome("dr", "selectivity", "m_product_49", "dr >20:1", {"metric": "dr", "ratio": ">20:1", "order": "major:minor", "major": "depicted product 49"}, evid(7)),
]
r["review"]["unresolved"] = ["Reagent amounts, duration, and all outcome measurement bases are unreported in Scheme 5b."]
add_label(r, evid(7), "C-O", 3, ("C", "O"))
reactions.append(r)

r = copied_reaction("page6_49_to_50", evid(7))
r["participants"] = [
    participant("diene49", "m_product_49", "reactant", e),
    participant("grubbsII", "m_grubbs_II", "catalyst", e),
    participant("DCM", "m_DCM", "solvent", e),
    participant("product", "m_product_50", "product", e),
]
r["solvent_conditions"]["composition"] = src("DCM", "dichloromethane", evid(7))
r["steps"] = [step("reflux", evid(7), ["diene49", "grubbsII", "DCM"], "perform ring-closing metathesis in refluxing DCM", temp=("DCM reflux", "reflux"))]
r["outcomes"] = [
    outcome("yield", "yield", "m_product_50", "84%", {"number": 84, "unit": "%"}, evid(7), yield_type="product yield"),
    outcome("ee", "selectivity", "m_product_50", "82% ee", {"metric": "ee", "value": {"number": 82, "unit": "%"}, "configuration": "(5R,6R)"}, evid(7)),
    outcome("dr", "selectivity", "m_product_50", "dr >20:1", {"metric": "dr", "ratio": ">20:1", "order": "major:minor", "major": "depicted product 50"}, evid(7)),
]
r["notes"] = [src("Page 6 prose rounds/reports the product as 83% ee, whereas Scheme 5b and the page-14 characterization state 82% ee.",
                  "Source conflict retained: 83% ee in page-6 prose versus 82% ee in Scheme 5b.", evid(6, 7))]
r["review"]["unresolved"] = ["Catalyst loading, duration, and all outcome measurement bases are unreported in Scheme 5b.", "The 82%/83% ee source discrepancy is not resolved by the model."]
add_label(r, evid(7), "DB-C-C", 13, ("C", "C"), after=2)
reactions.append(r)

# Scheme 5c gram-scale experiments: retain each row separately.
gram_specs = [
    ("page6_gram_allylation", "m_aldehyde_1", "m_B2", "m_product_2", "m_cyclohexane", 2, "room temperature", "allylation", 97, 99, None),
    ("page6_gram_crotylation", "m_aldehyde_for_8", "m_B3", "m_product_22", "m_toluene", 2, {"number": -30, "unit": "°C"}, "crotylation", 97, 93, ">20:1"),
    ("page6_gram_propargylation", "m_aldehyde_1", "m_B10", "m_product_36", "m_cyclohexane", 5, "room temperature", "propargylation", 97, 98, None),
]
for short_id, aldehyde, boronate, product, solvent, catpct, temp, kind, yp, eep, dr in gram_specs:
    r = copied_reaction(short_id, evid(7))
    p7 = evid(7)
    temp_raw = "−30 °C" if isinstance(temp, dict) else "room temperature"
    r["participants"] = [
        participant("aldehyde", aldehyde, "reactant", e, amount=("10 mmol", {"number": 10, "unit": "mmol"})),
        participant("boronate", boronate, "reactant", e),
        participant("catalyst", "m_PA1", "catalyst", e, molpct=(f"{catpct} mol%", {"number": catpct, "unit": "mol%"}), ref=aldehyde),
        participant("ms", "m_MS4A", "additive", e),
        participant("solvent", solvent, "solvent", e),
        participant("product", product, "product", e),
    ]
    r["reaction_scale"] = src("10 mmol aldehyde", {"number": 10, "unit": "mmol", "reference_molecule_id": aldehyde}, p7)
    r["solvent_conditions"]["composition"] = src("cyclohexane" if solvent == "m_cyclohexane" else "toluene", "cyclohexane" if solvent == "m_cyclohexane" else "toluene", p7)
    r["steps"] = [step("scale_up_hold", p7, ["aldehyde", "boronate", "catalyst", "ms", "solvent"], f"conduct the 10 mmol {kind} shown in Scheme 5c", temp=(temp_raw, temp), time=("12 h", {"number": 12, "unit": "h"}))]
    r["outcomes"] = [
        outcome("yield", "yield", product, f"{yp}%", {"number": yp, "unit": "%"}, p7, yield_type="product yield"),
        outcome("ee", "selectivity", product, f"{eep}% ee", {"metric": "ee", "value": {"number": eep, "unit": "%"}, "configuration": "source product configuration retained"}, p7),
    ]
    if dr:
        r["outcomes"].append(outcome("dr", "selectivity", product, f"dr {dr}", {"metric": "dr", "ratio": dr, "order": "major:minor", "major": "depicted product 52"}, p7))
    r["review"]["unresolved"] = ["Boronate/MS/solvent quantities and outcome measurement bases are not printed in Scheme 5c."]
    if kind == "crotylation":
        conflict = "Page 6 prose collectively says the 10 mmol reactions were conducted with benzaldehyde, but the specific Scheme 5c crotylation row depicts methyl 4-formylbenzoate and product 52 (the product-22 structure). Both source statements are retained."
        if not any(n.get("raw_value") == conflict for n in r["notes"]):
            r["notes"].append(src(conflict, conflict, evid(6, 7)))
        conflict_review = "Reactant conflict: generalized page-6 benzaldehyde prose versus the methyl 4-formylbenzoate structure in Scheme 5c."
        if conflict_review not in r["review"]["unresolved"]:
            r["review"]["unresolved"].append(conflict_review)
    classify_carbonyl_addition(r, p7)
    reactions.append(r)

# Page 8, Scheme 6a is proposed mechanism, not a performed experiment.
r = reaction("scheme6a_proposed_crotylboration_mechanism", evid(8), "proposed_mechanism", ["mechanistic_study"],
             "Scheme 6a", "proposed CPA-catalyzed crotylboration mechanism",
             ["Transition-state/intermediate structures are proposed, not isolated or independently verified experiments."])
r["reaction_name"] = src("Proposed bifunctional chair-like transition-state model", "Proposed bifunctional chair-like transition-state model", evid(8))
r["participants"] = [
    participant("benzaldehyde", "m_aldehyde_1", "reactant", evid(8)),
    participant("B3", "m_B3", "reactant", evid(8)),
    participant("PA1", "m_PA1", "catalyst", evid(8)),
    participant("product17", "m_product_17", "product", evid(8)),
]
r["steps"] = [step("proposed_transition_state", evid(8), ["benzaldehyde", "B3", "PA1"], "proposed axial-chair transition state followed by hydrolysis of intermediate 54")]
classify_carbonyl_addition(r, evid(8))
reactions.append(r)

# Page 8, Scheme 6b: six separate structure-selectivity probe rows.
scheme6_rows = [
    (1, "m_B11", "room temperature", "m_product_2_rac", 0, "racemic"),
    (2, "m_B11", {"number": -78, "unit": "°C"}, "m_product_2_S", 22, "S"),
    (3, "m_B12", "room temperature", "m_product_2_rac", 0, "racemic"),
    (4, "m_B12", {"number": -78, "unit": "°C"}, "m_product_2", 47, "R"),
    (5, "m_B13", "room temperature", "m_product_2_S", 20, "S"),
    (6, "m_B13", {"number": -78, "unit": "°C"}, "m_product_2", 32, "R"),
]
for row, boronate, temp, product, ee_value, config in scheme6_rows:
    e8 = evid(8)
    r = reaction(f"scheme6b_entry_{row}", e8, "specific_experiment", ["mechanistic_study"],
                 f"Scheme 6b entry {row}", "benzaldehyde allylboration backbone probe",
                 ["Reaction duration and product-isolation/yield details are not reported in Scheme 6b."])
    traw = "−78 °C" if isinstance(temp, dict) else "room temperature"
    r["participants"] = [
        participant("benzaldehyde", "m_aldehyde_1", "reactant", e8, amount=("0.10 mmol", {"number": 0.10, "unit": "mmol"})),
        participant("boronate", boronate, "reactant", e8, amount=("0.12 mmol", {"number": 0.12, "unit": "mmol"}), equiv=("1.2 equiv", {"number": 1.2, "unit": "equiv"}), ref="m_aldehyde_1"),
        participant("catalyst", "m_PA1", "catalyst", e8, molpct=("5 mol%", {"number": 5, "unit": "mol%"}), ref="m_aldehyde_1"),
        participant("ms", "m_MS4A", "additive", e8, mass=("50 mg", {"number": 50, "unit": "mg"})),
        participant("solvent", "m_toluene", "solvent", e8, volume=("1 mL", {"number": 1, "unit": "mL"})),
        participant("product", product, "product", e8),
    ]
    r["reaction_scale"] = src("0.10 mmol benzaldehyde", {"number": 0.10, "unit": "mmol", "reference_molecule_id": "m_aldehyde_1"}, e8)
    r["solvent_conditions"]["composition"] = src("toluene", "toluene", e8)
    r["solvent_conditions"]["total_volume"] = src("1 mL", {"number": 1, "unit": "mL"}, e8)
    r["steps"] = [step("backbone_probe", e8, ["benzaldehyde", "boronate", "catalyst", "ms", "solvent"], "perform the Scheme 6b allylboration comparison", temp=(traw, temp))]
    r["outcomes"] = [
        outcome("conversion", "conversion", "m_aldehyde_1", "quant.", "quantitative", e8, basis="1H NMR analysis"),
        outcome(
            "ee", "selectivity", product, f"{ee_value}% ee",
            {"metric": "ee", "value": {"number": ee_value, "unit": "%"}, "configuration": config},
            e8,
            basis=("chiral HPLC analysis; major configuration determined by chiral HPLC and optical-rotation data in the literature"
                   if config in ("R", "S") else "chiral HPLC analysis"),
        ),
    ]
    classify_carbonyl_addition(r, e8)
    reactions.append(r)


# Page 9 preparative experiments.
e9 = evid(9)
r = reaction("section4_2_1_diboron_ester", e9, "specific_experiment", ["other"], "Section 4.2.1",
             "boronic ester formation", ["Absolute reagent scale and yield measurement basis are unreported.", "Product connection table remains unresolved."])
r["participants"] = [
    participant("tetrahydroxydiboron", "m_tetrahydroxydiboron", "reactant", e9, equiv=("1.0 equiv", {"number": 1.0, "unit": "equiv"}), ref="m_tetrahydroxydiboron"),
    participant("diol", "m_bicyclopentyl_diol", "reactant", e9, equiv=("2.0 equiv", {"number": 2.0, "unit": "equiv"}), ref="m_tetrahydroxydiboron"),
    participant("toluene", "m_toluene", "solvent", e9),
    participant("product", "m_diboron_bicyclopentyl_ester", "product", e9),
]
r["solvent_conditions"]["composition"] = src("dry toluene", "toluene", e9)
r["steps"] = [
    step("charge", e9, ["tetrahydroxydiboron", "diol", "toluene"], "place reagents in a flame-dried two-neck flask with stir bar; add dry toluene under argon", temp=("room temperature", "room temperature")),
    step("heat", e9, ["tetrahydroxydiboron", "diol", "toluene"], "stir under argon", temp=("110 °C", {"number": 110, "unit": "°C"}), time=("24 h", {"number": 24, "unit": "h"})),
    step("cool", e9, ["tetrahydroxydiboron", "diol", "toluene"], "cool the reaction mixture to room temperature", temp=("room temperature", "room temperature")),
    step("concentrate", e9, ["tetrahydroxydiboron", "diol", "toluene"], "evaporate solvent and concentrate in vacuo to a white solid", stage="workup"),
]
r["steps"][0]["atmosphere"] = src("argon", "argon", e9)
r["steps"][1]["atmosphere"] = src("argon", "argon", e9)
r["outcomes"] = [outcome("yield", "yield", "m_diboron_bicyclopentyl_ester", "91%", {"number": 91, "unit": "%"}, e9, yield_type="product yield")]
add_label(r, e9, "other-B-O", 41, ("B", "O"))
reactions.append(r)

r = reaction("section4_2_2_B2", evid(3, 9), "specific_experiment", ["other"], "Section 4.2.2",
             "allylboronate ester formation", ["The cited literature procedure supplies conditions not reproduced in the article page; exact charges, solvent, temperature and duration remain unreported."])
r["participants"] = [participant("diol", "m_bicyclopentyl_diol", "reactant", evid(3, 9)), participant("allylboronic_acid", "m_allylboronic_acid", "reactant", evid(3, 9)), participant("product", "m_B2", "product", evid(3, 9))]
r["steps"] = [step("literature_procedure", evid(3, 9), ["diol", "allylboronic_acid"], "prepare B2 following the cited literature procedure for B-allyl-1,3,2-dioxaborinane; isolate a colorless oil")]
r["procedure"]["general_procedure_reference"] = src("literature procedure 16", "literature procedure 16", e9)
r["outcomes"] = [outcome("yield", "yield", "m_B2", "83%", {"number": 83, "unit": "%"}, e9, yield_type="product yield")]
add_label(r, e9, "other-B-O", 41, ("B", "O"))
reactions.append(r)

acetate_yields = {"AC1": (85, "distillation at 108–110 °C, 3 mmHg"), "AC2": (83, "distillation at 109–111 °C, 3 mmHg"), "AC3": (89, "distillation at 49–50 °C, 3 mmHg"), "AC4": (92, "distillation at 59–60 °C, 3 mmHg")}
for label, (yp, purification) in acetate_yields.items():
    r = reaction(f"section4_2_3_{label}", e9, "specific_experiment", ["other"], f"Section 4.2.3.{int(label[-1])}",
                 "allylic acetate preparation", ["The individual precursor class, acetylating reagent, charges, temperature, and duration are not reported on page 9."])
    r["participants"] = [participant("precursor", f"m_precursor_{label}", "reactant", e9), participant("product", f"m_{label}", "product", e9)]
    r["steps"] = [
        step("acetylation", e9, ["precursor"], "prepare by the stated general acetylation of the corresponding allylic alcohol or allylic chloride; individual conditions unreported"),
        step("purification", e9, ["precursor"], purification, stage="purification"),
    ]
    r["outcomes"] = [outcome("yield", "yield", f"m_{label}", f"{yp}%", {"number": yp, "unit": "%"}, e9, yield_type="product yield")]
    add_label(r, e9, "C-O", 3, ("C", "O"))
    reactions.append(r)

boronate_preps = {
    3: (81, "flash chromatography (pentane/ether 50:1)"),
    4: (86, "flash chromatography (pentane/ether 50:1)"),
    5: (81, "flash chromatography (pentane/ether 50:1), then distillation at 176–177 °C, 3 mmHg"),
    6: (78, "flash chromatography (pentane/ether 50:1), then distillation at 174–175 °C, 3 mmHg"),
    7: (81, "flash chromatography (pentane/ether 50:1), then distillation at 155–157 °C, 3 mmHg"),
    8: (72, "flash chromatography (pentane/ether 50:1), then distillation at 90–91 °C, 3 mmHg"),
}
for n, (yp, purification) in boronate_preps.items():
    r = reaction(f"section4_2_3_B{n}", e9, "specific_experiment", ["other"], f"Section 4.2.3.{n + 2}",
                 "allylic boronate preparation", ["The modified Morken procedure is cited but exact substrate/reagent charges, temperature, and duration are not reproduced on page 9."])
    r["participants"] = [participant("allylic_substrate", f"m_precursor_B{n}", "reactant", e9), participant("product", f"m_B{n}", "product", e9)]
    r["steps"] = [
        step("modified_literature_procedure", e9, ["allylic_substrate"], "prepare B3–B8 following the stated modified Morken procedure; exact individual conditions unreported"),
        step("purification", e9, ["allylic_substrate"], purification, stage="purification"),
    ]
    r["procedure"]["general_procedure_reference"] = src("modified Morken literature procedure 17", "literature procedure 17", e9)
    r["outcomes"] = [outcome("yield", "yield", f"m_B{n}", f"{yp}%", {"number": yp, "unit": "%"}, e9, yield_type="product yield")]
    add_label(r, e9, "C-B", 11, ("C", "B"))
    reactions.append(r)

r = reaction("section4_2_4_B10", e9, "specific_experiment", ["other"], "Section 4.2.4",
             "allenylboronate transesterification", ["The cited procedure is not reproduced on page 9; diol amount, solvent, temperature, time and workup remain unreported."])
r["participants"] = [participant("allenyl_Bpin", "m_allenyl_Bpin", "reactant", e9), participant("diol", "m_bicyclopentyl_diol", "reactant", e9), participant("product", "m_B10", "product", e9)]
r["steps"] = [step("literature_procedure", e9, ["allenyl_Bpin", "diol"], "prepare B10 following literature procedure 18 for propanediolboronates; isolate a colorless oil")]
r["procedure"]["general_procedure_reference"] = src("literature procedure 18 (Pei and Krauss)", "literature procedure 18", e9)
r["outcomes"] = [outcome("yield", "yield", "m_B10", "82%", {"number": 82, "unit": "%"}, e9, yield_type="product yield")]
add_label(r, e9, "other-B-O", 41, ("B", "O"))
reactions.append(r)


# Pages 10--12: characterization data supplement existing Scheme 3 records; no duplicate experiments.
hplc = {
    2: (10, "Chiralcel OD-H; hexane/iPrOH 99/1; 0.7 mL/min; tmajor 28.61 min, tminor 33.53 min", 95, 99),
    3: (10, "Chiralpak AD-H; hexane/iPrOH 98/2; 1.0 mL/min; tminor 13.29 min, tmajor 14.36 min", 96, 99),
    4: (10, "Chiralpak AD-H; hexane/iPrOH 98/2; 1.0 mL/min; tmajor 34.00 min, tminor 36.32 min", 96, 99),
    5: (10, "Chiralcel OD-H; hexane/iPrOH 98/2; 1.0 mL/min; tmajor 18.04 min, tminor 23.75 min", 97, 99),
    6: (10, "Chiralpak AD-H; hexane/iPrOH 98/2; 1.0 mL/min; tminor 23.81 min, tmajor 25.75 min", 98, 99),
    7: (10, "Chiralcel AD-H; hexane/iPrOH 95/5; 1.0 mL/min; tmajor 6.56 min, tminor 7.43 min", 94, 99),
    8: (10, "Chiralpak AD-H; hexane/iPrOH 95/5; 1.0 mL/min; tmajor 13.17 min, tminor 14.99 min", 93, 97),
    9: (10, "Chiralcel OD-H; hexane/iPrOH 98/2; 1.0 mL/min; tmajor 21.55 min, tminor 26.97 min", 97, 99),
    10: (10, "Chiralcel OJ-H; hexane/iPrOH 99/1; 1.0 mL/min; tminor 26.12 min, tmajor 30.45 min", 91, 97),
    11: (10, "Chiralcel OD-H; hexane/iPrOH 95/5; 1.0 mL/min; tmajor 13.13 min, tminor 22.64 min", 94, 99),
    12: (11, "Chiralcel OD-H; hexane/iPrOH 95/5; 1.0 mL/min; tmajor 9.21 min, tminor 26.06 min", 92, 99),
    13: (11, "Chiralcel OD-H; hexane/iPrOH 95/5; 1.0 mL/min; tmajor 8.77 min, tminor 13.95 min", 95, 95),
    14: (11, "Chiralcel OD-H; hexane/iPrOH 97/3; 0.5 mL/min; tminor 16.77 min, tmajor 22.07 min", 93, 95),
    15: (11, "Chiralpak AS-H; hexane/iPrOH 97/3; 0.5 mL/min; tminor 20.67 min, tmajor 23.57 min", 96, 93),
    16: (11, "after formation of the 3,5-dinitrobenzoate ester; Chiralcel OD-H; hexane/iPrOH 95/5; 1.0 mL/min; tmajor 9.970 min, tminor 10.666 min", 92, 87),
    17: (11, "Chiralcel OD-H; hexane/iPrOH 90/10; 1.0 mL/min; tminor 5.297 min, tmajor 6.192 min", 97, 95),
    18: (11, "Chiralcel OD-H; hexane/iPrOH 90/10; 1.0 mL/min; tminor 5.166 min, tmajor 5.455 min", 98, 99),
    19: (11, "Chiralcel OD-H; hexane/iPrOH 99/1; 1.0 mL/min; tminor 20.299 min, tmajor 20.823 min", 93, 99),
    20: (11, "Chiralpak AD-H; hexane/iPrOH 95/5; 1.0 mL/min; tmajor 14.132 min, tminor 20.221 min", 97, 99),
    21: (11, "Chiralcel OJ-H; hexane/iPrOH 95/5; 1.0 mL/min; tmajor 17.316 min, tminor 20.829 min", 91, 99),
    22: (11, "Chiralcel OD-H; hexane/iPrOH 95/5; 1.0 mL/min; tmajor 13.187 min, tminor 18.128 min", 94, 99),
    23: (12, "Chiralcel OD-H; hexane/iPrOH 90/10; 1.0 mL/min; tminor 11.641 min, tmajor 12.893 min", 94, 99),
    24: (12, "Chiralcel OD-H; hexane/iPrOH 90/10; 1.0 mL/min; tminor 8.326 min, tmajor 17.596 min", 93, 97),
    25: (12, "Chiralpak AD-H; hexane/iPrOH 95/5; 1.0 mL/min; tminor 20.571 min, tmajor 23.962 min", 92, 93),
    26: (12, "Chiralcel OJ-H; hexane/iPrOH 96/4; 1.0 mL/min; tminor 10.043 min, tmajor 10.569 min", 93, 99),
    27: (12, "Chiralcel OD-H; hexane/iPrOH 90/10; 1.0 mL/min; tmajor 7.074 min, tminor 9.823 min", 92, 95),
    28: (12, "Chiralpak AS-H; hexane/iPrOH 98/2; 1.0 mL/min; tmajor 12.498 min, tminor 15.472 min", 91, 99),
    29: (12, "determined according to ester derivative 30", 90, 91),
    31: (12, "Chiralcel OD-H; hexane/iPrOH 90/10; 1.0 mL/min; tminor 5.150 min, tmajor 6.786 min", 97, 99),
    32: (12, "Chiralpak AD-H; hexane/iPrOH 98/2; 1.0 mL/min; tminor 9.570 min, tmajor 10.661 min", 92, 77),
    33: (12, "Chiralpak AD-H; hexane/iPrOH 98/2; 1.0 mL/min; tminor 9.793 min, tmajor 11.241 min", 87, 95),
    34: (12, "Chiralpak AD-H; hexane/iPrOH 98/2; 1.0 mL/min; tminor 10.884 min, tmajor 11.443 min", 97, 91),
}


def scope_short_id(n: int) -> str:
    if n == 2:
        return "table1_entry_7"
    method = "a" if n <= 12 else "b"
    return f"scheme3_product_{n}_method_{method}"


method_a_text = "flame-dry and argon-backfill a screw-cap tube containing 4 Å MS (50 mg); add PA1, aldehyde and dry cyclohexane; add B2 dropwise over 30 s; stir 15 min at room temperature; load directly on silica and elute with ethyl acetate/hexane 1:9"
method_b_text = "flame-dry and argon-backfill a screw-cap tube containing 4 Å MS (50 mg); add catalyst, aldehyde and dry toluene; cool to −30 °C; add B2 dropwise over 30 s; stir 6 h at −30 °C; load directly on silica and elute with ethyl acetate/hexane 1:9"
substituted_text = "flame-dry and argon-backfill a screw-cap tube containing 4 Å MS (50 mg); add PA1, benzaldehyde and dry toluene; cool to −30 °C; add B3–B8 dropwise over 30 s; stir 6 h; add 1 mL of 1 M HCl and stir 15 min; record crude 1H NMR for dr; flash chromatograph with ethyl acetate/hexanes 1:9"

for n, (page, basis, yp, eep) in hplc.items():
    short = scope_short_id(n)
    r = copied_reaction(short, evid(page))
    pmid = f"m_product_{n}"
    for o in r["outcomes"]:
        if o["kind"] == "yield":
            o["value"] = src(f"{yp}%", {"number": yp, "unit": "%"}, evid(page))
            add_evidence_to_fact(o["yield_type"], evid(page))
        if o["kind"] == "selectivity" and o["outcome_id"] == "ee":
            o["value"] = src(f"{eep}% ee", {"metric": "ee", "value": {"number": eep, "unit": "%"}, "configuration": "source product configuration retained"}, evid(page))
            o["measurement_basis"] = src(basis, basis, evid(page))
    if n == 2:
        precision_note = "Scheme 3 reports >99% ee, while the characterization paragraph reports 99% ee for product 2."
        if not any(n.get("raw_value") == precision_note for n in r["notes"]):
            r["notes"].append(src(precision_note, "Source precision conflict retained: >99% ee versus 99% ee.", evid(4, 10)))
        precision_review = "Product 2 ee precision differs between Scheme 3 (>99%) and characterization (99%)."
        if precision_review not in r["review"]["unresolved"]:
            r["review"]["unresolved"].append(precision_review)
    if n <= 12:
        r["procedure"]["original_text"] = src(method_a_text, method_a_text, evid(10))
        r["procedure"]["addition_order"] = src("PA1, aldehyde and cyclohexane, then B2 dropwise over 30 s", ["PA1", "aldehyde", "cyclohexane", "B2 dropwise over 30 s"], evid(10))
        pids = [p["participant_id"] for p in r["participants"] if p["role"] != "product"]
        r["steps"] = [
            step("prepare_under_argon", evid(10), pids, "evacuate, flame-dry and back-fill the tube with argon", temp=("room temperature", "room temperature")),
            step("add_b2", evid(10), pids, "add B2 dropwise over 30 s", temp=("room temperature", "room temperature"), time=("30 s", {"number": 30, "unit": "s"})),
            step("hold", evid(10), pids, "stir until TLC shows complete aldehyde conversion", temp=("room temperature", "room temperature"), time=("15 min", {"number": 15, "unit": "min"})),
            step("purify", evid(10), pids, "load directly on silica and flash chromatograph with ethyl acetate/hexane 1:9", stage="purification"),
        ]
        r["steps"][0]["atmosphere"] = src("argon", "argon", evid(10))
    elif n in (13, 14, 15, 16):
        r["procedure"]["original_text"] = src(method_b_text, method_b_text, evid(10))
        r["procedure"]["addition_order"] = src("catalyst, aldehyde and toluene; cool to −30 °C; then B2 dropwise over 30 s", ["catalyst", "aldehyde", "toluene", "cool to −30 °C", "B2 dropwise over 30 s"], evid(10))
        pids = [p["participant_id"] for p in r["participants"] if p["role"] != "product"]
        r["steps"] = [
            step("prepare_under_argon", evid(10), pids, "evacuate, flame-dry and back-fill the tube with argon", temp=("room temperature", "room temperature")),
            step("cool", evid(10), pids, "cool the catalyst/aldehyde/toluene mixture", temp=("−30 °C", {"number": -30, "unit": "°C"})),
            step("add_b2", evid(10), pids, "add B2 dropwise over 30 s", temp=("−30 °C", {"number": -30, "unit": "°C"}), time=("30 s", {"number": 30, "unit": "s"})),
            step("hold", evid(10), pids, "stir until TLC shows complete aldehyde conversion", temp=("−30 °C", {"number": -30, "unit": "°C"}), time=("6 h", {"number": 6, "unit": "h"})),
            step("purify", evid(10), pids, "load directly on silica and flash chromatograph with ethyl acetate/hexane 1:9", stage="purification"),
        ]
        r["steps"][0]["atmosphere"] = src("argon", "argon", evid(10))
    elif n in (17, 31, 32, 33, 34):
        r["procedure"]["original_text"] = src(substituted_text, substituted_text, evid(11))
        pids = [p["participant_id"] for p in r["participants"] if p["role"] != "product"]
        r["steps"] = [
            step("prepare_under_argon", evid(11), pids, "evacuate, flame-dry and back-fill the tube with argon", temp=("room temperature", "room temperature")),
            step("cool", evid(11), pids, "cool PA1, benzaldehyde and dry toluene", temp=("−30 °C", {"number": -30, "unit": "°C"})),
            step("add_boronate", evid(11), pids, "add B3–B8 dropwise over 30 s", temp=("−30 °C", {"number": -30, "unit": "°C"}), time=("30 s", {"number": 30, "unit": "s"})),
            step("hold", evid(11), pids, "stir at −30 °C", temp=("−30 °C", {"number": -30, "unit": "°C"}), time=("6 h", {"number": 6, "unit": "h"})),
            step("acid_treatment", evid(11), pids, "add 1 mL of 1 M HCl and stir", time=("15 min", {"number": 15, "unit": "min"}), stage="workup"),
            step("crude_dr_analysis", evid(11), pids, "collect 1H NMR of the crude mixture for dr analysis", stage="characterization"),
            step("purify", evid(11), pids, "flash chromatograph with ethyl acetate/hexanes 1:9", stage="purification"),
        ]
        r["steps"][0]["atmosphere"] = src("argon", "argon", evid(11))
        for o in r["outcomes"]:
            if o["outcome_id"] == "dr":
                o["measurement_basis"] = src("1H NMR analysis of crude product", "1H NMR analysis of crude product", evid(11))
    reactions.append(r)

# Page 12 characterization of the independently retained 29 -> 30 derivatization.
r = copied_reaction("scheme3_29_to_30_derivatization", evid(12))
r["outcomes"] = [
    outcome("yield", "yield", "m_product_30", "71%", {"number": 71, "unit": "%"}, evid(12), yield_type="product yield"),
    outcome("ee", "selectivity", "m_product_30", "91% ee", {"metric": "ee", "value": {"number": 91, "unit": "%"}, "configuration": "(1S,2R)"}, evid(12), basis="Chiralcel IB N-5; hexane/iPrOH 99/1; 1.0 mL/min; tmajor 14.145 min, tminor 15.182 min"),
]
r["review"]["unresolved"] = ["Derivatization reagent, temperature, time and yield measurement basis are unreported in the inspected pages."]
reactions.append(r)


# Source context: characterization and literature-only content are explicitly separated.
bundle["source_context"] = [
    {"context_id": f"{DOC}:page7_47_to_48_literature", "document_id": DOC, "kind": "literature_followup", "evidence_ids": evid(7),
     "finding": "Scheme 5a's 47→48 arrow is attributed to reference 15 and is not a current-paper performed experiment."},
    {"context_id": f"{DOC}:pages9_12_characterization", "document_id": DOC, "kind": "characterization", "evidence_ids": evid(9, 10, 11, 12),
     "finding": "NMR, HRMS, optical rotation, melting-point and retention-time data are characterization context; only synthesis/preparation and analytical selectivity facts were attached to reaction records."},
    {"context_id": f"{DOC}:page9_instrumentation", "document_id": DOC, "kind": "analytical_methods", "evidence_ids": evid(9),
     "finding": "General TLC, chromatography, HPLC, polarimetry, NMR and HRMS instrumentation does not constitute a synthetic experiment."},
    {"context_id": f"{DOC}:page8_conclusion", "document_id": DOC, "kind": "discussion", "evidence_ids": evid(8),
     "finding": "Conclusion statements summarize scope/selectivity and are not duplicated as experiments."},
]

# Incoming validation requires every referenced molecule to be present in this increment.
for r in reactions:
    for p in r["participants"]:
        mid = p["molecule_id"]
        if mid not in molecules:
            molecules[mid] = copy.deepcopy(existing_molecules[mid])
bundle["molecules"] = list(molecules.values())
bundle["reactions"] = reactions

st.write(RUN / "submitted.json", bundle)
st.write(RUN / "raw_model/pages_0007_0012_standard_stereo_checks.json", stereo_diagnostics)
raw = {
    "batch_id": "allene-stereo-v31-200-20261017",
    "document_id": DOC,
    "title": TITLE,
    "visual_pages_inspected": [7, 8, 9, 10, 11, 12],
    "new_record_ids": new_ids,
    "updated_record_ids": updated_ids,
    "new_record_count": len(new_ids),
    "updated_record_count": len(updated_ids),
    "decisions": [
        "Scheme 5 contains three current application experiments and three separate gram-scale rows; 47→48 is literature-only.",
        "Scheme 6a is a proposed mechanism, while each of the six Scheme 6b table rows is a separate mechanistic-study experiment.",
        "Page-9 preparation subsections are independent partial records; cited-procedure conditions were not invented.",
        "Pages 10–12 characterization paragraphs supplement existing scope records instead of creating duplicates.",
        "Mass spectra/NMR peak lists are retained as source context, not synthetic experiments.",
        "The page-6/page-7 gram-crotylation substrate discrepancy is preserved as a source conflict.",
    ],
    "records": [{"reaction_id": r["reaction_id"], "record_type": r["record_type"], "location": r["location"]["raw_value"], "evidence_ids": r["evidence_ids"]} for r in reactions],
}
st.write(RUN / "raw_model/pages_0007_0012_raw_extraction.json", raw)

by_page = {p: [] for p in range(7, 13)}
for r in reactions:
    for page in by_page:
        if ev(page) in r["evidence_ids"] or ev(page, "text") in r["evidence_ids"]:
            by_page[page].append(r["reaction_id"])
coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {"page": 7, "state": "attempted", "regions": [
            {"region": "Scheme 5a", "finding": "19→47 is a current Mitsunobu experiment; 47→48 is explicitly a cited-literature follow-up.", "evidence_ids": evid(7)},
            {"region": "Scheme 5b", "finding": "19→49 acrylation and 49→50 ring-closing metathesis are separate experiments with their printed yields/ee/dr.", "evidence_ids": evid(7)},
            {"region": "Scheme 5c", "finding": "Three separate 10 mmol allylation, crotylation and propargylation rows were transcribed.", "evidence_ids": evid(7)},
        ], "record_ids": by_page[7], "unresolved": ["Most reagent quantities, reaction durations for Scheme 5a/b, and outcome measurement bases are not printed.", "The Scheme 5c crotylation structure conflicts with the generalized page-6 benzaldehyde prose."]},
        {"page": 8, "state": "attempted", "regions": [
            {"region": "Scheme 6a", "finding": "Proposed bifunctional chair-like transition-state pathway; classified as proposed_mechanism, not a performed experiment.", "evidence_ids": evid(8)},
            {"region": "Scheme 6b rows 1–6 and footnotes", "finding": "Six distinct boronate-backbone probe experiments, including two racemic room-temperature outcomes and four enantioenriched outcomes.", "evidence_ids": evid(8)},
            {"region": "conclusion text", "finding": "Summary discussion only; no additional experiment record.", "evidence_ids": evid(8)},
        ], "record_ids": by_page[8], "unresolved": ["Scheme 6b does not report reaction duration or isolated product yields."]},
        {"page": 9, "state": "attempted", "regions": [
            {"region": "general analytical methods", "finding": "Instrumentation and characterization context, not synthetic experiments.", "evidence_ids": evid(9)},
            {"region": "Sections 4.2.1–4.2.4", "finding": "Diboron ester, B2, AC1–AC4, B3–B8 and B10 preparations were recorded individually with printed yields and purification details.", "evidence_ids": evid(9)},
            {"region": "NMR/HRMS blocks", "finding": "Characterization data retained as source context and not converted into reactions.", "evidence_ids": evid(9)},
        ], "record_ids": by_page[9], "unresolved": ["Cited literature procedures omit many individual reagent charges and conditions; these remain missing.", "The diboron ester connection table is not unambiguously encoded from this page."]},
        {"page": 10, "state": "attempted", "regions": [
            {"region": "Sections 4.3 and 4.4", "finding": "Method A/B operations were split into setup, low-temperature/addition, hold and purification steps with paired temperature/time fields.", "evidence_ids": evid(10)},
            {"region": "Sections 4.4.1–4.4.10", "finding": "Products 2–11 characterization yields and chiral-HPLC methods supplement existing scope records.", "evidence_ids": evid(10)},
            {"region": "NMR/HRMS blocks", "finding": "Spectral peak lists are characterization context only.", "evidence_ids": evid(10)},
        ], "record_ids": by_page[10], "unresolved": ["Product 2 is >99% ee in Scheme 3 but 99% ee in the characterization paragraph; both are retained."]},
        {"page": 11, "state": "attempted", "regions": [
            {"region": "Sections 4.4.11–4.4.15", "finding": "Products 12–16 yields and HPLC methods supplement existing records.", "evidence_ids": evid(11)},
            {"region": "Section 4.5 general procedure", "finding": "Benzaldehyde+B3–B8 low-temperature addition, 6 h hold, HCl workup, 15 min treatment, crude-NMR dr analysis and purification were structured as ordered steps.", "evidence_ids": evid(11)},
            {"region": "Sections 4.5.1–4.5.6", "finding": "Products 17–22 yields and HPLC methods supplement existing Scheme 3 records.", "evidence_ids": evid(11)},
        ], "record_ids": by_page[11], "unresolved": ["The Section 4.5 procedure explicitly names benzaldehyde and is not automatically generalized to every other aldehyde scope member."]},
        {"page": 12, "state": "attempted", "regions": [
            {"region": "Sections 4.5.7–4.5.18", "finding": "Products 23–34 yields and HPLC/selectivity methods supplement existing scope records.", "evidence_ids": evid(12)},
            {"region": "Section 4.5.14 compound 30", "finding": "29→30 derivatization record was supplemented with 71% yield and 91% ee/HPLC basis; missing reaction conditions remain missing.", "evidence_ids": evid(12)},
            {"region": "NMR/HRMS blocks", "finding": "Spectral characterization is source context, not additional synthesis records.", "evidence_ids": evid(12)},
        ], "record_ids": by_page[12], "unresolved": ["Product 30 derivatization reagent, temperature and time are not reported on pages 7–12."]},
    ],
}
st.write(RUN / "coverage/pages_0007_0012.json", coverage)
print(json.dumps({"new_records": len(new_ids), "updated_records": len(updated_ids), "molecules_in_increment": len(molecules), "coverage_pages": len(coverage["pages"])}, ensure_ascii=False))
