#!/usr/bin/env python3
"""Build the page 1--6 v3.1 extraction bundle from visually checked source regions."""

from __future__ import annotations

import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st  # noqa: E402


DOC = "doc_60acc8016b4c0614"
DATASET = ROOT / "runs/current/dataset.json"
OUT = ROOT / "runs/current/submitted.json"
RAW = ROOT / "runs/current/raw_model/pages_0001_0006_extraction.json"
COVERAGE = ROOT / "runs/current/coverage/pages_0001_0006.json"


def ev(page: int) -> list[str]:
    return [f"{DOC}:page{page}:text", f"{DOC}:page{page}:image"]


def fact(raw=None, normalized=None, pages=(), method="original_extraction", derivation=None,
         status=None, alternatives=()):
    evidence = [eid for page in pages for eid in ev(page)]
    return st.fact(raw, normalized, evidence, method=method, status=status,
                   derivation=derivation, alternatives=alternatives)


def missing(pages=()):
    return st.fact(evidence=[eid for page in pages for eid in ev(page)])


dataset = st.read(DATASET)
bundle = st.empty(dataset["documents"])
molecules: dict[str, dict] = {}


def add_molecule(mid: str, label: str, name: str, pages: tuple[int, ...], *,
                 stereo="unresolved", smiles=None, template=None, identity_note=None):
    evidence = [eid for page in pages for eid in ev(page)]
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = fact(label, label, pages)
    m["names"]["original"] = fact(name, name, pages)
    m["names"]["english"] = fact(name, name, pages)
    m["names"]["abbreviations"] = missing(pages)
    m["names"]["aliases"] = fact(label, [label], pages)
    m["cas_numbers"] = missing(pages)
    m["external_properties"] = []
    if smiles:
        diag = st.stereo_inspect(smiles)
        m["representation_type"] = "smiles"
        m["smiles"] = fact(
            {"source": "structure/name shown in the article", "label": label},
            diag["canonical_isomeric_smiles"], pages, method="structure_derivation",
            derivation="Candidate SMILES was transcribed from the source drawing/name and canonicalized with RDKit; parsing does not certify source identity.")
        m["generic_template"] = None
        m["stereochemistry"]["atom_centers"] = fact(
            {"source": "drawn structure"}, diag["atom_centers"], pages,
            method="structure_derivation", derivation="RDKit diagnostics on the candidate source-derived SMILES.")
        m["stereochemistry"]["double_bonds"] = fact(
            {"source": "drawn structure"}, diag["double_bonds"], pages,
            method="structure_derivation", derivation="RDKit diagnostics; carbonyl/aromatic double bonds are not E/Z assignments.")
        m["stereochemistry"]["axial"] = fact(
            {"source": "drawn structure"}, [], pages, method="structure_derivation",
            derivation="No allene/biaryl axis occurs in the candidate graph.")
        m["stereochemistry"]["mixture"] = missing(pages)
        m["stereochemistry"]["assignment"] = fact(
            {"source": "drawn structure"}, "achiral candidate graph", pages,
            method="structure_derivation", derivation="No potential stereogenic element was returned for the candidate graph.")
        m["stereochemistry"]["chirality_status"] = fact(
            {"source": "drawn structure"}, "achiral", pages, method="structure_derivation",
            derivation="Candidate graph has no potential atom or axial stereogenic element.")
        m["stereochemistry"]["geometry_status"] = fact(
            {"source": "drawn structure"}, "not_applicable", pages, method="structure_derivation",
            derivation="No alkene for which E/Z is applicable occurs in the candidate graph.")
        m["stereochemistry"]["representation_status"] = fact(
            {"source": "candidate SMILES"}, "encoded", pages, method="structure_derivation",
            derivation="All stereochemical content of this achiral candidate graph is represented in SMILES.")
        m["stereochemistry"]["isomeric_smiles"] = fact(
            {"source": "candidate SMILES"}, diag["canonical_isomeric_smiles"], pages,
            method="structure_derivation", derivation="RDKit canonical isomeric SMILES for the candidate graph.")
        m["stereochemistry"]["identity_key"] = fact(
            None, diag["stereo_graph_key"], pages, method="structure_derivation",
            derivation="SHA-256 stereo graph key calculated by the bundled standard stereo helper.")
    elif template:
        m["representation_type"] = "generic_template"
        m["smiles"] = missing(pages)
        m["generic_template"] = template
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = missing(pages)
        m["generic_template"] = None

    if not smiles:
        assignments = {
            "specified": ("specified", "not_applicable", "not_encoded"),
            "partial": ("partially_specified", "unresolved", "not_encoded"),
            "mixture": ("mixture", "unresolved", "not_encoded"),
            "achiral_template": ("achiral", "not_applicable", "not_encoded"),
            "unresolved": ("unresolved", "unresolved", "unresolved"),
        }
        chiral, geom, rep = assignments[stereo]
        drawing = {"label": label, "source_region": f"PDF page(s) {','.join(map(str, pages))} structure/name"}
        m["stereochemistry"]["chirality_status"] = fact(
            drawing, chiral, pages, method="inference",
            derivation="Status preserves visible wedge/dash, optical-sign, R/S, wavy-bond, or missing assignment information without converting omission into achiral/racemic.")
        m["stereochemistry"]["geometry_status"] = fact(
            drawing, geom, pages, method="inference",
            derivation="E/Z applicability or assignment cannot be fully encoded from the bounded-page drawing.")
        m["stereochemistry"]["representation_status"] = fact(
            drawing, rep, pages, method="inference",
            derivation="No complete source-checked isomeric SMILES is asserted for this complex/generic structure.")
        if stereo in {"specified", "partial", "mixture"}:
            m["stereochemistry"]["assignment"] = fact(
                drawing, "source drawing retains explicit relative/absolute stereochemical marks; full atomwise assignment unresolved",
                pages, method="inference",
                derivation="The source drawing and/or optical/R descriptor was inspected, but atomwise R/S derivation was not completed.")
        if mid in {"m2", "m18", "m23"}:
            m["stereochemistry"]["axial"] = fact(
                drawing, "substituted allene present; no R_a/S_a or P/M assignment stated on pages 1-6",
                pages, method="inference",
                derivation="Allene axis is visible, while the bounded pages provide no explicit axial descriptor.")
        if stereo == "mixture":
            m["stereochemistry"]["mixture"] = fact(
                drawing, "stereoisomeric mixture or unresolved wavy-bond configuration; composition not reported",
                pages, method="inference",
                derivation="The source uses a wavy bond or explicitly calls the outcome a mixture; no 1:1 ratio is assumed.")
    m["identity_note"] = fact(
        identity_note or "Identity retained by source compound label; no CAS number was reported on pages 1-6.",
        identity_note or "Identity retained by source compound label; no CAS number was reported on pages 1-6.", pages)
    molecules[mid] = m
    bundle["molecules"].append(m)
    return m


def generic_template(label: str, fixed: str, pages: tuple[int, ...], variable="R"):
    return {
        "fixed_structure": fact({"source": f"Scheme 2 structure {label}"}, fixed, pages),
        "variables": [{
            "label": variable,
            "attachment_points": fact(f"{variable} attached through the drawn O-{variable} bond", "drawn oxygen substituent", pages),
            "definition": fact(f"{variable} = Ac, CF3CO, or H", "Ac, CF3CO, or H", pages),
            "allowed_values": fact("Ac, CF3CO, or H", ["Ac", "CF3CO", "H"], pages),
            "constraints": fact("identity within each NMR-yield pool was not resolved", "do not collapse to one concrete molecule", pages),
        }],
        "instance_mapping": [],
        "reason": fact("Scheme 2 explicitly labels R = Ac, CF3CO, or H.",
                       "A generic template is required because the reported NMR-yield pool is not a single concrete structure.", pages),
    }


# Source-labelled substrates, intermediates, products, and explicitly used agents.
complex_molecules = [
    ("m1", "1", "(+)-chloropupukeananin", (1, 2, 5, 6), "specified"),
    ("m2", "2", "(+)-iso-A82775C", (1, 2, 4, 5, 6), "partial"),
    ("m3", "3", "pestheic acid", (1, 2), "unresolved"),
    ("m4", "4", "(−)-(R)-maldoxin", (1, 2, 3, 5, 6), "specified"),
    ("m5", "5", "chloropestolide A", (2,), "specified"),
    ("m6", "6", "(+)-chloropestolide B", (2, 5, 6), "specified"),
    ("m7", "7", "chloropestolide C (Si-syn)", (2, 5), "specified"),
    ("m8", "8", "chloropestolide D", (2,), "specified"),
    ("m9", "9", "chloropestolide E", (2,), "specified"),
    ("m10", "10", "chloropestolide F", (2,), "specified"),
    ("m11", "11", "(−)-chloropupukeanolide D", (2, 5, 6), "specified"),
    ("m12", "12", "chloropupukeanolide C", (2,), "specified"),
    ("m13", "13", "chloropestolide H (Re-syn)", (2, 3, 6), "specified"),
    ("m14", "14", "(−)-chloropestolide I (Si-syn)", (2, 3, 6), "specified"),
    ("m15", "15", "(−)-chloropestolide J (Si-anti)", (2, 3, 6), "specified"),
    ("m16", "16", "chloropestolide K (Re-anti)", (2, 3, 6), "specified"),
    ("m17", "17", "siccayne", (2, 3, 6), "unresolved"),
    ("m18", "18", "(+)-model allene 18", (2, 3), "partial"),
    ("m19", "19", "(−)-model diene 19", (2, 3), "specified"),
    ("m20", "20", "Si-syn cycloadduct 20", (2, 3), "specified"),
    ("m21", "21", "(−)-tricyclic model compound 21", (2, 3), "specified"),
    ("m22", "22", "propargylic chloride 22", (4, 6), "specified"),
    ("m23", "23", "protected allene 23", (4, 6), "mixture"),
]
for row in complex_molecules:
    add_molecule(*row[:4], stereo=row[4])

add_molecule("m24", "24", "4-chloro-2,5-dimethoxybenzaldehyde", (4,),
             smiles="O=Cc1c(OC)cc(Cl)c(OC)c1")
add_molecule("m25", "25", "methyl 4-chloro-2,5-dimethoxybenzoate", (4,),
             smiles="COC(=O)c1c(OC)cc(Cl)c(OC)c1")
add_molecule("m26", "26", "methyl 4-chloro-3-hydroxy-2,5-dimethoxybenzoate", (4,),
             smiles="COC(=O)c1c(OC)c(O)c(Cl)c(OC)c1")
add_molecule("m27", "27", "C3 oxidative substitution regioisomer pool 27", (4,),
             stereo="achiral_template", template=generic_template("27", "C3-substituted aromatic O-R regioisomer", (4,)))
add_molecule("m28", "28", "C6 oxidative substitution regioisomer pool 28", (4,),
             stereo="achiral_template", template=generic_template("28", "C6-substituted aromatic O-R regioisomer", (4,)))
add_molecule("mix_chlorohydrins", "chlorohydrins", "unidentified chlorohydrin isomer mixture", (4,),
             stereo="mixture", identity_note="Table 2 footnote: mixture of isomers; structures were not identified.")


agents = {
    "isopropenyl_mgbr": "isopropenylmagnesium bromide",
    "cucn": "copper(I) cyanide",
    "zncl2": "zinc chloride",
    "incl3": "indium(III) chloride",
    "pdpph3": "tetrakis(triphenylphosphine)palladium(0)",
    "organoindium_reagent": "freshly prepared organoindium reagent solution",
    "tbaf": "tetrabutylammonium fluoride",
    "pida": "(diacetoxyiodo)benzene (PIDA)",
    "tfa": "trifluoroacetic acid (TFA)",
    "acoh": "acetic acid",
    "k2co3": "potassium carbonate",
    "h2so4": "sulfuric acid",
    "kotbu": "potassium tert-butoxide",
    "thf": "tetrahydrofuran",
    "etoac": "ethyl acetate",
    "h2o": "water",
    "hfip": "hexafluoroisopropanol",
    "benzene": "benzene",
    "ch2cl2": "dichloromethane",
    "mecn": "acetonitrile",
    "acetone": "acetone",
    "meoh": "methanol",
    "et2o": "diethyl ether",
    "etoh": "ethanol",
    "iproh": "2-propanol",
    "dce": "1,2-dichloroethane",
    "dmf": "N,N-dimethylformamide",
    "nh4cl": "ammonium chloride",
    "nhexane": "n-hexane",
    "mgso4": "magnesium sulfate",
    "brine": "brine",
    "silica": "silica gel",
}
agent_pages = {
    "isopropenyl_mgbr": (4, 6), "cucn": (4,), "zncl2": (4,), "incl3": (4, 6),
    "pdpph3": (4, 6), "organoindium_reagent": (4, 6), "tbaf": (4, 6),
    "pida": (4,), "tfa": (4,), "acoh": (4,), "k2co3": (4,), "h2so4": (4,),
    "kotbu": (5, 6), "thf": (4, 6), "etoac": (3, 6), "h2o": (3,), "hfip": (3,),
    "benzene": (3, 6), "ch2cl2": (3, 4, 5, 6), "mecn": (3, 6), "acetone": (3,),
    "meoh": (3, 4, 6), "et2o": (3,), "etoh": (3,), "iproh": (3,), "dce": (5,),
    "dmf": (5, 6), "nh4cl": (6,), "nhexane": (6,), "mgso4": (6,), "brine": (6,),
    "silica": (6,),
}
for mid, name in agents.items():
    add_molecule(mid, name, name, agent_pages[mid], stereo="unresolved",
                 identity_note="Participant identity is retained by the unambiguous source name; no independent structure/CAS lookup was performed.")


def participant(r, pid, mid, role, pages, *, stage="reaction", mass=None, amount=None,
                volume=None, equiv=None, molpct=None, reference=None, product_role=None):
    p = st.new_participant(pid, mid, role, [eid for page in pages for eid in ev(page)], stage=stage)
    if product_role:
        p["product_role"] = fact(product_role, product_role, pages)
    values = {"mass": mass, "amount": amount, "volume": volume,
              "equivalents": equiv, "mol_percent": molpct}
    for key, value in values.items():
        if value is not None:
            p["charges"][key] = fact(value[0], value[1], pages)
    if reference:
        p["equivalence_reference"] = fact(reference, reference, pages)
    r["participants"].append(p)
    return p


def add_step(r, sid, pages, operation, *, temperature=None, time=None, pressure=None,
             stage="reaction", participants=(), atmosphere=None):
    s = st.new_step(sid, [eid for page in pages for eid in ev(page)], stage=stage)
    s["participant_ids"] = list(participants)
    s["operation"] = fact(operation, operation, pages)
    if temperature is not None:
        s["temperature"] = fact(temperature[0], temperature[1], pages)
    if time is not None:
        s["time"] = fact(time[0], time[1], pages)
    if pressure is not None:
        s["pressure"] = fact(pressure[0], pressure[1], pages)
    if atmosphere is not None:
        s["atmosphere"] = fact(atmosphere, atmosphere, pages)
    r["steps"].append(s)
    return s


def outcome(r, oid, kind, targets, pages, raw, normalized, *, basis=None,
            yield_type=None, role=None, scope=None):
    o = st.new_outcome(oid, kind, targets, [eid for page in pages for eid in ev(page)])
    o["value"] = fact(raw, normalized, pages)
    if basis is not None:
        o["measurement_basis"] = fact(basis, basis, pages)
    if yield_type is not None:
        o["yield_type"] = fact(yield_type, yield_type, pages)
    if role is not None:
        o["product_role"] = fact(role, role, pages)
    if scope is not None:
        o["scope"] = fact(scope, scope, pages)
    r["outcomes"].append(o)
    return o


def classification(r, mode="unresolved"):
    pages = tuple(sorted({int(e.split(":page")[1].split(":")[0]) for e in r["evidence_ids"] if ":page" in e}))
    c = r["classification"]
    c["mapped_reaction_smiles"] = missing(pages)
    c["model_votes"] = []
    if mode == "unresolved":
        c["coverage"] = fact(None, "unresolved", pages, method="inference",
                             derivation="No source-checked atom mapping was available for this branched/complex transformation.")
        c["labels"] = missing(pages)
        c["bond_changes"] = []
        c["limitations"] = ["No atom-mapped reaction SMILES; no bond-formation label is asserted."]
        return
    changes = []
    if mode in {"da", "da_ene", "da_ene_migration", "previous_da_ene"}:
        changes.extend([
            {"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
            {"change_id": "b2", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
            {"change_id": "b3", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 1, "after_order": 2, "change_type": "order_increased"},
        ])
    if mode in {"da_ene", "da_ene_migration", "previous_da_ene"}:
        changes.extend([
            {"change_id": "b4", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
            {"change_id": "b5", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 1, "after_order": 2, "change_type": "order_increased"},
            {"change_id": "b6", "atoms": [{"atom_map": None, "element": "O"}, {"atom_map": None, "element": "H"}], "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain", "hydrogen_count": 1},
        ])
    if mode == "da_ene_migration":
        changes.append({"change_id": "b7", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "O"}], "before_order": 0, "after_order": 1, "change_type": "formed"})
    if mode == "sn2prime":
        changes = [
            {"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
            {"change_id": "b2", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 1, "after_order": 2, "change_type": "order_increased"},
            {"change_id": "b3", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 3, "after_order": 2, "change_type": "order_decreased"},
            {"change_id": "b4", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "Cl"}], "before_order": 1, "after_order": 0, "change_type": "broken"},
        ]
    if mode in {"hydroxylation", "migration"}:
        changes = [{"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "O"}], "before_order": 0, "after_order": 1, "change_type": "formed"}]
    labels = {}
    for change in changes:
        if change["change_type"] not in {"formed", "order_increased", "hydrogen_gain"}:
            continue
        label = st.bond_label(change["atoms"][0]["element"], change["atoms"][1]["element"], change["after_order"])
        labels.setdefault(label["code"], {**label, "supporting_change_ids": []})["supporting_change_ids"].append(change["change_id"])
    c["bond_changes"] = [fact({"source": "reaction drawing and named transformation"}, ch, pages,
                                      method="structure_derivation",
                                      derivation="Net bond-type change inferred from the named/drawn transformation; atom identities are not mapped.") for ch in changes]
    c["labels"] = fact({"source": "derived from partial bond changes"}, list(labels.values()), pages,
                       method="structure_derivation",
                       derivation="Taxonomy labels were calculated from the partial source-derived bond-change list.")
    c["coverage"] = fact(None, "partial", pages, method="inference",
                         derivation="Connectivity class is clear, but complex structures were not atom-mapped and spectator/coproduct changes may be absent.")
    c["limitations"] = ["No atom-mapped reaction SMILES.", "Bond changes are transformation-level and partial, not atom-specific.", "Implicit-H/proton source is unresolved where present."]


def reaction(suffix, pages, record_type, categories, name, transform, *, mode="unresolved", location=None):
    evidence = [eid for page in pages for eid in ev(page)]
    r = st.new_reaction(f"{DOC}:{suffix}", DOC, evidence)
    r["record_type"] = record_type
    r["experiment_categories"] = fact(categories, categories, pages)
    r["location"] = fact(location or {"pdf_pages": list(pages)}, location or {"pdf_pages": list(pages)}, pages)
    r["reaction_name"] = fact(name, name, pages)
    r["transformation_type"] = fact(None, transform, pages, method="inference",
                                     derivation="Normalized from the source reaction title/text; not an independent mechanism claim.")
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                   "checked_fields": [], "unresolved": []}
    r["completeness"] = "partial"
    classification(r, mode)
    bundle["reactions"].append(r)
    return r


page_records: dict[int, list[str]] = {i: [] for i in range(1, 7)}


def index(r, *pages):
    for page in pages:
        page_records[page].append(r["reaction_id"])


# Page 2: proposed biosynthetic network, and page 3: clearly labelled previous work.
r = reaction("figure2_proposed_biosynthesis", (2,), "proposed_mechanism", ["other"],
             "proposed biosynthetic pathway of the chloropupukeananin family",
             "proposed oxidative dearomatization / Diels–Alder / carbonyl-ene / acyl migration network",
             location={"page": 2, "figure": "Figure 2"})
for pid, mid, role in [("p3", "m3", "reactant"), ("p2", "m2", "reactant"), ("p4", "m4", "reactant"),
                       ("p1", "m1", "product"), ("p5", "m5", "product"), ("p6", "m6", "product"),
                       ("p7", "m7", "product"), ("p8", "m8", "product"), ("p9", "m9", "product"),
                       ("p10", "m10", "product"), ("p11", "m11", "product"), ("p12", "m12", "product")]:
    participant(r, pid, mid, role, (2,), product_role="proposed pathway product" if role == "product" else None)
r["procedure"]["original_text"] = fact("Figure 2. Proposed biosynthetic pathway of the chloropupukeananin family.",
                                          "proposed pathway only; not a reported experiment", (2,))
r["review"]["unresolved"] = ["Branched proposed pathway is not atom-mapped and must not be treated as experimental evidence."]
index(r, 2)

r = reaction("scheme1_previous_da_ene", (2, 3), "literature_example", ["other"],
             "Diels–Alder/carbonyl-ene cascade reaction in our previous study",
             "intermolecular Diels–Alder/carbonyl-ene cascade", mode="previous_da_ene",
             location={"page": 3, "scheme": "Scheme 1", "caption": "Our Previous Study"})
participant(r, "s18", "m18", "reactant", (3,))
participant(r, "s19", "m19", "reactant", (3,))
participant(r, "i20", "m20", "product", (3,), product_role="proposed/depicted cycloadduct intermediate")
participant(r, "p21", "m21", "product", (3,), product_role="reported product")
participant(r, "solv", "ch2cl2", "solvent", (3,))
add_step(r, "high_pressure_hold", (3,), "subjected to high pressure in CH2Cl2",
         temperature=("25 °C", {"number": 25, "unit": "°C"}),
         time=("96 h", {"number": 96, "unit": "h"}),
         pressure=("1.0 GPa", "1.0 GPa"), participants=("s18", "s19", "solv"))
outcome(r, "yield_21", "yield", ["m21"], (2, 3), "70%", {"number": 70, "unit": "%"},
        role="tricyclic product", scope={"source": "previous study, ref 14a"})
r["review"]["unresolved"] = ["Prior-work example; measurement basis and full charges are not reported on pages 2-3.",
                              "Complex structures and complete atom mapping remain unresolved."]
index(r, 2, 3)


# Table 1: every current-paper entry is an independent experiment. The literature isolation ratio is source context, not an entry.
table1 = [
    (1, "neat", None, 24, [9, 32, 51, 8]),
    (2, "EtOAc", "etoac", 48, [6, 16, 71, 7]),
    (3, "MeOH/H2O (4/1)", ("meoh", "h2o"), 48, [4, 11, 74, 11]),
    (4, "HFIP", "hfip", 48, [14, 37, 34, 15]),
    (5, "benzene", "benzene", 48, [9, 42, 44, 5]),
    (6, "CH2Cl2", "ch2cl2", 48, [14, 44, 36, 6]),
    (7, "CH3CN", "mecn", 48, [10, 28, 53, 9]),
    (8, "acetone", "acetone", 48, [6, 20, 64, 10]),
    (9, "MeOH", "meoh", 48, [4, 11, 75, 10]),
    (10, "Et2O", "et2o", 48, [6, 17, 69, 7]),
    (11, "EtOH", "etoh", 48, [3, 7, 80, 10]),
    (12, "i-PrOH", "iproh", 48, [3, 5, 81, 11]),
    (13, "MeOH/CH2Cl2 (10/1), 1.0 GPa", ("meoh", "ch2cl2"), 1, [3, 11, 75, 11]),
    (14, "CH2Cl2, 1.0 GPa", "ch2cl2", 1, [15, 46, 33, 6]),
]
for entry, cond, solvents, hours, ratios in table1:
    pages = (3, 6) if 2 <= entry <= 12 else (3,)
    r = reaction(f"table1_entry_{entry:02d}", pages, "specific_experiment", ["optimization"],
                 "intermolecular Diels–Alder reaction between (−)-4 and 17",
                 "intermolecular Diels–Alder cycloaddition", mode="da",
                 location={"page": 3, "table": "Table 1", "entry": entry})
    eq17 = 1.1 if entry in (13, 14) else 2.0
    participant(r, "s17", "m17", "reactant", pages,
                mass=("34.8 mg", {"number": 34.8, "unit": "mg"}) if 2 <= entry <= 12 else None,
                amount=("0.200 mmol", {"number": 0.200, "unit": "mmol"}) if 2 <= entry <= 12 else None,
                equiv=(f"{eq17} equiv", {"number": eq17, "unit": "equiv"}), reference="m4")
    participant(r, "s4", "m4", "reactant", pages,
                mass=("38.0 mg", {"number": 38.0, "unit": "mg"}) if 2 <= entry <= 12 else None,
                amount=("0.100 mmol", {"number": 0.100, "unit": "mmol"}) if 2 <= entry <= 12 else None,
                equiv=("1.0 equiv", {"number": 1.0, "unit": "equiv"}), reference="m4")
    for pid, mid, role in [("p13", "m13", "minor product"), ("p14", "m14", "product"),
                           ("p15", "m15", "product"), ("p16", "m16", "minor product")]:
        participant(r, pid, mid, "product", pages, product_role=role)
    solvent_pids = []
    if isinstance(solvents, tuple):
        for i, smid in enumerate(solvents, 1):
            pid = f"solv{i}"
            participant(r, pid, smid, "solvent", pages)
            solvent_pids.append(pid)
    elif solvents:
        participant(r, "solv", solvents, "solvent", pages,
                    volume=("0.50 mL", {"number": 0.50, "unit": "mL"}) if 2 <= entry <= 12 else None)
        solvent_pids.append("solv")
    r["solvent_conditions"]["composition"] = fact(cond, cond, pages)
    if "/" in cond and "(" in cond:
        ratio = cond.split("(")[-1].split(")")[0]
        r["solvent_conditions"]["mixture_ratio"] = fact(ratio, ratio, pages)
    if 2 <= entry <= 12:
        r["solvent_conditions"]["total_volume"] = fact("0.50 mL", {"number": 0.50, "unit": "mL"}, pages)
        r["reaction_scale"] = fact("(−)-4: 0.100 mmol", {"reference_molecule_id": "m4", "amount": {"number": 0.100, "unit": "mmol"}}, pages)
    pressure = ("1.0 GPa", "1.0 GPa") if entry in (13, 14) else None
    add_step(r, "reaction", pages, f"stirred under {cond}",
             temperature=("rt", "r.t."), time=(f"{hours} h", {"number": hours, "unit": "h"}),
             pressure=pressure, participants=("s17", "s4", *solvent_pids))
    if 2 <= entry <= 12:
        add_step(r, "concentration", pages, "reaction mixture was concentrated in vacuo",
                 stage="workup", participants=())
    outcome(r, "product_ratio", "selectivity", ["m13", "m14", "m15", "m16"], pages,
            f"13/14/15/16 = {':'.join(map(str, ratios))}",
            {"metric": "product_ratio", "order": ["13", "14", "15", "16"], "ratio": ratios},
            basis="1H NMR", scope={"stereochemical_order": ["Re-syn", "Si-syn", "Si-anti", "Re-anti"]})
    if entry == 1:
        outcome(r, "mixture_yield", "yield", ["m13", "m14", "m15", "m16"], pages,
                "quantitative yield", "quantitative", role="combined cycloadduct mixture")
    if entry == 6:
        outcome(r, "isolated_yield_14", "yield", ["m14"], pages, "39%", {"number": 39, "unit": "%"},
                basis="isolated after preparative TLC", yield_type="isolated")
        outcome(r, "isolated_yield_15", "yield", ["m15"], pages, "34%", {"number": 34, "unit": "%"},
                basis="isolated after preparative TLC", yield_type="isolated")
        outcome(r, "isolated_yield_13_16", "yield", ["m13", "m16"], pages, "17%", {"number": 17, "unit": "%"},
                basis="isolated inseparable mixture after preparative TLC", yield_type="isolated")
        outcome(r, "isolated_mixture_mass", "product_amount", ["m13", "m14", "m15", "m16"], pages,
                "54.8 mg (99%) before preparative TLC", {"mass": {"number": 54.8, "unit": "mg"}, "yield": {"number": 99, "unit": "%"}},
                basis="isolated mixture after silica gel chromatography")
        for pid, mid, role in [("pur_hex", "nhexane", "purification_agent"), ("pur_ea", "etoac", "purification_agent"),
                               ("pur_bz", "benzene", "purification_agent"), ("pur_acn", "mecn", "purification_agent"),
                               ("pur_silica", "silica", "purification_agent")]:
            participant(r, pid, mid, role, pages, stage="purification")
        add_step(r, "purification", pages,
                 "silica gel chromatography (hexane/EtOAc 2/1 to 1/1), then preparative TLC (benzene/acetonitrile 9/1)",
                 stage="purification", participants=("pur_hex", "pur_ea", "pur_bz", "pur_acn", "pur_silica"))
    r["procedure"]["original_text"] = fact(
        f"Table 1 entry {entry}: {cond}; {hours} h; ratio {':'.join(map(str, ratios))}.",
        f"Table 1 entry {entry}", pages)
    r["review"]["unresolved"] = ["Complex reactant/product structures are retained by source labels; atomwise stereochemistry and mapping remain unresolved."]
    if entry not in (1, 6):
        r["review"]["unresolved"].append("Only product ratio is reported; absolute product yield and product-specific isolation are unreported.")
    index(r, 3)
    if 2 <= entry <= 12:
        index(r, 6)


# Table 2: entry 1 is explicitly previous work; entries 2--6 are current optimization experiments.
table2 = [
    (1, "cucn", 5.0, 2.0, None, 30, 59, None, "literature_example"),
    (2, "zncl2", 2.0, 2.0, None, 48, 18, 16, "specific_experiment"),
    (3, "incl3", 5.5, 2.0, 10.0, 78, None, 15, "specific_experiment"),
    (4, "incl3", 4.0, 2.0, 10.0, 85, None, 3, "specific_experiment"),
    (5, "incl3", 2.0, 2.0, 10.0, 21, 65, 2, "specific_experiment"),
    (6, "incl3", 2.4, 1.2, 5.0, 93, None, 3, "specific_experiment"),
]
for entry, salt, grignard_eq, salt_eq, pdmol, y23, rec22, chlor, rtype in table2:
    pages = (4, 6) if entry == 6 else (4,)
    r = reaction(f"table2_entry_{entry:02d}", pages, rtype,
                 ["optimization"] if rtype == "specific_experiment" else ["other"],
                 "Pd-catalyzed anti-SN2′ reaction using organometallic reagent",
                 "anti-SN2′ substitution / propargyl-to-allene rearrangement", mode="sn2prime",
                 location={"page": 4, "table": "Table 2", "entry": entry})
    participant(r, "s22", "m22", "reactant", pages,
                mass=("282 mg", {"number": 282, "unit": "mg"}) if entry == 6 else None,
                amount=("0.58 mmol", {"number": 0.58, "unit": "mmol"}) if entry == 6 else None)
    participant(r, "grignard", "isopropenyl_mgbr", "reagent", pages,
                amount=("1.6 mmol", {"number": 1.6, "unit": "mmol"}) if entry == 6 else None,
                volume=("3.2 mL of 0.5 M solution", {"number": 3.2, "unit": "mL"}) if entry == 6 else None,
                equiv=(f"{grignard_eq} equiv", {"number": grignard_eq, "unit": "equiv"}), reference="m22")
    participant(r, "salt", salt, "reagent", pages,
                mass=("177 mg", {"number": 177, "unit": "mg"}) if entry == 6 else None,
                amount=("0.80 mmol", {"number": 0.80, "unit": "mmol"}) if entry == 6 else None,
                equiv=(f"{salt_eq} equiv", {"number": salt_eq, "unit": "equiv"}), reference="m22")
    if pdmol is not None:
        participant(r, "pd", "pdpph3", "catalyst", pages,
                    mass=("35 mg", {"number": 35, "unit": "mg"}) if entry == 6 else None,
                    amount=("0.030 mmol", {"number": 0.030, "unit": "mmol"}) if entry == 6 else None,
                    molpct=(f"{pdmol} mol %", {"number": pdmol, "unit": "mol%"}), reference="m22")
    participant(r, "p23", "m23", "product", pages, product_role="allene product")
    participant(r, "byprod", "mix_chlorohydrins", "product", pages, product_role="unidentified byproduct mixture")
    participant(r, "thf", "thf", "solvent", pages,
                volume=("4.0 mL", {"number": 4.0, "unit": "mL"}) if entry == 6 else None)
    temp = ("−40 °C", {"number": -40, "unit": "°C"}) if entry == 1 else ("rt", "r.t.")
    time_value = ("5 h", {"number": 5, "unit": "h"}) if entry == 6 else None
    if entry != 6:
        add_step(r, "anti_sn2prime", pages, "anti-SN2′ substitution in THF", temperature=temp,
                 time=time_value, participants=tuple(p["participant_id"] for p in r["participants"] if p["stage"] == "reaction"))
    outcome(r, "nmr_yield_23", "yield", ["m23"], pages, f"{y23}%", {"number": y23, "unit": "%"},
            basis="NMR", yield_type="NMR", role="allene 23")
    if rec22 is not None:
        outcome(r, "recovery_22", "recovery", ["m22"], pages, f"{rec22}%", {"number": rec22, "unit": "%"},
                basis="NMR", role="recovered starting material 22")
    if chlor is not None:
        outcome(r, "nmr_yield_chlorohydrins", "yield", ["mix_chlorohydrins"], pages, f"{chlor}%", {"number": chlor, "unit": "%"},
                basis="NMR", yield_type="NMR", role="unidentified chlorohydrin isomers")
    if entry == 6:
        participant(r, "thf_prep", "thf", "solvent", pages,
                    volume=("2.2 mL", {"number": 2.2, "unit": "mL"}))
        participant(r, "organoindium", "organoindium_reagent", "reagent", pages,
                    amount=("0.71 mmol", {"number": 0.71, "unit": "mmol"}),
                    volume=("4.7 mL of 0.15 M solution", {"number": 4.7, "unit": "mL"}),
                    equiv=("1.2 equiv", {"number": 1.2, "unit": "equiv"}), reference="m22")
        add_step(r, "organoindium_addition", pages,
                 "isopropenylmagnesium bromide solution was added to dried InCl3 in THF",
                 temperature=("0 °C", {"number": 0, "unit": "°C"}),
                 participants=("grignard", "salt", "thf_prep"))
        add_step(r, "organoindium_hold", pages,
                 "mixture was warmed to room temperature and stirred",
                 temperature=("room temperature", "r.t."), time=("30 min", {"number": 30, "unit": "min"}),
                 participants=("grignard", "salt", "thf_prep", "organoindium"))
        add_step(r, "coupling_addition", pages,
                 "freshly prepared organoindium reagent solution was added to 22 and Pd(PPh3)4 in THF",
                 temperature=("0 °C", {"number": 0, "unit": "°C"}),
                 participants=("s22", "pd", "thf", "organoindium"))
        add_step(r, "coupling_hold", pages,
                 "coupling mixture was stirred at room temperature",
                 temperature=("room temperature", "r.t."), time=("5 h", {"number": 5, "unit": "h"}),
                 participants=("s22", "pd", "thf", "organoindium", "p23", "byprod"))
        participant(r, "tbaf", "tbaf", "reagent", pages,
                    volume=("1.60 mL of 1.0 M solution", {"number": 1.60, "unit": "mL"}),
                    amount=("1.60 mmol", {"number": 1.60, "unit": "mmol"}),
                    equiv=("3.0 equiv", {"number": 3.0, "unit": "equiv"}), reference="m22")
        participant(r, "thf_desilylation", "thf", "solvent", pages,
                    volume=("7.0 mL", {"number": 7.0, "unit": "mL"}))
        participant(r, "p2", "m2", "product", pages, product_role="two-step deprotected product")
        add_step(r, "desilylation", pages, "TBAF was added to crude 23 in THF at 0 °C and stirred",
                 temperature=("0 °C", {"number": 0, "unit": "°C"}),
                 time=("30 min", {"number": 30, "unit": "min"}), participants=("p23", "tbaf", "thf_desilylation", "p2"))
        outcome(r, "isolated_yield_2", "yield", ["m2"], pages, "140 mg, 92% (2 steps)", {"number": 92, "unit": "%"},
                basis="isolated after flash silica-gel chromatography", yield_type="isolated, two steps",
                role="(+)-iso-A82775C", scope={"start": "22", "end": "2", "product_mass": {"number": 140, "unit": "mg"}})
        for pid, mid, role in [("quench1", "nh4cl", "workup_agent"), ("hex", "nhexane", "workup_agent"),
                               ("dry", "mgso4", "workup_agent"), ("brine", "brine", "workup_agent"),
                               ("ea", "etoac", "workup_agent"), ("silica", "silica", "purification_agent")]:
            participant(r, pid, mid, role, pages, stage="workup" if role == "workup_agent" else "purification")
        add_step(r, "workup_purification", pages,
                 "NH4Cl quench/extraction, brine workup, MgSO4 drying, and flash silica-gel chromatography (n-hexane/EtOAc 3/1)",
                 stage="workup", participants=("quench1", "hex", "dry", "brine", "ea", "silica"))
        r["reaction_scale"] = fact("22: 0.58 mmol", {"reference_molecule_id": "m22", "amount": {"number": 0.58, "unit": "mmol"}}, pages)
        r["notes"].append(fact(
            "Organoindium preparation: InCl3 177 mg (0.80 mmol) in THF 2.2 mL; isopropenylmagnesium bromide 0.5 M, 3.2 mL (1.6 mmol); 4.7 mL of the 0.15 M organoindium solution (0.71 mmol, 1.2 equiv) was then added to 22.",
            "The preparation batch charges and the aliquot added to 22 have different scopes and are retained separately; no silent equivalent recalculation was made.",
            pages, method="inference",
            derivation="Interpretation of the explicit preparation-batch and reaction-aliquot quantities on page 6."))
    r["review"]["unresolved"] = ["Complex source structures are not represented as source-checked atom-mapped SMILES."]
    if entry != 6:
        r["review"]["unresolved"].append("Reaction time and full procedure are unreported for this optimization row on pages 1-6.")
    index(r, 4)
    if entry == 6:
        index(r, 6)


# Scheme 2: two-step aldehyde-to-ester sequence.
r = reaction("scheme2_24_to_25", (4,), "specific_experiment", ["other"],
             "improved synthesis of phenol 26: preparation of salicylate 25",
             "Pinnick oxidation followed by Fischer esterification", location={"page": 4, "scheme": "Scheme 2", "arrow": "24 to 25"})
participant(r, "s24", "m24", "reactant", (4,))
participant(r, "acid", "h2so4", "reagent", (4,))
participant(r, "meoh", "meoh", "solvent", (4,))
participant(r, "p25", "m25", "product", (4,), product_role="two-step product")
add_step(r, "pinnick_oxidation", (4,), "Pinnick oxidation", participants=("s24",))
add_step(r, "fischer_esterification", (4,), "H2SO4/MeOH esterification", participants=("s24", "acid", "meoh", "p25"))
outcome(r, "yield_25", "yield", ["m25"], (4,), "81% (2 steps)", {"number": 81, "unit": "%"},
        yield_type="two-step yield", scope={"start": "24", "end": "25"})
r["review"]["unresolved"] = ["Pinnick reagents, quantities, times, temperatures, and yield measurement basis are not given in Scheme 2."]
index(r, 4)


# Scheme 2 regioselectivity screen, each row independent.
screen = [("TFA", "tfa", None, 50, 23), ("TFA/AcOH (20/1)", ("tfa", "acoh"), "20/1", 70, 17),
          ("TFA/AcOH (10/1)", ("tfa", "acoh"), "10/1", 84, 12),
          ("TFA/AcOH (5/1)", ("tfa", "acoh"), "5/1", 72, 12)]
for i, (cond, solvents, ratio, y27, y28) in enumerate(screen, 1):
    r = reaction(f"scheme2_screen_{i:02d}", (4,), "specific_experiment", ["optimization"],
                 "Kita oxidative nucleophilic substitution regioselectivity screen",
                 "oxidative aromatic C–H oxygenation", mode="hydroxylation",
                 location={"page": 4, "scheme": "Scheme 2 inset table", "row": cond})
    participant(r, "s25", "m25", "reactant", (4,))
    participant(r, "pida", "pida", "reagent", (4,), equiv=("1.2 equiv", {"number": 1.2, "unit": "equiv"}), reference="m25")
    mids = solvents if isinstance(solvents, tuple) else (solvents,)
    spids = []
    for j, smid in enumerate(mids, 1):
        pid = f"solv{j}"
        participant(r, pid, smid, "solvent", (4,))
        spids.append(pid)
    participant(r, "p27", "m27", "product", (4,), product_role="C3 regioisomer pool")
    participant(r, "p28", "m28", "product", (4,), product_role="C6 regioisomer pool")
    r["solvent_conditions"]["composition"] = fact(cond, cond, (4,))
    if ratio:
        r["solvent_conditions"]["mixture_ratio"] = fact(ratio, ratio, (4,))
    add_step(r, "oxidative_substitution", (4,), "PIDA oxidative nucleophilic substitution",
             temperature=("0 °C", {"number": 0, "unit": "°C"}), participants=("s25", "pida", *spids, "p27", "p28"))
    outcome(r, "nmr_yield_27", "yield", ["m27"], (4,), f"{y27}%", {"number": y27, "unit": "%"},
            basis="NMR", yield_type="NMR", role="regioisomer pool 27")
    outcome(r, "nmr_yield_28", "yield", ["m28"], (4,), f"{y28}%", {"number": y28, "unit": "%"},
            basis="NMR", yield_type="NMR", role="regioisomer pool 28")
    r["review"]["unresolved"] = ["27 and 28 each cover R = Ac, CF3CO, or H; individual composition and structures are unresolved.",
                                  "Reaction time, scale, and workup are not stated in the page-4 screen."]
    index(r, 4)


r = reaction("scheme2_preparative_25_to_26", (4,), "specific_experiment", ["scale_up"],
             "one-pot oxidative substitution/hydrolysis to phenol 26",
             "oxidative aromatic C–H oxygenation followed by hydrolysis", mode="hydroxylation",
             location={"page": 4, "scheme": "Scheme 2", "arrow": "25 to 26"})
participant(r, "s25", "m25", "reactant", (4,))
participant(r, "pida", "pida", "reagent", (4,), equiv=("1.2 equiv", {"number": 1.2, "unit": "equiv"}), reference="m25")
participant(r, "tfa", "tfa", "solvent", (4,))
participant(r, "acoh", "acoh", "solvent", (4,))
participant(r, "k2co3", "k2co3", "reagent", (4,))
participant(r, "meoh", "meoh", "solvent", (4,))
participant(r, "p26", "m26", "product", (4,), product_role="phenol product")
r["solvent_conditions"]["composition"] = fact("TFA/AcOH (10/1), then MeOH", "TFA/AcOH then MeOH", (4,))
r["solvent_conditions"]["mixture_ratio"] = fact("10/1", "10/1", (4,))
add_step(r, "oxidative_substitution", (4,), "PIDA in TFA/AcOH (10/1)",
         temperature=("0 °C", {"number": 0, "unit": "°C"}), participants=("s25", "pida", "tfa", "acoh"))
add_step(r, "hydrolysis", (4,), "then K2CO3/MeOH hydrolysis", participants=("k2co3", "meoh", "p26"))
outcome(r, "yield_26", "yield", ["m26"], (4,), "68%", {"number": 68, "unit": "%"}, role="phenol 26")
r["review"]["unresolved"] = ["The page-4 scheme does not state scale, reaction times, hydrolysis temperature, or yield measurement basis."]
index(r, 4)


# Table 3: seven independent experiments.
table3 = [
    (1, 1.3, 1.0, "1.0 GPa, CH2Cl2", "ch2cl2", 25, [1], [17, 67, 3], False),
    (2, 1.3, 1.0, "1.0 GPa, CH2Cl2", "ch2cl2", 25, [4], [17, 62, 13], False),
    (3, 1.3, 1.0, "1.0 GPa, CH2Cl2", "ch2cl2", 25, [16], [17, 42, 40], False),
    (4, 1.3, 1.0, "1.0 GPa, CH2Cl2", "ch2cl2", 25, [64], [17, 5, 71], False),
    (5, 1.3, 1.0, "(i) 1.0 GPa, CH2Cl2; (ii) 0.1 MPa, neat", "ch2cl2", [25, 60], [64, 9], [21, None, 69], True),
    (6, 1.0, 1.5, "0.1 MPa, (CH2Cl)2", "dce", 70, [30], [21, 13, 19], False),
    (7, 1.0, 1.0, "(i) 0.1 MPa, neat; (ii) 0.1 MPa, neat", None, [25, 60], [120, 68], [25, None, 57], True),
]
for entry, eq2, eq4, cond, solvent, temps, times, yields, isolated in table3:
    r = reaction(f"table3_entry_{entry:02d}", (5,), "specific_experiment", ["optimization"],
                 "intermolecular Diels–Alder/carbonyl-ene cascade reaction of (+)-2 and (−)-4",
                 "Diels–Alder/carbonyl-ene cascade", mode="da_ene",
                 location={"page": 5, "table": "Table 3", "entry": entry})
    participant(r, "s2", "m2", "reactant", (5,), equiv=(f"{eq2} equiv", {"number": eq2, "unit": "equiv"}), reference="m4")
    participant(r, "s4", "m4", "reactant", (5,), equiv=(f"{eq4} equiv", {"number": eq4, "unit": "equiv"}), reference="m4")
    participant(r, "p6", "m6", "product", (5,), product_role="Si-anti cycloadduct/byproduct")
    participant(r, "p7", "m7", "product", (5,), product_role="Si-syn cycloadduct/intermediate")
    participant(r, "p11", "m11", "product", (5,), product_role="carbonyl-ene cascade product")
    if solvent:
        participant(r, "solv", solvent, "solvent", (5,))
        r["solvent_conditions"]["composition"] = fact(cond, solvent, (5,))
    temp_list = temps if isinstance(temps, list) else [temps]
    for j, (temp, hours) in enumerate(zip(temp_list, times), 1):
        if entry in (5, 7) and j == 2:
            pressure = ("0.1 MPa", "0.1 MPa")
            operation = "concentrated/kept neat and heated for carbonyl-ene conversion"
        else:
            pressure = (("1.0 GPa", "1.0 GPa") if entry <= 5 else ("0.1 MPa", "0.1 MPa"))
            operation = "intermolecular Diels–Alder/cascade reaction under listed conditions"
        add_step(r, f"reaction_{j}", (5,), operation,
                 temperature=(f"{temp} °C", {"number": temp, "unit": "°C"}),
                 time=(f"{hours} h", {"number": hours, "unit": "h"}), pressure=pressure,
                 participants=tuple(p["participant_id"] for p in r["participants"] if p["stage"] == "reaction"))
    for mid, value, role in zip(("m6", "m7", "m11"), yields, ("chloropestolide B", "chloropestolide C", "chloropupukeanolide D")):
        if value is None:
            continue
        outcome(r, f"yield_{mid}", "yield", [mid], (5,), f"{value}%", {"number": value, "unit": "%"},
                basis="isolated" if isolated else "NMR", yield_type="isolated" if isolated else "NMR", role=role)
    if entry == 6:
        outcome(r, "decomposition_note", "decomposition", ["m2", "m4"], (5,),
                "low yields owing to decomposition of the starting materials", "decomposition observed",
                scope={"quantity": "not reported"})
    r["procedure"]["original_text"] = fact(
        f"Table 3 entry {entry}: ratio 2/4 = {eq2}/{eq4}; {cond}; temperatures {temp_list}; times {times} h.",
        f"Table 3 entry {entry}", (5,))
    r["review"]["unresolved"] = ["Complex reactant/product structures and complete atom mapping are unresolved."]
    index(r, 5)


# Scheme 3a migration and Scheme 3b one-pot synthesis.
r = reaction("scheme3a_migration_11_to_1", (5, 6), "specific_experiment", ["derivatization"],
             "migration of the p-orsellinate group of 11",
             "base-promoted intramolecular acyl migration", mode="migration",
             location={"page": 6, "scheme": "Scheme 3a"})
participant(r, "s11", "m11", "reactant", (5, 6))
participant(r, "base", "kotbu", "reagent", (5, 6), equiv=("1.5 equiv", {"number": 1.5, "unit": "equiv"}), reference="m11")
participant(r, "dmf", "dmf", "solvent", (5, 6))
participant(r, "p1", "m1", "product", (5, 6), product_role="migration product")
add_step(r, "migration", (5, 6), "treatment with KOt-Bu in DMF",
         temperature=("0 °C", {"number": 0, "unit": "°C"}), participants=("s11", "base", "dmf", "p1"))
outcome(r, "yield_1", "yield", ["m1"], (5, 6), "83%", {"number": 83, "unit": "%"}, role="chloropupukeananin")
r["review"]["unresolved"] = ["Reaction time and yield measurement basis are not stated on pages 5-6.",
                              "Complex rearrangement lacks atom mapping; classification is partial."]
index(r, 5, 6)

r = reaction("scheme3b_one_pot_2_4_to_1_6", (1, 5, 6), "specific_experiment", ["scale_up"],
             "one-pot biomimetic synthesis of chloropupukeananin",
             "Diels–Alder/carbonyl-ene cascade followed by acyl migration", mode="da_ene_migration",
             location={"page": 6, "scheme": "Scheme 3b"})
participant(r, "s2", "m2", "reactant", (1, 5, 6), equiv=("1.0 equiv", {"number": 1.0, "unit": "equiv"}), reference="m4")
participant(r, "s4", "m4", "reactant", (1, 5, 6), equiv=("1.0 equiv", {"number": 1.0, "unit": "equiv"}), reference="m4")
participant(r, "base", "kotbu", "reagent", (5, 6))
participant(r, "dmf", "dmf", "solvent", (5, 6))
participant(r, "p1", "m1", "product", (1, 5, 6), product_role="desired one-pot product")
participant(r, "p6", "m6", "product", (5, 6), product_role="cycloadduct byproduct")
add_step(r, "neat_cycloaddition", (5, 6), "neat atmospheric-pressure cascade reaction",
         temperature=("5 °C", {"number": 5, "unit": "°C"}), time=("168 h", {"number": 168, "unit": "h"}),
         pressure=("0.1 MPa (atmospheric pressure)", "0.1 MPa"), participants=("s2", "s4"))
add_step(r, "carbonyl_ene_heating", (5, 6), "then heated neat",
         temperature=("60 °C", {"number": 60, "unit": "°C"}), time=("9 h", {"number": 9, "unit": "h"}),
         pressure=("0.1 MPa (atmospheric pressure)", "0.1 MPa"), participants=("s2", "s4"))
add_step(r, "migration", (5, 6), "KOt-Bu/DMF migration step",
         temperature=("0 °C", {"number": 0, "unit": "°C"}), time=("30 min", {"number": 30, "unit": "min"}),
         participants=("base", "dmf", "p1", "p6"))
outcome(r, "yield_1", "yield", ["m1"], (1, 5, 6), "64%", {"number": 64, "unit": "%"},
        role="chloropupukeananin")
outcome(r, "yield_6", "yield", ["m6"], (5, 6), "20%", {"number": 20, "unit": "%"}, role="chloropestolide B")
r["review"]["unresolved"] = ["KOt-Bu amount and yield measurement basis are not stated in Scheme 3b on pages 5-6.",
                              "Complex cascade/migration atom mapping is unresolved; classification records only transformation-level bond classes."]
index(r, 1, 5, 6)


# Non-synthetic and contextual observations retained separately.
bundle["source_context"] = [
    {"context_id": f"{DOC}:page1_article_metadata", "document_id": DOC, "kind": "article_metadata",
     "evidence_ids": ev(1), "text": "Title, authors, DOI 10.1021/acs.joc.1c02108, graphical abstract, and Figure 1 precursor structures."},
    {"context_id": f"{DOC}:page1_figure1", "document_id": DOC, "kind": "background_figure",
     "evidence_ids": ev(1), "text": "Figure 1 identifies chloropupukeananin and proposed biosynthetic precursors; it is not a current-paper experiment."},
    {"context_id": f"{DOC}:page2_figure3", "document_id": DOC, "kind": "mechanistic_orientation_model",
     "evidence_ids": ev(2), "text": "Figure 3 depicts four Re/Si and syn/anti approach orientations; no standalone experiment is reported."},
    {"context_id": f"{DOC}:page3_natural_isolation_ratio", "document_id": DOC, "kind": "cited_literature_comparison",
     "evidence_ids": ev(3), "text": "The 'isolation (ref 12e)' ratio is cited natural-product data, not a Table 1 experiment of this paper."},
    {"context_id": f"{DOC}:page5_general_instrumentation", "document_id": DOC, "kind": "general_information",
     "evidence_ids": ev(5), "text": "General IR/HRMS/NMR instrumentation begins; these measurements are not synthetic experiments."},
    {"context_id": f"{DOC}:page6_characterization_14_15_13_16", "document_id": DOC, "kind": "characterization",
     "evidence_ids": ev(6), "text": "Independent NMR/IR/HRMS/optical rotation data for 14, 15, and the 13/16 mixture are characterization, not new syntheses."},
    {"context_id": f"{DOC}:page6_general_reagents", "document_id": DOC, "kind": "general_information",
     "evidence_ids": ev(6), "text": "General atmosphere, apparatus, chromatography, and commercial-reagent statements apply only where explicitly inherited; they are not separate experiments."},
]


coverage_regions = {
    1: [
        {"region": "title, abstract, and graphical abstract", "finding": "article metadata and a summary of the one-pot synthesis; the specific one-pot experiment is linked to Scheme 3b rather than duplicated", "evidence_ids": ev(1)},
        {"region": "Figure 1 and introductory text", "finding": "background structures and proposed precursor relation; no current-paper experimental row", "evidence_ids": ev(1)},
    ],
    2: [
        {"region": "Figure 2", "finding": "proposed branched biosynthetic mechanism, recorded as proposed_mechanism", "evidence_ids": ev(2)},
        {"region": "Figure 3", "finding": "stereochemical approach-orientation model, retained as source context", "evidence_ids": ev(2)},
        {"region": "Results and Discussion opening", "finding": "description of Scheme 1 as prior work, not a current experiment", "evidence_ids": ev(2)},
    ],
    3: [
        {"region": "Scheme 1", "finding": "previously published Diels–Alder/carbonyl-ene example, classified as literature_example", "evidence_ids": ev(3)},
        {"region": "Table 1 rows 1-14", "finding": "14 independent current solvent/pressure experiments extracted; cited isolation row excluded from experiment count", "evidence_ids": ev(3)},
        {"region": "Table 1 footnotes", "finding": "1H NMR ratio basis, isolated yields for entry 6, and 1.1 equiv/1.0 GPa conditions for entries 13-14 captured", "evidence_ids": ev(3)},
    ],
    4: [
        {"region": "Table 2 rows 1-6 and footnotes", "finding": "one previous-work row plus five current optimization experiments; chlorohydrin isomer mixture and NMR basis retained", "evidence_ids": ev(4)},
        {"region": "Table 2 reaction arrow", "finding": "row 6 linked to the two-step deprotection outcome rather than duplicated", "evidence_ids": ev(4)},
        {"region": "Scheme 2 main arrows", "finding": "24→25 two-step sequence and preparative 25→26 one-pot sequence extracted separately", "evidence_ids": ev(4)},
        {"region": "Scheme 2 inset screen", "finding": "four solvent-ratio rows extracted independently; 27/28 kept as generic regioisomer pools", "evidence_ids": ev(4)},
    ],
    5: [
        {"region": "Table 3 rows 1-7 and footnotes", "finding": "seven independent Diels–Alder/carbonyl-ene experiments extracted with staged temperature/time and NMR versus isolated yield basis", "evidence_ids": ev(5)},
        {"region": "migration and one-pot discussion", "finding": "Scheme 3a and 3b outcomes/conditions linked to the page-6 arrows", "evidence_ids": ev(5)},
        {"region": "general instrumentation opening", "finding": "spectroscopy setup retained as source context, not a synthesis", "evidence_ids": ev(5)},
    ],
    6: [
        {"region": "Scheme 3a-b", "finding": "migration and one-pot synthesis arrows extracted; graphical abstract duplicate merged into the same one-pot record", "evidence_ids": ev(6)},
        {"region": "General Procedure for Table 1 entries 2-12", "finding": "shared charges, solvent volume, room-temperature/48 h conditions and entry-6 purification inherited into the corresponding rows", "evidence_ids": ev(6)},
        {"region": "characterization of 14, 15, and 13/16 mixture", "finding": "spectral/optical-rotation data retained in source_context and not made into reactions", "evidence_ids": ev(6)},
        {"region": "(+)-iso-A82775C (Table 2, entry 6) procedure", "finding": "detailed charges, staged organoindium preparation/coupling/desilylation/workup and two-step 92% outcome merged with Table 2 entry 6", "evidence_ids": ev(6)},
    ],
}

coverage = {"document_id": DOC, "semantic_accuracy_verified": False, "pages": []}
for page in range(1, 7):
    coverage["pages"].append({
        "page": page,
        "state": "attempted",
        "regions": coverage_regions[page],
        "record_ids": page_records[page],
        "unresolved": [
            "Complex source structures are not converted to source-checked atom-mapped/isomeric SMILES.",
            "Only PDF pages 1-6 were claimed as visually covered; unavailable SI was not associated.",
        ],
    })

raw_inventory = {
    "document_id": DOC,
    "bounded_pages": [1, 2, 3, 4, 5, 6],
    "visual_inspection": "all six supplied page JPGs inspected at original detail",
    "reaction_ids": [r["reaction_id"] for r in bundle["reactions"]],
    "molecule_ids": [m["molecule_id"] for m in bundle["molecules"]],
    "source_context_ids": [c["context_id"] for c in bundle["source_context"]],
    "limitations": [
        "Machine extraction remains needs_review; schema validation is not chemistry verification.",
        "No SI file was explicitly provided or associated.",
        "Full atomwise stereochemical and bond mapping of the complex natural products is unresolved.",
    ],
}

OUT.parent.mkdir(parents=True, exist_ok=True)
RAW.parent.mkdir(parents=True, exist_ok=True)
COVERAGE.parent.mkdir(parents=True, exist_ok=True)
st.write(OUT, bundle)
st.write(RAW, raw_inventory)
st.write(COVERAGE, coverage)
print(json.dumps({"molecules": len(bundle["molecules"]), "reactions": len(bundle["reactions"]),
                  "source_context": len(bundle["source_context"]), "coverage_pages": len(coverage["pages"])}, ensure_ascii=False))
