#!/usr/bin/env python3
"""Build the evidence-linked v3.1 submission for doc_8a0a9af017c80023.

This is a source-transcription helper.  It deliberately leaves compact-drawing
structures unresolved and never upgrades machine records to source-reviewed.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

DOC = "doc_8a0a9af017c80023"
RUN = ROOT / "runs/current"
DATASET = json.loads((RUN / "dataset.json").read_text())
EVT = {p: f"{DOC}:page{p}:text" for p in range(1, 6)}
EVI = {p: f"{DOC}:page{p}:image" for p in range(1, 6)}


def f(raw=None, normalized=None, pages=(), method="original_extraction", derivation=None,
      status=None):
    ev = []
    for page in pages:
        ev.extend([EVI[page], EVT[page]])
    return st.fact(raw, normalized, ev, method=method, derivation=derivation, status=status)


def fi(normalized, pages, derivation):
    return f(None, normalized, pages, method="inference", derivation=derivation)


def mol(mid, label, pages, *, name=None, chirality="unresolved", geometry="unresolved",
        representation="not_encoded", assignment=None, double_bonds=None, mixture=None,
        note=None, smiles=None):
    ev = []
    for page in pages:
        ev.extend([EVI[page], EVT[page]])
    m = st.new_molecule(mid, DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f(label, label, pages)
    m["names"]["original"] = f(name or label, name or label, pages)
    m["names"]["english"] = f(name, name, pages) if name else st.fact()
    m["names"]["abbreviations"] = st.fact()
    m["names"]["aliases"] = st.fact()
    m["cas_numbers"] = st.fact()
    m["generic_template"] = None
    m["stereochemistry"]["atom_centers"] = (
        f("stereochemical wedge/dash shown in source drawing",
          "source drawing contains central stereochemical information", pages)
        if chirality in {"specified", "partially_specified"} else st.fact()
    )
    m["stereochemistry"]["double_bonds"] = (
        f(double_bonds, double_bonds, pages) if double_bonds is not None else st.fact()
    )
    m["stereochemistry"]["axial"] = st.fact()
    m["stereochemistry"]["mixture"] = (
        f(mixture, mixture, pages) if mixture is not None else st.fact()
    )
    m["stereochemistry"]["assignment"] = (
        f(assignment, assignment, pages) if assignment is not None else st.fact()
    )
    m["stereochemistry"]["chirality_status"] = f(chirality, chirality, pages)
    m["stereochemistry"]["geometry_status"] = f(geometry, geometry, pages)
    m["stereochemistry"]["representation_status"] = f(representation, representation, pages)
    m["stereochemistry"]["isomeric_smiles"] = st.fact()
    m["stereochemistry"]["identity_key"] = st.fact()
    m["identity_note"] = f(note or "Compact main-text identity retained; no SI was supplied.",
                             note or "Compact main-text identity retained; no SI was supplied.", pages)
    m["external_properties"] = []
    if smiles:
        diag = st.stereo_inspect(smiles)
        m["representation_type"] = "smiles"
        m["smiles"] = f(label, smiles, pages, method="structure_derivation",
                          derivation="Candidate SMILES derived from the legible source structure; graph parsing is not source verification.")
        m["stereochemistry"]["chirality_status"] = fi(
            "achiral", pages, "The complete candidate graph has no stereogenic center/axis; checked with standard stereo.")
        m["stereochemistry"]["geometry_status"] = fi(
            "not_applicable", pages, "The candidate graph has no configurable alkene geometry.")
        m["stereochemistry"]["representation_status"] = fi(
            "encoded", pages, "All stereochemical information applicable to the candidate graph is encoded.")
        m["stereochemistry"]["isomeric_smiles"] = f(
            None, diag["canonical_isomeric_smiles"], pages, method="structure_derivation",
            derivation="Canonical isomeric SMILES returned by standard stereo for the source-derived candidate graph.")
        m["stereochemistry"]["identity_key"] = f(
            None, diag["stereo_graph_key"], pages, method="structure_derivation",
            derivation="stereo_graph_key returned by standard stereo; not a sample identity.")
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = st.fact()
    return m


def template_mol(mid, label, pages, fixed, variables, note):
    m = mol(mid, label, pages, name=label, chirality="partially_specified",
            geometry="unresolved", representation="not_encoded", note=note)
    m["representation_type"] = "generic_template"
    m["generic_template"] = {
        "fixed_structure": f(fixed, fixed, pages),
        "variables": [
            {
                "label": v,
                "attachment_points": f(f"source-marked attachment for {v}", f"source-marked attachment for {v}", pages),
                "definition": f(desc, desc, pages),
                "allowed_values": st.fact(),
                "constraints": st.fact(),
            }
            for v, desc in variables
        ],
        "instance_mapping": [],
        "reason": f(note, note, pages),
    }
    return m


MOLECULES = {}


def addm(m):
    MOLECULES[m["molecule_id"]] = m
    return m["molecule_id"]


def ensure(mid, label, pages, **kwargs):
    if mid not in MOLECULES:
        addm(mol(mid, label, pages, **kwargs))
    return mid


def charge(p, key, raw, normalized, pages, *, method="original_extraction", derivation=None):
    p["charges"][key] = f(raw, normalized, pages, method=method, derivation=derivation)


def participant(pid, mid, role, pages, *, stage="reaction", product_role=None):
    p = st.new_participant(pid, mid, role, [EVI[x] for x in pages], stage=stage)
    if product_role:
        p["product_role"] = f(product_role, product_role, pages)
    return p


def classification(kind, pages):
    ev = [EVI[p] for p in pages] + [EVT[p] for p in pages]
    c = {
        "standard_id": st.STANDARD,
        "mapped_reaction_smiles": st.fact(),
        "bond_changes": [],
        "labels": None,
        "coverage": None,
        "limitations": [],
        "model_votes": [],
    }
    changes = []
    if kind == "hydroalkylation_alkyne":
        changes = [
            {"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
            {"change_id": "b2", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "H"}], "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain", "hydrogen_count": 1},
            {"change_id": "b3", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 3, "after_order": 2, "change_type": "order_decreased"},
        ]
    elif kind == "hydroalkylation_allene":
        changes = [
            {"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
            {"change_id": "b2", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "H"}], "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain", "hydrogen_count": 1},
            {"change_id": "b3", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 2, "after_order": 1, "change_type": "order_decreased"},
        ]
    elif kind == "methylation":
        changes = [{"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "N"}], "before_order": 0, "after_order": 1, "change_type": "formed"}]
    elif kind == "esterification":
        changes = [{"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "O"}], "before_order": 0, "after_order": 1, "change_type": "formed"}]
    elif kind in {"amidation", "cyclization"}:
        changes = [{"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "N"}], "before_order": 0, "after_order": 1, "change_type": "formed"}]
    elif kind == "deuterium_exchange":
        changes = [{"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "H"}], "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain", "hydrogen_count": 1}]
    elif kind in {"not_applicable", "unresolved"}:
        coverage = kind
        c["labels"] = fi([], pages, "No observed bond-forming transformation is asserted for this record.")
        c["coverage"] = fi(coverage, pages, "Classification is not applicable to a proposed mechanism/no-reaction record." if kind == "not_applicable" else "Transient species or incomplete structures prevent a defensible bond inventory.")
        c["limitations"] = ["No atom mapping or complete structure set is available from the supplied main-text pages."]
        return c
    labels = {}
    for ch in changes:
        if ch["change_type"] in {"formed", "hydrogen_gain", "order_increased"}:
            lab = st.bond_label(ch["atoms"][0]["element"], ch["atoms"][1]["element"], ch["after_order"])
            labels.setdefault(lab["code"], {**lab, "supporting_change_ids": []})["supporting_change_ids"].append(ch["change_id"])
    deriv = "Partial net bond inventory inferred from the source transformation drawing; atom identities are unresolved and no complete atom map is claimed."
    c["bond_changes"] = [f(None, ch, pages, method="inference", derivation=deriv) for ch in changes]
    c["labels"] = f(None, list(labels.values()), pages, method="inference", derivation=deriv)
    c["coverage"] = fi("partial", pages, deriv)
    c["limitations"] = [
        "Compact main-text drawings and absent SI do not support a complete atom map.",
        "C-H includes the net implicit-hydrogen gain; the proton/deuteron source is not atom-mapped.",
    ]
    return c


def reaction(rid, pages, *, category, name, transform, class_kind, record_type="specific_experiment",
             general=()):
    ev = [EVI[p] for p in pages] + [EVT[p] for p in pages]
    r = st.new_reaction(f"{DOC}:{rid}", DOC, ev)
    r["record_type"] = record_type
    r["experiment_categories"] = f(category, category if isinstance(category, list) else [category], pages)
    r["location"] = f(f"PDF page(s) {', '.join(map(str, pages))}", {"pages": list(pages)}, pages)
    r["general_scheme_ids"] = list(general)
    r["reaction_name"] = f(name, name, pages)
    r["transformation_type"] = fi(transform, pages, "Normalized transformation description inferred from the source reaction name/drawing.")
    r["classification"] = classification(class_kind, pages)
    r["review"] = {
        "status": "needs_review",
        "method": "model",
        "reviewer": None,
        "checked_fields": [],
        "unresolved": ["No explicitly associated SI was supplied; machine record remains needs_review."],
    }
    r["completeness"] = "partial"
    return r


def step(r, sid, pages, pids, *, operation, temperature=None, time=None, atmosphere=None,
         pressure=None, stage="reaction"):
    s = st.new_step(sid, [EVI[p] for p in pages], stage=stage)
    s["participant_ids"] = list(pids)
    s["operation"] = f(operation, operation, pages)
    s["temperature"] = f(temperature[0], temperature[1], pages) if temperature else st.fact()
    s["time"] = f(time[0], time[1], pages) if time else st.fact()
    s["atmosphere"] = f(atmosphere, atmosphere, pages) if atmosphere else st.fact()
    s["pressure"] = f(pressure, pressure, pages) if pressure else st.fact()
    r["steps"].append(s)


def outcome(r, oid, kind, targets, pages, raw, normalized, *, product_role=None,
            yield_type=None, basis=None, scope=None):
    o = st.new_outcome(oid, kind, targets, [EVI[p] for p in pages])
    o["value"] = f(raw, normalized, pages)
    if product_role:
        o["product_role"] = f(product_role, product_role, pages)
    if yield_type:
        o["yield_type"] = f(yield_type, yield_type, pages)
    if basis:
        o["measurement_basis"] = f(basis, basis, pages)
    if scope is not None:
        o["scope"] = f(scope, scope, pages)
    r["outcomes"].append(o)


def set_solvent(r, pages, composition, *, volume=None, concentration=None, basis=None):
    r["solvent_conditions"]["composition"] = f(composition, composition, pages)
    if volume:
        r["solvent_conditions"]["total_volume"] = f(volume[0], volume[1], pages)
    if concentration:
        r["solvent_conditions"]["concentration"] = f(concentration[0], concentration[1], pages)
    if basis:
        r["solvent_conditions"]["concentration_basis"] = f(basis, basis, pages)


def note(r, pages, text):
    r["notes"].append(f(text, text, pages))


REACTIONS = []

# Generic Scheme 1(b), with R groups retained as templates.
addm(template_mol("m_generic_imidazolidinone_1", "generic imidazolidinone derivative 1", [1],
                  "imidazolidinone core with source-defined R substituent and drawn stereochemistry",
                  [("R", "amino-acid-derived substituent")],
                  "Scheme 1(b) is generic and does not define a single concrete R group."))
addm(template_mol("m_generic_unsaturated_2_4", "generic aryl alkyne or allene", [1],
                  "Ar-substituted alkyne or allene coupling partner",
                  [("Ar", "aryl substituent")],
                  "Scheme 1(b) shows alternative alkyne/allene classes, not a concrete molecule."))
addm(template_mol("m_generic_product_3", "generic α-alkenyl imidazolidinone product", [1],
                  "quaternary α-alkenyl imidazolidinone core with E-alkene drawn",
                  [("R", "amino-acid-derived substituent"), ("Ar", "aryl substituent")],
                  "Scheme 1(b) is a generic transformation."))
g = reaction("scheme1b_generic", [1], category="other", name="intermolecular hydroalkylation of alkynes or allenes",
             transform="generic α-alkenylation", class_kind="hydroalkylation_allene", record_type="generic_scheme")
g["participants"] = [participant("generic_1", "m_generic_imidazolidinone_1", "reactant", [1]),
                     participant("generic_2_4", "m_generic_unsaturated_2_4", "reactant", [1]),
                     participant("generic_3", "m_generic_product_3", "product", [1], product_role="desired product")]
note(g, [1], "Scheme 1(a) is cited prior work (Clayden et al.) and is retained in source_context, not counted as this-work experiment.")
REACTIONS.append(g)

# Core molecules and common reagents.
ensure("m_1a", "1a", [1, 2, 3, 4], name="imidazolidinone 1a", chirality="partially_specified",
       geometry="not_applicable", assignment="source drawing contains fixed relative stereochemistry")
ensure("m_2a", "2a", [1, 2, 4], name="1-phenyl-1-propyne (source label 2a)", smiles="CC#Cc1ccccc1")
ensure("m_allene_4", "4 / 4a (allene substrate)", [3], name="1-phenylallene",
       smiles="C=C=Cc1ccccc1", note="This allene substrate is distinct from amino-acid product 4a despite the caption/prose label overlap.")
MOLECULES["m_allene_4"]["stereochemistry"]["geometry_status"] = fi(
    "unresolved", [3], "standard stereo reports an unspecified potential double-bond element for this allene graph; the main text gives no E/Z assignment for the allene substrate.")
MOLECULES["m_allene_4"]["stereochemistry"]["representation_status"] = fi(
    "partial", [3], "The candidate connectivity is encoded, but toolkit allene/double-bond stereo support does not justify a complete encoding claim.")
for mid, label, pages, name in [
    ("m_LiHMDS", "LiHMDS", [2, 3, 4], "lithium bis(trimethylsilyl)amide"),
    ("m_NaHMDS", "NaHMDS", [2], "sodium bis(trimethylsilyl)amide"),
    ("m_KHMDS", "KHMDS", [2], "potassium bis(trimethylsilyl)amide"),
    ("m_THF", "THF", [2, 3, 4], "tetrahydrofuran"),
    ("m_toluene", "toluene", [2], "toluene"),
]:
    ensure(mid, label, pages, name=name, chirality="unresolved", geometry="not_applicable")

# Product registry 3a-ai from Scheme 2.
scope = [
    ("3a", "R = H", 83, ">99:1", ">99:1"), ("3b", "R = 3-F", 85, ">99:1", ">99:1"),
    ("3c", "R = 3-Me", 80, ">99:1", ">99:1"), ("3d", "R = 3-OMe", 77, ">99:1", ">99:1"),
    ("3e", "R = 4-Br", 85, "92:8", ">99:1"), ("3f", "R = 4-Cl", 89, ">99:1", ">99:1"),
    ("3g", "R = 3-CF3", 84, ">99:1", "96:4"), ("3h", "naphthyl product", 90, ">99:1", None),
    ("3i", "naphthyl product", 83, ">99:1", None), ("3j", "thiophenyl product", 78, "96:4", ">99:1"),
    ("3k", "thiophenyl product", 82, ">99:1", "99:1"), ("3l", "bis-alkyne-derived monoaddition product", 78, "95:5", None),
    ("3m", "polyfluoro product", 76, ">99:1", "97:3"), ("3n", "conjugated enyne-derived product", 82, "97:3", None),
    ("3o", "R = 3-F", 65, "97:3", "97:3"), ("3p", "R = 4-Cl", 68, "98:2", "97:3"),
    ("3q", "R = H", 84, "98:2", None), ("3r", "R = 3-F", 87, ">99:1", ">99:1"),
    ("3s", "R = 4-Cl", 90, ">99:1", ">99:1"), ("3t", "R = H", 87, ">99:1", "95:5"),
    ("3u", "R = 4-Cl", 88, ">99:1", "95:5"), ("3v", "R = H", 73, ">99:1", ">99:1"),
    ("3w", "R = 4-Cl", 76, ">99:1", ">99:1"), ("3x", "R = 3-F", 78, ">99:1", ">99:1"),
    ("3y", "R = H", 76, ">99:1", ">99:1"), ("3z", "R = 3-F", 78, ">99:1", ">99:1"),
    ("3aa", "R = 4-C≡CCH3", 74, ">99:1", ">99:1"), ("3ab", "methoxy/naphthyl product", 75, ">99:1", None),
    ("3ac", "R = H", 85, ">99:1", None), ("3ad", "R = 4-Me", 80, ">99:1", "86:14"),
    ("3ae", "R = H", 88, ">99:1", None), ("3af", "R = H", 78, ">99:1", None),
    ("3ag", "R = H", 77, ">99:1", "87:13"), ("3ah", "R = H", 82, ">99:1", None),
    ("3ai", "R = H", 80, ">99:1", None),
]
for label, desc, yld, dr, er in scope:
    ensure(f"m_{label}", label, [2, 3] if label in {"3a", "3s", "3q", "3v", "3y", "3ac"} else [2],
           name=f"source product {label} ({desc})", chirality="partially_specified", geometry="specified",
           representation="not_encoded", double_bonds="E alkene drawn/reported",
           assignment=("X-ray confirmed relative and absolute configuration" if label == "3ac" else "source drawing and selectivity reported; absolute R/S label not transcribed"))
    if label != "3a":
        ensure(f"m_1_for_{label}", f"imidazolidinone precursor shown for {label}", [2],
               chirality="partially_specified", geometry="not_applicable",
               note="Specific compact drawing recognized, but not expanded to an unverified SMILES.")
        ensure(f"m_2_for_{label}", f"alkyne precursor shown for {label}", [2],
               chirality="unresolved", geometry="not_applicable",
               note="Specific alkyne identity is retained by source/product label; no SI-quality SMILES asserted.")

# Table 1 optimization, one record per row.
table_rows = [
    (1, "m_KHMDS", 2.0, "80 °C", {"number": 80, "unit": "°C"}, 73, 79, ">99:1", "68:32", "THF", 20),
    (2, "m_NaHMDS", 2.0, "80 °C", {"number": 80, "unit": "°C"}, 75, 82, ">99:1", "56:44", "THF", 20),
    (3, "m_LiHMDS", 2.0, "80 °C", {"number": 80, "unit": "°C"}, 83, 91, ">99:1", ">99:1", "THF", 20),
    (4, "m_LiHMDS", 2.0, "80 °C", {"number": 80, "unit": "°C"}, 78, None, ">99:1", "94:6", "toluene", 20),
    (5, "m_LiHMDS", 2.0, "50 °C", {"number": 50, "unit": "°C"}, 48, None, ">99:1", ">99:1", "THF", 20),
    (6, "m_LiHMDS", 2.0, "r.t.", "r.t.", 25, None, ">99:1", ">99:1", "THF", 48),
    (7, "m_LiHMDS", 1.0, "80 °C", {"number": 80, "unit": "°C"}, 23, None, ">99:1", ">99:1", "THF", 20),
    (8, "m_LiHMDS", 0.5, "80 °C", {"number": 80, "unit": "°C"}, "<5", None, ">99:1", "ND", "THF", 20),
]
for entry, base_mid, eq, traw, tnorm, yiso, ynmr, dr, er, solvent, hours in table_rows:
    r = reaction(f"table1_entry{entry}", [2], category="optimization", name="optimization of alkyne hydroalkylation",
                 transform="α-alkenylation by alkyne hydroalkylation", class_kind="hydroalkylation_alkyne",
                 general=(f"{DOC}:scheme1b_generic",))
    ps = [participant("p_1a", "m_1a", "reactant", [2]), participant("p_2a", "m_2a", "reactant", [2]),
          participant("p_base", base_mid, "reagent", [2]),
          participant("p_solvent", "m_toluene" if solvent == "toluene" else "m_THF", "solvent", [2]),
          participant("p_3a", "m_3a", "product", [2], product_role="desired product")]
    charge(ps[0], "amount", "0.2 mmol", {"number": 0.2, "unit": "mmol"}, [2])
    charge(ps[1], "amount", "0.2 mmol", {"number": 0.2, "unit": "mmol"}, [2])
    charge(ps[2], "volume", "0.4 mL of 1.0 M HMDS solution", {"number": 0.4, "unit": "mL"}, [2])
    charge(ps[2], "equivalents", f"{eq} eq.", {"number": eq, "unit": "equiv"}, [2])
    ps[2]["equivalence_reference"] = f("relative to 1a", "m_1a", [2])
    r["participants"] = ps
    r["reaction_scale"] = f("1a (0.2 mmol)", {"reference_molecule_id": "m_1a", "amount": {"number": 0.2, "unit": "mmol"}}, [2])
    composition = "toluene as reaction solvent; HMDS reagent supplied in THF" if solvent == "toluene" else "THF"
    set_solvent(r, [2], composition, volume=("2.0 mL", {"number": 2.0, "unit": "mL"}),
                concentration=("0.08 M", {"number": 0.08, "unit": "M"}), basis="reported reaction concentration")
    step(r, "reaction", [2], [p["participant_id"] for p in ps[:-1]], operation="combine under Table 1 conditions",
         temperature=(traw, tnorm), time=(f"{hours} h", {"number": hours, "unit": "h"}))
    yval = ({"number": 5, "unit": "%", "qualifier": "lt", "raw_text": "<5%"} if yiso == "<5" else {"number": yiso, "unit": "%"})
    outcome(r, "isolated_yield", "yield", ["m_3a"], [2], f"{yiso}%", yval,
            product_role="desired product", yield_type="isolated yield", basis="isolated product")
    if ynmr is not None:
        outcome(r, "nmr_yield", "yield", ["m_3a"], [2], f"({ynmr})", {"number": ynmr, "unit": "%"},
                product_role="desired product", yield_type="NMR yield", basis="1H NMR with 1,3,5-trimethoxybenzene (TMB) as reference")
    outcome(r, "dr", "selectivity", ["m_3a"], [2], dr, dr, basis="1H NMR", scope={"metric": "dr", "isomer_order": "major:minor"})
    outcome(r, "er", "selectivity", ["m_3a"], [2], er, er, basis="chiral HPLC", scope={"metric": "er", "isomer_order": "major:minor"})
    if entry == 4:
        note(r, [2], "Table footnote f: toluene was used as the solvent; the HMDS solution carrier is still reported as THF.")
    REACTIONS.append(r)

# Scheme 2 scope, one record per displayed product.
for label, desc, yld, dr, er in scope:
    r = reaction(f"scheme2_{label}", [2], category="substrate_scope", name="hydroalkylation of an alkyne with an imidazolidinone derivative",
                 transform="α-alkenylation by alkyne hydroalkylation", class_kind="hydroalkylation_alkyne",
                 general=(f"{DOC}:scheme1b_generic",))
    mid1 = "m_1a" if label == "3a" else f"m_1_for_{label}"
    mid2 = "m_2a" if label == "3a" else f"m_2_for_{label}"
    ps = [participant("p_imidazolidinone", mid1, "reactant", [2]), participant("p_alkyne", mid2, "reactant", [2]),
          participant("p_LiHMDS", "m_LiHMDS", "reagent", [2]), participant("p_THF", "m_THF", "solvent", [2]),
          participant("p_product", f"m_{label}", "product", [2], product_role="desired product")]
    charge(ps[0], "amount", "0.2 mmol", {"number": 0.2, "unit": "mmol"}, [2])
    charge(ps[1], "amount", "0.2 mmol", {"number": 0.2, "unit": "mmol"}, [2])
    charge(ps[2], "volume", "0.4 mL, 1.0 M in THF", {"number": 0.4, "unit": "mL"}, [2])
    charge(ps[2], "mol_percent", "200 mol%", {"number": 200, "unit": "mol%"}, [2])
    ps[2]["equivalence_reference"] = f("relative to imidazolidinone substrate", mid1, [2])
    r["participants"] = ps
    r["reaction_scale"] = f("imidazolidinone (0.2 mmol)", {"reference_molecule_id": mid1, "amount": {"number": 0.2, "unit": "mmol"}}, [2])
    set_solvent(r, [2], "THF", volume=("2.0 mL", {"number": 2.0, "unit": "mL"}),
                concentration=("0.08 M", {"number": 0.08, "unit": "M"}), basis="reported reaction concentration")
    step(r, "reaction", [2], [p["participant_id"] for p in ps[:-1]], operation="Scheme 2 hydroalkylation",
         temperature=("80 °C", {"number": 80, "unit": "°C"}), time=("20 h", {"number": 20, "unit": "h"}))
    outcome(r, "yield", "yield", [f"m_{label}"], [2], f"{yld}%", {"number": yld, "unit": "%"},
            product_role="desired product", yield_type="isolated yield", basis="isolated product")
    outcome(r, "dr", "selectivity", [f"m_{label}"], [2], dr, dr, basis="1H NMR", scope={"metric": "dr", "isomer_order": "major:minor"})
    if er is not None:
        outcome(r, "er", "selectivity", [f"m_{label}"], [2], er, er,
                scope={"metric": "er", "isomer_order": "major:minor"})
    else:
        note(r, [2], "Scheme 2 footnote a: er not determined.")
        r["review"]["unresolved"].append("er was explicitly not determined for this product.")
    REACTIONS.append(r)

# Scheme 3(a): starting-material preparation and stereodivergent hydroalkylation.
for mid, label, name in [
    ("m_LPhe", "L-Phe", "L-phenylalanine"), ("m_LPheOMeHCl", "L-Phe-OMe·HCl", "L-phenylalanine methyl ester hydrochloride"),
    ("m_LPheNHMe", "L-Phe-NHMe", "L-phenylalanine N-methylamide"), ("m_1h_prime", "1′h", "cis/trans imidazolidinone diastereomer 1′h"),
    ("m_1i_prime", "1′i", "cis/trans imidazolidinone diastereomer 1′i"), ("m_1h", "1h", "N-methyl imidazolidinone 1h"),
    ("m_1i", "1i", "N-methyl imidazolidinone 1i"), ("m_3s_prime", "3s′", "stereodivergent product 3s′"),
    ("m_SOCl2", "SOCl2", "thionyl chloride"), ("m_MeOH", "MeOH", "methanol"), ("m_MeNH2", "MeNH2", "methylamine"),
    ("m_EtOH", "EtOH", "ethanol"), ("m_pivaldehyde", "t-BuCHO", "pivaldehyde"), ("m_YbOTf3", "Yb(OTf)3", "ytterbium triflate"),
    ("m_CHCl3", "CHCl3", "chloroform"), ("m_MeI", "MeI", "methyl iodide"), ("m_K2CO3", "K2CO3", "potassium carbonate"),
    ("m_DMF", "DMF", "N,N-dimethylformamide"),
]:
    kw = {}
    if mid in {"m_LPhe", "m_LPheOMeHCl", "m_LPheNHMe", "m_1h_prime", "m_1i_prime", "m_1h", "m_1i", "m_3s_prime"}:
        kw = {"chirality": "partially_specified", "geometry": "not_applicable" if mid != "m_3s_prime" else "specified",
              "double_bonds": "E alkene drawn/reported" if mid == "m_3s_prime" else None}
    ensure(mid, label, [3], name=name, **kw)

def simple_arrow(rid, input_mid, output_mids, reagents, solvent_mid, temp, hrs, yields, class_kind, name):
    r = reaction(rid, [3], category="other", name=name, transform=name, class_kind=class_kind)
    ps = [participant("p_input", input_mid, "reactant", [3])]
    for i, mid in enumerate(reagents): ps.append(participant(f"p_reagent_{i+1}", mid, "reagent", [3]))
    if solvent_mid: ps.append(participant("p_solvent", solvent_mid, "solvent", [3]))
    for i, mid in enumerate(output_mids): ps.append(participant(f"p_product_{i+1}", mid, "product", [3], product_role="product"))
    r["participants"] = ps
    set_solvent(r, [3], MOLECULES[solvent_mid]["names"]["original"]["normalized_value"] if solvent_mid else "unreported")
    step(r, "reaction", [3], [p["participant_id"] for p in ps if p["role"] != "product"], operation=name,
         temperature=temp, time=hrs)
    for i, (mid, y) in enumerate(zip(output_mids, yields), 1):
        if y is not None: outcome(r, f"yield_{i}", "yield", [mid], [3], f"{y}%", {"number": y, "unit": "%"}, product_role="product")
    return r

REACTIONS.append(simple_arrow("scheme3a_esterification", "m_LPhe", ["m_LPheOMeHCl"], ["m_SOCl2"], "m_MeOH",
                              ("r.t.", "r.t."), ("24 h", {"number": 24, "unit": "h"}), [None], "esterification", "formation of L-Phe-OMe·HCl"))
REACTIONS.append(simple_arrow("scheme3a_amidation", "m_LPheOMeHCl", ["m_LPheNHMe"], ["m_MeNH2"], "m_EtOH",
                              ("r.t.", "r.t."), ("48 h", {"number": 48, "unit": "h"}), [None], "amidation", "formation of L-Phe-NHMe"))
cyc = simple_arrow("scheme3a_cyclization", "m_LPheNHMe", ["m_1h_prime", "m_1i_prime"], ["m_pivaldehyde", "m_YbOTf3"], "m_CHCl3",
                   ("r.t.", "r.t."), ("8 h", {"number": 8, "unit": "h"}), [58, 11], "cyclization", "diastereomer-forming imidazolidinone cyclization")
step(cyc, "purification", [3], ["p_product_1", "p_product_2"], operation="compounds 1′h and 1′i were separated using column chromatography", stage="purification")
REACTIONS.append(cyc)
REACTIONS.append(simple_arrow("scheme3a_methylation_1h", "m_1h_prime", ["m_1h"], ["m_MeI", "m_K2CO3"], "m_DMF",
                              ("60 °C", {"number": 60, "unit": "°C"}), ("12 h", {"number": 12, "unit": "h"}), [69], "methylation", "N-methylation of 1′h"))
REACTIONS.append(simple_arrow("scheme3a_methylation_1i", "m_1i_prime", ["m_1i"], ["m_MeI", "m_K2CO3"], "m_DMF",
                              ("60 °C", {"number": 60, "unit": "°C"}), ("12 h", {"number": 12, "unit": "h"}), [70], "methylation", "N-methylation of 1′i"))

for rid, imid, prod, y in [("scheme3a_hydroalkylation_1h", "m_1h", "m_3s", 90), ("scheme3a_hydroalkylation_1i", "m_1i", "m_3s_prime", 89)]:
    r = reaction(rid, [3], category="substrate_scope", name="stereodivergent hydroalkylation", transform="α-alkenylation by alkyne hydroalkylation", class_kind="hydroalkylation_alkyne")
    ps = [participant("p_imid", imid, "reactant", [3]), participant("p_alkyne", "m_2_for_3s", "reactant", [3]),
          participant("p_LiHMDS", "m_LiHMDS", "reagent", [3]), participant("p_THF", "m_THF", "solvent", [3]),
          participant("p_product", prod, "product", [3], product_role="desired product")]
    charge(ps[2], "mol_percent", "200 mol%", {"number": 200, "unit": "mol%"}, [3]); ps[2]["equivalence_reference"] = f("relative to imidazolidinone", imid, [3])
    r["participants"] = ps; set_solvent(r, [3], "THF")
    step(r, "reaction", [3], [p["participant_id"] for p in ps[:-1]], operation="Scheme 3(a) hydroalkylation",
         temperature=("80 °C", {"number": 80, "unit": "°C"}), time=("20 h", {"number": 20, "unit": "h"}))
    outcome(r, "yield", "yield", [prod], [3], f"{y}%", {"number": y, "unit": "%"}, product_role="desired product")
    outcome(r, "dr", "selectivity", [prod], [3], ">99:1", ">99:1", scope={"metric": "dr", "isomer_order": "major:minor"})
    outcome(r, "er", "selectivity", [prod], [3], ">99:1", ">99:1", scope={"metric": "er", "isomer_order": "major:minor"})
    REACTIONS.append(r)

# Scheme 3(b,c): scale-up and hydrolysis products.  Amino-acid 4a is not allene 4/4a.
for label, y in [("4a", 70), ("4b", 77), ("4c", 72), ("4d", 75), ("4e", 63)]:
    ensure(f"m_aminoacid_{label}", label, [3], name=f"α-alkenyl amino acid derivative {label}",
           chirality="partially_specified", geometry="specified", representation="not_encoded",
           double_bonds="E alkene drawn/reported", assignment="source drawing preserves the quaternary stereocenter and alkene geometry",
           note=("Amino-acid product 4a; distinct from 1-phenylallene 4/4a." if label == "4a" else "Hydrolysis product drawn in Scheme 3(c)."))
ensure("m_HCl", "HCl", [3], name="hydrochloric acid")
scale = reaction("scheme3b_scaleup_3a", [3], category="scale_up", name="large-scale alkyne hydroalkylation", transform="α-alkenylation by alkyne hydroalkylation", class_kind="hydroalkylation_alkyne")
ps = [participant("p_1a", "m_1a", "reactant", [3]), participant("p_2a", "m_2a", "reactant", [3]), participant("p_LiHMDS", "m_LiHMDS", "reagent", [3]), participant("p_THF", "m_THF", "solvent", [3]), participant("p_3a", "m_3a", "product", [3], product_role="desired product")]
charge(ps[0], "amount", "5.5 mmol", {"number": 5.5, "unit": "mmol"}, [3]); charge(ps[2], "mol_percent", "200 mol%", {"number": 200, "unit": "mol%"}, [3]); ps[2]["equivalence_reference"] = f("relative to 1a", "m_1a", [3])
scale["participants"] = ps; scale["reaction_scale"] = f("1a, 5.5 mmol", {"reference_molecule_id": "m_1a", "amount": {"number": 5.5, "unit": "mmol"}}, [3]); set_solvent(scale, [3], "THF")
step(scale, "reaction", [3], [p["participant_id"] for p in ps[:-1]], operation="large-scale hydroalkylation", temperature=("80 °C", {"number": 80, "unit": "°C"}), time=("20 h", {"number": 20, "unit": "h"}))
outcome(scale, "yield", "yield", ["m_3a"], [3], "74% isolated yield", {"number": 74, "unit": "%"}, product_role="desired product", yield_type="isolated yield", basis="isolated product")
outcome(scale, "amount", "product_amount", ["m_3a"], [3], "1.22 g", {"number": 1.22, "unit": "g"}, product_role="desired product")
outcome(scale, "dr", "selectivity", ["m_3a"], [3], ">99:1", ">99:1", scope={"metric": "dr", "isomer_order": "major:minor"})
outcome(scale, "er", "selectivity", ["m_3a"], [3], ">99:1", ">99:1", scope={"metric": "er", "isomer_order": "major:minor"})
REACTIONS.append(scale)

hyd = reaction("scheme3b_hydrolysis_4a", [3], category="derivatization", name="acidic hydrolysis of imidazolidinone 3a", transform="release of α-alkenyl amino acid derivative", class_kind="unresolved")
ps = [participant("p_3a", "m_3a", "reactant", [3]), participant("p_HCl", "m_HCl", "reagent", [3]), participant("p_EtOH", "m_EtOH", "solvent", [3]), participant("p_4a", "m_aminoacid_4a", "product", [3], product_role="desired product")]
hyd["participants"] = ps; set_solvent(hyd, [3], "6 M aqueous HCl/EtOH (10:1)")
step(hyd, "reaction", [3], [p["participant_id"] for p in ps[:-1]], operation="microwave acidic hydrolysis", temperature=("160 °C", {"number": 160, "unit": "°C"}), time=("4 h", {"number": 4, "unit": "h"}))
hyd["special_conditions"]["reactor"] = f("MW", "microwave reactor", [3])
outcome(hyd, "yield", "yield", ["m_aminoacid_4a"], [3], "70%", {"number": 70, "unit": "%"}, product_role="desired product")
outcome(hyd, "amount", "product_amount", ["m_aminoacid_4a"], [3], "0.63 g", {"number": 0.63, "unit": "g"}, product_role="desired product")
REACTIONS.append(hyd)

for label, precursor, y in [("4b", "m_3ac", 77), ("4c", "m_3s", 72), ("4d", "m_3q", 75), ("4e", "m_3t", 63)]:
    r = reaction(f"scheme3c_hydrolysis_{label}", [3], category="derivatization", name="hydrolysis of an imidazolidinone product", transform="release of α-alkenyl amino acid derivative", class_kind="unresolved")
    r["participants"] = [participant("p_precursor", precursor, "reactant", [3]), participant("p_product", f"m_aminoacid_{label}", "product", [3], product_role="desired product")]
    step(r, "reaction", [3], ["p_precursor"], operation="hydrolysis; detailed conditions not reported on supplied main-text page")
    outcome(r, "yield", "yield", [f"m_aminoacid_{label}"], [3], f"{y}%", {"number": y, "unit": "%"}, product_role="desired product")
    r["review"]["unresolved"].append("Scheme 3(c) gives product/yield but not the full hydrolysis conditions or yield measurement basis.")
    REACTIONS.append(r)

# Scheme 4 allene scope.
ensure("m_allene_internal_me", "internal allene 4, R2 = CH3", [3], name="internal phenylallene with R2 = CH3")
ensure("m_allene_internal_benzyl", "internal allene 4, R2 = CH2Ph", [3], name="internal phenylallene with R2 = CH2Ph")
for label in ["3aj", "3ak"]:
    ensure(f"m_{label}", label, [3], name=f"allene hydroalkylation product {label}", chirality="partially_specified", geometry="specified", representation="not_encoded", double_bonds="E alkene drawn/reported")
allene_scope = [("3a", "m_1a", "m_allene_4", 85, ">99:1", ">99:1"), ("3q", "m_1_for_3q", "m_allene_4", 88, ">99:1", ">99:1"),
                ("3v", "m_1_for_3v", "m_allene_4", 78, ">99:1", ">99:1"), ("3y", "m_1_for_3y", "m_allene_4", 81, ">99:1", ">99:1"),
                ("3ac", "m_1_for_3ac", "m_allene_4", 86, ">99:1", None), ("3aj", "m_1a", "m_allene_internal_me", 74, ">99:1", None),
                ("3ak", "m_1a", "m_allene_internal_benzyl", 71, ">99:1", None)]
for label, imid, allene, y, dr, er in allene_scope:
    r = reaction(f"scheme4_{label}", [3], category="substrate_scope", name="hydroalkylation of an allene with an imidazolidinone derivative", transform="α-alkenylation by allene hydroalkylation", class_kind="hydroalkylation_allene", general=(f"{DOC}:scheme1b_generic",))
    ps = [participant("p_imid", imid, "reactant", [3]), participant("p_allene", allene, "reactant", [3]), participant("p_LiHMDS", "m_LiHMDS", "reagent", [3]), participant("p_THF", "m_THF", "solvent", [3]), participant("p_product", f"m_{label}", "product", [3], product_role="desired product")]
    charge(ps[0], "amount", "0.2 mmol", {"number": 0.2, "unit": "mmol"}, [3]); charge(ps[1], "amount", "0.2 mmol", {"number": 0.2, "unit": "mmol"}, [3]); charge(ps[2], "volume", "0.2 mL, 1.0 M in THF", {"number": 0.2, "unit": "mL"}, [3]); charge(ps[2], "mol_percent", "100 mol%", {"number": 100, "unit": "mol%"}, [3]); ps[2]["equivalence_reference"] = f("relative to imidazolidinone", imid, [3])
    r["participants"] = ps; r["reaction_scale"] = f("imidazolidinone (0.2 mmol)", {"reference_molecule_id": imid, "amount": {"number": 0.2, "unit": "mmol"}}, [3])
    set_solvent(r, [3], "THF", volume=("2.0 mL", {"number": 2.0, "unit": "mL"}), concentration=("0.09 M", {"number": 0.09, "unit": "M"}), basis="reported reaction concentration")
    step(r, "reaction", [3], [p["participant_id"] for p in ps[:-1]], operation="Scheme 4 allene hydroalkylation", temperature=("80 °C", {"number": 80, "unit": "°C"}), time=("20 h", {"number": 20, "unit": "h"}))
    outcome(r, "yield", "yield", [f"m_{label}"], [3], f"{y}%", {"number": y, "unit": "%"}, product_role="desired product", yield_type="isolated yield", basis="isolated product")
    outcome(r, "dr", "selectivity", [f"m_{label}"], [3], dr, dr, scope={"metric": "dr", "isomer_order": "major:minor"})
    if er: outcome(r, "er", "selectivity", [f"m_{label}"], [3], er, er, scope={"metric": "er", "isomer_order": "major:minor"})
    else: r["review"]["unresolved"].append("Scheme 4 footnote a: er not determined.")
    REACTIONS.append(r)

# Page 4 mechanistic studies and proposed mechanism.
for mid, label, name in [
    ("m_allenyl_II", "II", "dilithiated allenyl species II"), ("m_enolate_III", "III", "lithium enolate III"),
    ("m_2a_D65", "2a, 65% D", "deuterium-labelled 2a mixture"), ("m_nBuLi", "n-BuLi", "n-butyllithium"),
    ("m_CD3OD", "CD3OD", "methanol-d4"), ("m_1d", "1d", "imidazolidinone 1d"),
    ("m_1d_D", "1d-D, >99% D", "deuterated imidazolidinone 1d-D"), ("m_2a_D3", "2a-D3, >99% D", "methyl-d3 alkyne 2a-D3"),
    ("m_D2O", "D2O", "deuterium oxide"),
]: ensure(mid, label, [4], name=name, mixture=("isotopomer mixture" if mid == "m_2a_D65" else None))

mech = reaction("figure1a_proposed_mechanism", [4], category="mechanistic_study", name="proposed hydroalkylation reaction mechanism", transform="author-proposed mechanism", class_kind="not_applicable", record_type="proposed_mechanism")
mech["participants"] = [participant("p_1a", "m_1a", "reactant", [4]), participant("p_2a", "m_2a", "reactant", [4]), participant("p_LiHMDS", "m_LiHMDS", "reagent", [4]), participant("p_II", "m_allenyl_II", "product", [4]), participant("p_III", "m_enolate_III", "product", [4]), participant("p_3a", "m_3a", "product", [4])]
note(mech, [4], "Species I–IV and arrows are proposed mechanistic intermediates/pathways, not independently isolated syntheses.")
REACTIONS.append(mech)

ir = reaction("figure1b_in_situ_ir", [4], category="mechanistic_study", name="in situ IR experiment on the hydroalkylation reaction", transform="in situ monitoring of hydroalkylation", class_kind="unresolved")
ir["participants"] = [participant("p_1a", "m_1a", "reactant", [4]), participant("p_2a", "m_2a", "reactant", [4]), participant("p_LiHMDS", "m_LiHMDS", "reagent", [4]), participant("p_THF", "m_THF", "solvent", [4]), participant("p_II", "m_allenyl_II", "product", [4]), participant("p_III", "m_enolate_III", "product", [4]), participant("p_3a", "m_3a", "product", [4])]
charge(ir["participants"][2], "equivalents", "2 eq.", {"number": 2, "unit": "equiv"}, [4]); ir["participants"][2]["equivalence_reference"] = f("relative to substrates", "m_1a", [4]); set_solvent(ir, [4], "THF")
step(ir, "reaction_monitoring", [4], ["p_1a", "p_2a", "p_LiHMDS", "p_THF"], operation="monitor reaction by in situ IR", temperature=("80 °C", {"number": 80, "unit": "°C"}))
for oid, target, peak, assignment in [("peak_1926", "m_allenyl_II", 1926, "dilithiated allenyl species"), ("peak_1666", "m_enolate_III", 1666, "enolate of 1a"), ("peak_1712", "m_3a", 1712, "product carbonyl")]:
    outcome(ir, oid, "other", [target], [4], f"~{peak} cm-1", f"~{peak} cm-1", basis="in situ IR", scope={"assignment": assignment})
REACTIONS.append(ir)

lih_quench = reaction("figure1a_alkyne_lihmds_cd3od", [4], category=["control", "mechanistic_study"], name="LiHMDS treatment and CD3OD quench of alkyne 2a", transform="reversible alkyne–allenyl lithiation followed by deuteration", class_kind="deuterium_exchange")
lih_quench["participants"] = [participant("p_2a", "m_2a", "reactant", [4]), participant("p_LiHMDS", "m_LiHMDS", "reagent", [4]), participant("p_CD3OD", "m_CD3OD", "reagent", [4]), participant("p_product", "m_2a_D65", "product", [4], product_role="regenerated deuterium-labelled alkyne")]
step(lih_quench, "base_treatment", [4], ["p_2a", "p_LiHMDS"], operation="treat 2a with LiHMDS to form equilibrating allenyl species", temperature=("80 °C", {"number": 80, "unit": "°C"}))
step(lih_quench, "quench", [4], ["p_CD3OD"], operation="quench with CD3OD")
outcome(lih_quench, "deuterium_incorporation", "other", ["m_2a_D65"], [4], "65% D in the methyl group", {"number": 65, "unit": "%"}, basis="NMR and mass analysis", scope={"metric": "D incorporation", "site": "methyl group"})
REACTIONS.append(lih_quench)

nbu = reaction("figure1b_nBuLi_ir_control", [4], category=["control", "mechanistic_study"], name="n-BuLi IR control with alkyne 2a", transform="formation of an allenyl lithium species", class_kind="unresolved")
nbu["participants"] = [participant("p_2a", "m_2a", "reactant", [4]), participant("p_nBuLi", "m_nBuLi", "reagent", [4]), participant("p_II", "m_allenyl_II", "product", [4])]
charge(nbu["participants"][1], "equivalents", "2 eq.", {"number": 2, "unit": "equiv"}, [4]); nbu["participants"][1]["equivalence_reference"] = f("relative to 2a", "m_2a", [4])
step(nbu, "ir_control", [4], ["p_2a", "p_nBuLi"], operation="add n-BuLi to 2a and observe IR band")
outcome(nbu, "peak_1926", "other", ["m_allenyl_II"], [4], "~1926 cm-1", "~1926 cm-1", basis="IR", scope={"assignment": "allenyl lithium species"})
nbu["review"]["unresolved"].append("The prose sequence does not unambiguously assign the subsequent 65% D quench result to the n-BuLi or LiHMDS IR mixture; it is not duplicated here.")
REACTIONS.append(nbu)

# Figure 1(c), four table rows.
d_rows = [(1, "m_1d", "m_2a_D3", 68, 82, False), (2, "m_1d", "m_2a", 14, 0, False),
          (3, "m_1d_D", "m_2a_D3", 82, 69, False), (4, "m_1d_D", "m_2a_D3", 82, 68, True)]
for entry, imid, alkyne, x, y, d2o in d_rows:
    prod = ensure(f"m_D_product_entry{entry}", f"deuterium-labelled product, Fig. 1c entry {entry}", [4],
                  name=f"isotopomer mixture from Fig. 1(c), entry {entry}", chirality="partially_specified", geometry="specified", representation="not_encoded", double_bonds="E alkene drawn", mixture="site-specific deuterium incorporation percentages x and y")
    r = reaction(f"figure1c_entry{entry}", [4], category="mechanistic_study", name="deuterium-labeling hydroalkylation experiment", transform="isotope-labelled alkyne hydroalkylation", class_kind="hydroalkylation_alkyne")
    ps = [participant("p_imid", imid, "reactant", [4]), participant("p_alkyne", alkyne, "reactant", [4]), participant("p_LiHMDS", "m_LiHMDS", "reagent", [4]), participant("p_THF", "m_THF", "solvent", [4])]
    if d2o: ps.append(participant("p_D2O", "m_D2O", "reagent", [4]))
    ps.append(participant("p_product", prod, "product", [4], product_role="deuterium-labelled product")); r["participants"] = ps
    charge(ps[2], "equivalents", "2 eq.", {"number": 2, "unit": "equiv"}, [4]); ps[2]["equivalence_reference"] = f("relative to imidazolidinone", imid, [4]); set_solvent(r, [4], "THF")
    step(r, "reaction", [4], [p["participant_id"] for p in ps if p["role"] != "product" and p["molecule_id"] != "m_D2O"], operation="deuterium-labeling hydroalkylation", temperature=("80 °C", {"number": 80, "unit": "°C"}), time=("20 h", {"number": 20, "unit": "h"}))
    if d2o: step(r, "quench", [4], ["p_D2O"], operation="quench with D2O")
    outcome(r, "D_site_x", "other", [prod], [4], f"x = {x}% D", {"number": x, "unit": "%"}, scope={"metric": "D incorporation", "site": "source-marked vinylic position x"})
    outcome(r, "D_site_y", "other", [prod], [4], f"y = {y}% D", {"number": y, "unit": "%"}, scope={"metric": "D incorporation", "site": "source-marked methyl position y"})
    REACTIONS.append(r)

mei = reaction("figure1c_methyl_iodide_quench", [4], category=["control", "mechanistic_study"], name="methyl iodide quench control", transform="quench test for organolithium intermediate", class_kind="not_applicable")
mei["participants"] = [participant("p_1d", "m_1d", "reactant", [4]), participant("p_2a", "m_2a", "reactant", [4]), participant("p_LiHMDS", "m_LiHMDS", "reagent", [4]), participant("p_MeI", "m_MeI", "reagent", [4])]
step(mei, "reaction", [4], ["p_1d", "p_2a", "p_LiHMDS"], operation="run hydroalkylation mixture")
step(mei, "quench", [4], ["p_MeI"], operation="quench with methyl iodide")
outcome(mei, "no_alkylated_product", "no_reaction", [], [4], "did not result in the alkylated product", "no alkylated product observed", scope={"target": "methyl-iodide-trapped product"})
REACTIONS.append(mei)

# Page 5 footnote 13: two additive controls, no loading/yield basis supplied.
for mid, label, name, ee in [("m_HMPA", "HMPA", "hexamethylphosphoramide", 73), ("m_12crown4", "[12]crown-4", "12-crown-4", 51)]:
    ensure(mid, label, [5], name=name)
    r = reaction(f"footnote13_{label.replace('[','').replace(']','').replace('-','_')}", [5], category=["control", "optimization"], name="lithium-chelation additive control", transform="additive-modified alkyne hydroalkylation", class_kind="hydroalkylation_alkyne")
    r["participants"] = [participant("p_1a", "m_1a", "reactant", [5]), participant("p_2a", "m_2a", "reactant", [5]), participant("p_LiHMDS", "m_LiHMDS", "reagent", [5]), participant("p_additive", mid, "additive", [5]), participant("p_3a", "m_3a", "product", [5], product_role="desired product")]
    step(r, "reaction", [5], ["p_1a", "p_2a", "p_LiHMDS", "p_additive"], operation=f"use {label} as an additive to capture lithium")
    outcome(r, "ee", "selectivity", ["m_3a"], [5], f"ee decreased to {ee}%", {"number": ee, "unit": "%"}, scope={"metric": "ee", "major_enantiomer": "not assigned in footnote"})
    r["review"]["unresolved"].append("Additive loading, reaction yield, temperature/time, and ee measurement method are not reported in footnote 13.")
    REACTIONS.append(r)

# Non-reaction observations and cited context.
SOURCE_CONTEXT = [
    {"context_id": f"{DOC}:scheme1a_prior_work", "document_id": DOC, "type": "cited_prior_work",
     "location": {"page": 1, "region": "Scheme 1(a)"}, "evidence_ids": [EVI[1], EVT[1]],
     "note": "Clayden et al. intramolecular N→C migration is cited prior work, not a this-paper experiment."},
    {"context_id": f"{DOC}:figure1d_dft", "document_id": DOC, "type": "computational_mechanism",
     "location": {"page": 4, "region": "Figure 1(d)"}, "evidence_ids": [EVI[4], EVT[4]],
     "note": "DFT transition-state structures and relative free energies; not synthetic experiments."},
    {"context_id": f"{DOC}:xray_3ac", "document_id": DOC, "type": "characterization",
     "location": {"page": 3, "region": "Scheme 4 X-ray inset"}, "evidence_ids": [EVI[3], EVT[3]],
     "note": "Single-crystal X-ray analysis of 3ac confirmed relative and absolute configuration; CCDC 1950247 is stated on page 1."},
    {"context_id": f"{DOC}:page5_bibliography", "document_id": DOC, "type": "bibliography_and_conclusion",
     "location": {"page": 5}, "evidence_ids": [EVI[5], EVT[5]],
     "note": "Conclusion, authorship, acknowledgements, and references; no additional synthetic arrow beyond footnote 13 controls."},
]

# Preserve raw visual transcription and deterministic stereo diagnostics separately.
raw = {
    "document_id": DOC,
    "pages_inspected": [1, 2, 3, 4, 5],
    "machine_generated": True,
    "semantic_accuracy_verified": False,
    "inventory": {
        "page1": ["Scheme 1(a) cited prior work", "Scheme 1(b) this-work generic transformation", "article metadata and prose"],
        "page2": ["Table 1 entries 1-8", "Scheme 2 products 3a-ai (35 scope arrows under one general arrow)", "scope discussion"],
        "page3": ["Scheme 3(a) seven arrows", "Scheme 3(b) scale-up plus hydrolysis arrow", "Scheme 3(c) products 4a-e", "Scheme 4 products 3a,3q,3v,3y,3ac,3aj,3ak", "X-ray inset for 3ac"],
        "page4": ["Figure 1(a) proposed mechanism and LiHMDS/CD3OD control", "Figure 1(b) in situ IR", "Figure 1(c) four deuterium-labeling rows", "methyl iodide quench prose", "Figure 1(d) DFT"],
        "page5": ["conclusion and references", "footnote 13 HMPA and 12-crown-4 controls"],
    },
    "reaction_ids": [r["reaction_id"] for r in REACTIONS],
    "limitations": ["No SI supplied or associated.", "Compact structures remain unresolved where a stereochemically safe SMILES could not be established."],
}
(RUN / "raw_model/page1_5_visual_transcription.json").write_text(json.dumps(raw, ensure_ascii=False, indent=2) + "\n")
diagnostics = {mid: st.stereo_inspect(st.value(m["smiles"])) for mid, m in MOLECULES.items() if st.value(m["smiles"])}
(RUN / "raw_model/standard_stereo_diagnostics.json").write_text(json.dumps(diagnostics, ensure_ascii=False, indent=2) + "\n")

bundle = st.empty(DATASET["documents"])
bundle["molecules"] = list(MOLECULES.values())
bundle["reactions"] = REACTIONS
bundle["source_context"] = SOURCE_CONTEXT
(RUN / "submitted.json").write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")
print(json.dumps({"molecules": len(MOLECULES), "reactions": len(REACTIONS), "source_context": len(SOURCE_CONTEXT)}, ensure_ascii=False))
