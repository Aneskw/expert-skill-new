#!/usr/bin/env python3
"""Build the source-bounded v3.1 submission for doc_27ee161674b4e4ee.

This is a source-specific transcription aid.  It deliberately leaves compact
structure drawings, SI-only procedures, and unreported measurement methods
unresolved.  It does not alter schemas or extraction tools.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

RUN = ROOT / "runs/current"
RAW = RUN / "raw_model"
DOC = "doc_27ee161674b4e4ee"
P = {n: {"t": f"{DOC}:page{n}:text", "i": f"{DOC}:page{n}:image"} for n in range(1, 6)}
MISSING = object()


def f(raw=None, normalized=MISSING, ev=(), method="original_extraction", status=None,
      derivation=None, alternatives=()):
    if normalized is MISSING:
        normalized = raw
    return st.fact(raw, normalized, ev, method=method, status=status,
                   derivation=derivation, alternatives=alternatives)


bundle = json.loads((RAW / "standard_init.json").read_text())
molecules = {}
stereo_diagnostics = {}


def ensure_molecule(mid, label, ev, *, smiles=None, original_name=None,
                    chirality="unresolved", geometry="unresolved",
                    representation=None, axial=None, atom_centers=None,
                    double_bonds=None, mixture=None, assignment=None,
                    identity_note=None, template=None):
    ev = list(dict.fromkeys(ev))
    if mid in molecules:
        occ = molecules[mid]["occurrences"][0]
        occ["evidence_ids"] = list(dict.fromkeys(occ["evidence_ids"] + ev))
        return molecules[mid]
    m = st.new_molecule(mid, DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"] = [{"document_id": DOC,
                          "compound_label": f(label, label, ev),
                          "evidence_ids": ev}]
    name = original_name or label
    m["names"]["original"] = f(name, name, ev)
    m["names"]["english"] = f()
    m["names"]["abbreviations"] = f()
    m["names"]["aliases"] = f()
    m["cas_numbers"] = f()
    if template is not None:
        m["representation_type"] = "generic_template"
        m["smiles"] = f()
        m["generic_template"] = template
        representation = representation or "not_encoded"
    elif smiles:
        diag = st.stereo_inspect(smiles)
        stereo_diagnostics[mid] = diag
        m["representation_type"] = "smiles"
        m["smiles"] = f({"source": "structure drawing", "candidate_smiles": smiles}, smiles, ev,
                         method="structure_derivation",
                         derivation="Candidate SMILES transcribed from the source drawing; parsing does not verify identity.")
        m["generic_template"] = None
        representation = representation or ("partial" if diag["potential_stereo"] else "encoded")
        m["stereochemistry"]["isomeric_smiles"] = f(
            {"candidate_smiles": diag["canonical_isomeric_smiles"]}, diag["canonical_isomeric_smiles"], ev,
            method="structure_derivation", derivation="RDKit canonical isomeric SMILES for the candidate graph; source configuration is not inferred.")
        m["stereochemistry"]["identity_key"] = f(
            None, diag["stereo_graph_key"], ev, method="structure_derivation",
            derivation="Stereo graph key calculated from the candidate SMILES; not a sample identity.")
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = f()
        m["generic_template"] = None
        representation = representation or "not_encoded"
    m["stereochemistry"]["atom_centers"] = f(atom_centers, atom_centers, ev) if atom_centers is not None else f()
    m["stereochemistry"]["double_bonds"] = f(double_bonds, double_bonds, ev) if double_bonds is not None else f()
    m["stereochemistry"]["axial"] = f(axial, axial, ev) if axial is not None else f()
    m["stereochemistry"]["mixture"] = f(mixture, mixture, ev) if mixture is not None else f()
    m["stereochemistry"]["assignment"] = f(assignment, assignment, ev) if assignment is not None else f()
    for key, value in (("chirality_status", chirality), ("geometry_status", geometry),
                       ("representation_status", representation)):
        m["stereochemistry"][key] = f(value, value, ev, method="inference",
                                        derivation="Conservative status from the inspected drawing and explicit source descriptors; absence of a descriptor is not achirality or racemic composition.")
    if not smiles:
        m["stereochemistry"]["isomeric_smiles"] = f()
        m["stereochemistry"]["identity_key"] = f()
    m["identity_note"] = f(identity_note, identity_note, ev) if identity_note else f()
    m["external_properties"] = []
    molecules[mid] = m
    return m


def generic_template(description, variables, ev):
    return {
        "fixed_structure": f(description, {"source_description": description}, ev),
        "variables": [{
            "label": v,
            "attachment_points": f(f"Source drawing attachment for {v}", [v], ev),
            "definition": f(f"Variable {v} as drawn; individual allowed values are not exhaustively defined here.",
                            f"Variable {v} as drawn", ev),
            "allowed_values": f(),
            "constraints": f(),
        } for v in variables],
        "instance_mapping": [],
        "reason": f("R-group/general drawing", "Source is a general family, not one concrete molecule.", ev),
    }


def participant(pid, mid, role, ev, *, stage="reaction", product_role=None,
                mass=None, amount=None, volume=None, equivalents=None,
                mol_percent=None, reference=None):
    p = st.new_participant(pid, mid, role, ev, stage=stage)
    if product_role:
        p["product_role"] = f(product_role, product_role, ev)
    for key, item in (("mass", mass), ("amount", amount), ("volume", volume),
                      ("equivalents", equivalents), ("mol_percent", mol_percent)):
        if item is not None:
            raw, norm = item
            p["charges"][key] = f(raw, norm, ev)
    if reference is not None:
        p["equivalence_reference"] = f(reference, reference, ev)
    return p


def reaction(rid, ev, record_type="specific_experiment", categories=()):
    r = st.new_reaction(f"{DOC}:{rid}", DOC, ev)
    r["record_type"] = record_type
    r["experiment_categories"] = f(list(categories), list(categories), ev)
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                   "checked_fields": [], "unresolved": []}
    r["completeness"] = "partial"
    return r


def set_location(r, raw, page, figure):
    r["location"] = f(raw, {"page": page, "figure": figure}, r["evidence_ids"])


def set_fact_name(r, name, transform=None, inferred=False):
    r["reaction_name"] = f(name, name, r["evidence_ids"])
    if transform is None:
        transform = name
    if inferred:
        r["transformation_type"] = f(None, transform, r["evidence_ids"], method="inference",
                                      derivation="Class-level transformation inferred from the source arrow and prose.")
    else:
        r["transformation_type"] = f(transform, transform, r["evidence_ids"])


def add_step(r, sid, ev, participants=(), *, operation=None, temperature=None,
             time=None, atmosphere=None, pressure=None, stage="reaction"):
    s = st.new_step(sid, ev, stage=stage)
    s["participant_ids"] = list(participants)
    if operation is not None:
        s["operation"] = f(operation, operation, ev)
    if temperature is not None:
        s["temperature"] = f(temperature[0], temperature[1], ev)
    if time is not None:
        s["time"] = f(time[0], time[1], ev)
    if atmosphere is not None:
        s["atmosphere"] = f(atmosphere, atmosphere, ev)
    if pressure is not None:
        s["pressure"] = f(pressure[0], pressure[1], ev)
    r["steps"].append(s)


def add_outcome(r, oid, kind, targets, ev, raw, norm=MISSING, *, product_role=None,
                yield_type=None, basis=None, scope=None):
    o = st.new_outcome(oid, kind, targets, ev)
    if norm is MISSING:
        norm = raw
    o["value"] = f(raw, norm, ev)
    if product_role:
        o["product_role"] = f(product_role, product_role, ev)
    if yield_type:
        o["yield_type"] = f(yield_type, yield_type, ev)
    if basis:
        o["measurement_basis"] = f(basis, basis, ev)
    if scope:
        o["scope"] = f(scope, scope, ev)
    r["outcomes"].append(o)
    return o


def classification(r, kind):
    schemes = {
        "alkynylalkenylation": [
            ("b1", "C", "C", 0, 1, "formed"),
            ("b2", "C", "C", 0, 1, "formed"),
            ("b3", "C", "O", 0, 1, "formed"),
            ("b4", "C", "C", 2, 3, "order_increased"),
            ("b5", "C", "Br", 1, 0, "broken"),
            ("b6", "C", "S", 1, 0, "broken"),
        ],
        "enynylalkenylation": [
            ("b1", "C", "C", 0, 1, "formed"),
            ("b2", "C", "C", 0, 1, "formed"),
            ("b3", "C", "C", 2, 3, "order_increased"),
            ("b4", "C", "Br", 1, 0, "broken"),
            ("b5", "C", "S", 1, 0, "broken"),
        ],
        "cc_only": [("b1", "C", "C", 0, 1, "formed")],
        "epoxidation": [
            ("b1", "C", "O", 0, 1, "formed"),
            ("b2", "C", "O", 0, 1, "formed"),
            ("b3", "C", "C", 2, 1, "order_decreased"),
            ("b4", "O", "O", 1, 0, "broken"),
        ],
        "hydrogenation": [
            ("b1", "C", "C", 3, 2, "order_decreased"),
            ("b2", "C", "H", 0, 1, "hydrogen_gain"),
            ("b3", "C", "H", 0, 1, "hydrogen_gain"),
        ],
        "deoxygenation": [
            ("b1", "C", "O", 1, 0, "broken"),
            ("b2", "C", "H", 0, 1, "hydrogen_gain"),
        ],
        "allene_formation": [
            ("b1", "C", "O", 1, 0, "broken"),
            ("b2", "C", "C", 3, 2, "order_decreased"),
            ("b3", "C", "C", 1, 2, "order_increased"),
            ("b4", "C", "H", 0, 1, "hydrogen_gain"),
        ],
        "cyanoalkylation": [
            ("b1", "C", "C", 0, 1, "formed"),
            ("b2", "C", "C", 1, 2, "order_increased"),
        ],
    }
    if kind in ("not_applicable", "unresolved"):
        r["classification"]["bond_changes"] = []
        r["classification"]["labels"] = f(None, [], r["evidence_ids"], method="inference",
                                              derivation="No observed synthetic bond formation is asserted for this record.")
        r["classification"]["coverage"] = f(None, kind, r["evidence_ids"], method="inference",
                                                derivation="Mechanistic measurement/proposal or insufficient source detail.")
        r["classification"]["limitations"] = ["No expected-product bond labels were invented for a non-transformative or unresolved record."]
        r["classification"]["model_votes"] = []
        return
    changes = []
    labels = {}
    for cid, a, b, before, after, ctype in schemes[kind]:
        change = {"change_id": cid, "atoms": [{"atom_map": None, "element": a},
                                                {"atom_map": None, "element": b}],
                  "before_order": before, "after_order": after, "change_type": ctype}
        changes.append(f(None, change, r["evidence_ids"], method="inference",
                         derivation="Class-level bond change inferred from the inspected source drawings; exact atom mapping is unresolved."))
        if ctype in ("formed", "order_increased", "hydrogen_gain"):
            lab = st.bond_label(a, b, after)
            labels.setdefault(lab["code"], {**lab, "supporting_change_ids": []})["supporting_change_ids"].append(cid)
    r["classification"]["bond_changes"] = changes
    r["classification"]["labels"] = f(None, list(labels.values()), r["evidence_ids"], method="inference",
                                          derivation="All class-level formed/promoted bonds listed above determine the taxonomy labels.")
    r["classification"]["coverage"] = f(None, "partial", r["evidence_ids"], method="inference",
                                            derivation="Exact atom identities, complete hydrogen changes, coproducts, and atom mapping are unresolved.")
    r["classification"]["limitations"] = [
        "Class-level inference from source drawings; no atom-mapped reaction is asserted.",
        "Hydrogen source, coproducts, and all elementary bond changes are not established on the supplied pages.",
    ]
    r["classification"]["model_votes"] = []


def set_common_photo(r, ev, base_mid, substrate_mid, product_mid, *, base_amount="0.2 mmol",
                     substrate_eq=1.5, catalyst_pct=1, kh_eq=2, solvents=("MeCN", "MeOH"),
                     volumes=(0.5, 1.5), power=5, time_h=None, catalyst_mid="fac_irppy3"):
    solvent_ids = {"MeCN": "m_me_c_n", "MeOH": "m_meoh", "DMF": "m_dmf",
                   "cyclohexanol": "m_cyclohexanol", "5-hexyn-1-ol": "m_5_hexyn_1_ol",
                   "water": "m_water"}
    parts = [
        participant("p_base", base_mid, "reactant", ev, amount=(base_amount, {"number": float(base_amount.split()[0]), "unit": "mmol"})),
        participant("p_substrate", substrate_mid, "reactant", ev,
                    equivalents=(f"{substrate_eq} equiv", {"number": substrate_eq, "unit": "equiv"}), reference=base_mid),
        participant("p_photocatalyst", catalyst_mid, "catalyst", ev,
                    mol_percent=(f"{catalyst_pct} mol%", {"number": catalyst_pct, "unit": "mol%"}), reference=base_mid),
        participant("p_base_additive", "m_kh2po4", "additive", ev,
                    equivalents=(f"{kh_eq} equiv", {"number": kh_eq, "unit": "equiv"}), reference=base_mid),
    ]
    for i, (name, vol) in enumerate(zip(solvents, volumes), 1):
        parts.append(participant(f"p_solvent_{i}", solvent_ids[name], "solvent", ev,
                                 volume=(f"{vol} mL", {"number": vol, "unit": "mL"})))
    parts.append(participant("p_product", product_mid, "product", ev, product_role="main product"))
    r["participants"] = parts
    r["reaction_scale"] = f(base_amount, {"reference_molecule_id": base_mid,
                                           "amount": {"number": float(base_amount.split()[0]), "unit": "mmol"}}, ev)
    r["solvent_conditions"]["composition"] = f("/".join(solvents), list(solvents), ev)
    r["solvent_conditions"]["mixture_ratio"] = f(
        "/".join(str(x) for x in volumes) + " mL", {"component_volumes_mL": list(volumes)}, ev)
    r["solvent_conditions"]["total_volume"] = f(
        f"{sum(volumes)} mL (calculated from listed component volumes)",
        {"number": sum(volumes), "unit": "mL"}, ev, method="inference",
        derivation="Arithmetic sum of individually reported solvent volumes; not an independently reported total volume.")
    pids = [p["participant_id"] for p in parts if p["role"] != "product"]
    add_step(r, "reaction", ev, pids, operation="irradiated with blue LEDs",
             temperature=("rt", "rt"),
             time=(f"{time_h} h", {"number": time_h, "unit": "h"}) if time_h is not None else None,
             atmosphere="N2")
    r["special_conditions"]["light_source"] = f("blue LEDs", "blue LEDs", ev)
    r["special_conditions"]["power"] = f(f"{power} W", {"number": power, "unit": "W"}, ev)
    r["procedure"]["original_text"] = f("Figure caption reaction conditions", "Figure caption reaction conditions", ev)


# Core molecules and condition agents.
p12 = [P[1]["i"], P[2]["i"]]
ensure_molecule("m_1a", "1a", p12, original_name="dual-function reagent 1a",
                smiles="CC(C)=C=C(Br)S(=O)(=O)C=Cc1ccccc1",
                chirality="achiral", geometry="unspecified", representation="partial",
                axial={"allene": "not stereogenic on candidate graph because one terminal carbon bears two methyl groups"},
                identity_note="Candidate graph is source-derived and unverified; vinyl geometry is not explicitly assigned in the supplied text.")
ensure_molecule("m_2a", "2a", [P[2]["i"], P[3]["i"]], original_name="styrene",
                smiles="C=Cc1ccccc1", chirality="achiral", geometry="not_applicable", representation="encoded")
ensure_molecule("m_3a", "3a", [P[2]["i"], P[4]["i"]], original_name="1,5-enyne product 3a",
                smiles="COC(C)(C)C#CC(c1ccccc1)CC=Cc1ccccc1",
                chirality="unspecified", geometry="unspecified", representation="partial",
                atom_centers=[{"candidate_smiles_atom_index": 7, "configuration": "unassigned"}],
                identity_note="The candidate graph has an unassigned tetrahedral center and unassigned alkene geometry; no racemic claim is made.")

for mid, label, ev in [
    ("fac_irppy3", "fac-Ir(ppy)3", [P[2]["i"], P[3]["i"]]),
    ("m_kh2po4", "KH2PO4", [P[2]["i"], P[3]["i"]]),
    ("m_me_c_n", "MeCN", [P[2]["i"], P[3]["i"]]),
    ("m_meoh", "MeOH", [P[2]["i"], P[3]["i"]]),
    ("m_dmf", "DMF", [P[3]["i"]]),
    ("m_water", "H2O", [P[2]["i"], P[4]["i"]]),
    ("m_cyclohexanol", "cyclohexanol", [P[2]["i"]]),
    ("m_5_hexyn_1_ol", "5-hexyn-1-ol", [P[2]["i"]]),
]:
    ensure_molecule(mid, label, ev, original_name=label)


# Scheme 1: two cited generic families and one present-work generic scheme.
lit_a_sub = "m_lit_a_substrate"
lit_a_prod = "m_lit_a_product"
lit_b_sub = "m_lit_b_substrate"
lit_b_prod = "m_lit_b_product"
generic_dfr = "m_generic_dfr"
generic_alkene = "m_generic_alkene"
generic_product = "m_generic_product"
for mid, label, desc, vars_ in [
    (lit_a_sub, "Scheme 1A(a) generic alkene", "Generic alkene R-substituted substrate", ["R"]),
    (lit_a_prod, "Scheme 1A(a) generic product", "Generic prior-work product bearing one unsaturated group and X", ["R", "R2", "X"]),
    (lit_b_sub, "Scheme 1A(b) generic tertiary alcohol", "Generic tertiary allylic/propargylic alcohol", ["R", "R2"]),
    (lit_b_prod, "Scheme 1A(b) generic migration product", "Generic prior-work migration product", ["R", "R2", "X"]),
    (generic_dfr, "Scheme 1B dual-function reagent", "Alkenyl and bromoallenyl groups tethered by sulfone", ["R1", "R2", "R3"]),
    (generic_alkene, "Scheme 1B alkene 2", "Generic alkene", ["R"]),
    (generic_product, "Scheme 1B products", "Generic 1,2-difunctionalization products", ["R", "R1", "R2", "R3", "R4"]),
]:
    ensure_molecule(mid, label, [P[1]["i"]], template=generic_template(desc, vars_, [P[1]["i"]]),
                    chirality="unresolved", geometry="unresolved", representation="not_encoded")

for rid, sub, prod, loc, name in [
    ("scheme1a_prior_cross_coupling", lit_a_sub, lit_a_prod, "Scheme 1A(a)", "Prior radical addition/cross-coupling or trapping family"),
    ("scheme1b_prior_migration", lit_b_sub, lit_b_prod, "Scheme 1A(b)", "Prior intramolecular alkenyl/alkynyl migration family"),
]:
    r = reaction(rid, [P[1]["i"], P[1]["t"]], record_type="literature_example")
    set_location(r, loc, 1, "Scheme 1A")
    set_fact_name(r, name, "cited prior-work generic transformation")
    r["participants"] = [participant("p_substrate", sub, "reactant", r["evidence_ids"]),
                         participant("p_product", prod, "product", r["evidence_ids"], product_role="illustrative generic product")]
    add_step(r, "conceptual_arrow", r["evidence_ids"], ["p_substrate"], operation="generic cited transformation; conditions not supplied")
    classification(r, "cc_only")
    r["review"]["unresolved"] = ["Generic cited literature family only; individual cited experiments and conditions are not represented on this page."]
    bundle["reactions"].append(r)

GENERIC_ID = f"{DOC}:scheme1_current_generic"
r = reaction("scheme1_current_generic", [P[1]["i"], P[1]["t"], P[2]["i"], P[3]["i"]], record_type="generic_scheme")
set_location(r, "Scheme 1B and repeated generic headers in Figures 2–4", 1, "Scheme 1B")
set_fact_name(r, "1,2-alkynylalkenylation / 1,2-enynylalkenylation", "radical docking-migration cascade")
r["participants"] = [participant("p_dfr", generic_dfr, "reactant", r["evidence_ids"]),
                     participant("p_alkene", generic_alkene, "reactant", r["evidence_ids"]),
                     participant("p_product", generic_product, "product", r["evidence_ids"], product_role="generic product family")]
add_step(r, "generic_arrow", r["evidence_ids"], ["p_dfr", "p_alkene"], operation="photocatalytic radical docking-migration cascade")
classification(r, "enynylalkenylation")
r["review"]["unresolved"] = ["R groups and product nucleophile/deprotonation branch remain generic; this is not a concrete experiment."]
bundle["reactions"].append(r)


# Figure 1/2 scope.  3a is the optimized example repeated in Figure 1 and Figure 2.
scope2 = [
    ("3a", 72, 47, None, None),
    ("3b", 65, 75, None, None), ("3c", 75, 36, None, None),
    ("3d", 71, 36, None, None), ("3e", 73, 36, None, None),
    ("3f", 71, 30, None, None), ("3g", 67, 38, None, None),
    ("3h", 72, 38, None, None), ("3i", 61, 36, None, None),
    ("3j", 66, 72, None, None), ("3k", 72, 44, None, None),
    ("3l", 70, 44, None, None), ("3m", 62, 75, None, None),
    ("3n", 70, 36, None, None), ("3o", 68, 72, None, None),
    ("3p", 52, 70, None, None), ("3q", 56, 84, None, None),
    ("3r", 69, 35, "cis", None), ("3s", 74, 35, "cis", None),
    ("3t", 43, 40, "dr >20:1", "a"), ("3u", 75, 40, None, None),
    ("3v", 58, 35, None, None), ("3w", 24, 36, None, None),
    ("3x", 36, 24, None, "a"), ("3y", 32, 36, "cis", "a"),
    ("3z", 47, 68, None, None), ("3aa", 65, 48, None, None),
    ("3ab", 50, 48, None, "a"), ("3ac", 50, 48, None, "a"),
    ("3ad", 43, 48, None, "a"), ("3ae", 49, 48, None, "a"),
    ("3af", 58, 24, "dr >20:1", "a"), ("3ag", 32, 48, "cis", "a"),
    ("3ah", 46, 36, "cis", "a"), ("3ai", 49, 62, "dr = 1:1", "a"),
    ("3aj", 55, 36, "dr = 1:1", "b"),
    ("4a", 41, 17, None, "c"), ("4b", 48, 26, None, "c"),
    ("4c", 60, 12, None, "ac"), ("4d", 46, 17, None, "ac"),
    ("4e", 46, 16, None, "ac"), ("4f", 35, 17, "dr = 1:1", "ac"),
]


def stereo_for_product(label, selectivity):
    if selectivity == "cis":
        return dict(chirality="partially_specified", geometry="not_applicable",
                    assignment={"relative_configuration": "cis", "absolute_configuration": "unassigned"})
    if selectivity and "dr" in selectivity:
        return dict(chirality="mixture", geometry="unspecified",
                    mixture={"type": "diastereomeric mixture", "reported_ratio": selectivity,
                             "major_configuration": "unassigned"})
    return dict(chirality="unspecified", geometry="unspecified")


for label, yld, hours, sel, foot in scope2:
    ev = [P[2]["i"], P[2]["t"]]
    if label in {"3ab", "3ac", "3ad", "3ae", "3af", "3ag", "3ah", "3ai", "3aj", "4a", "4b", "4c", "4d", "4e", "4f"}:
        ev.append(P[3]["t"])
    product_mid = f"m_product_{label}"
    if label == "3a":
        product_mid = "m_3a"
        substrate_mid = "m_2a"
    else:
        substrate_mid = f"m_alkene_for_{label}"
        ensure_molecule(substrate_mid, f"alkene substrate drawn for {label}", [P[2]["i"]],
                        original_name=f"unnumbered alkene substrate corresponding to {label}",
                        chirality="unresolved", geometry="unresolved", representation="not_encoded",
                        identity_note="Compact source drawing was not converted to a complete candidate SMILES.")
        kw = stereo_for_product(label, sel)
        ensure_molecule(product_mid, label, [P[2]["i"]], original_name=f"product {label}",
                        representation="not_encoded", identity_note="Structure is retained by source label and drawing; complete SMILES not asserted.", **kw)
    r = reaction(f"scope_{label}", ev, categories=("optimization", "substrate_scope") if label == "3a" else ("substrate_scope",))
    set_location(r, "Figure 1 optimized reaction and Figure 2 scope" if label == "3a" else f"Figure 2, product {label}", 2, "Figure 1/Figure 2" if label == "3a" else "Figure 2")
    set_fact_name(r, "1,2-alkynylalkenylation of alkenes")
    r["general_scheme_ids"] = [GENERIC_ID]
    substrate_eq = 3 if foot and "a" in foot else 2 if foot == "b" else 1.5
    catalyst_pct = 2 if foot and ("a" in foot or foot == "b") else 1
    if foot and "c" in foot:
        cosolvent = "cyclohexanol" if label == "4a" else "5-hexyn-1-ol" if label == "4b" else "water"
        solvents = ("MeCN", cosolvent)
        volumes = (1.5, 0.3)
    else:
        solvents, volumes = ("MeCN", "MeOH"), (0.5, 1.5)
    set_common_photo(r, ev, "m_1a", substrate_mid, product_mid, substrate_eq=substrate_eq,
                     catalyst_pct=catalyst_pct, kh_eq=2, solvents=solvents, volumes=volumes,
                     power=5, time_h=hours)
    add_outcome(r, "yield", "yield", [product_mid], ev, f"{yld}%", {"number": yld, "unit": "%"},
                product_role="main product")
    if sel:
        metric = "cis relative stereochemistry" if sel == "cis" else "diastereomeric ratio"
        add_outcome(r, "selectivity", "selectivity", [product_mid], ev, sel, sel,
                    scope={"metric": metric, "order": "unreported for dr" if "dr" in sel else "cis product", "major_configuration": "unassigned"})
    classification(r, "alkynylalkenylation")
    r["notes"] = [f("Yield measurement basis is not stated in Figure 1/2 or adjacent prose.",
                    "Yield measurement basis is not stated in Figure 1/2 or adjacent prose.", ev)]
    r["review"]["unresolved"] = [
        "Yield measurement basis is unreported on pages 1–5.",
        "Complete experimental procedure and characterization are in advertised SI, which was not supplied.",
        "Exact atom mapping and full product stereochemical assignment are unresolved.",
    ]
    bundle["reactions"].append(r)


# Additional optimization/scope experiments explicitly reported in page-2 prose
# but not assigned numbered structures in Figure 2.
for mid, label in [("m_na2co3", "Na2CO3"), ("m_k2co3", "K2CO3"),
                   ("m_et3n", "Et3N"), ("m_dipea", "DIPEA")]:
    ensure_molecule(mid, label, [P[2]["t"]], original_name=label)

page2_prose_specs = [
    ("optimization_insufficient_meoh", None, "Insufficient methanol resulted in a lower yield of 3a", "lower yield of 3a", "unresolved"),
    ("optimization_meoh_only", None, "Methanol as sole solvent entirely impeded the transformation", "transformation entirely impeded", "not_applicable"),
    ("optimization_na2co3", "m_na2co3", "Na2CO3 was not effective for the reaction", "not effective", "unresolved"),
    ("optimization_k2co3", "m_k2co3", "K2CO3 was not effective for the reaction", "not effective", "unresolved"),
    ("optimization_et3n", "m_et3n", "Et3N was not effective for the reaction", "not effective", "unresolved"),
    ("optimization_dipea", "m_dipea", "DIPEA was not effective for the reaction", "not effective", "unresolved"),
]
for rid, alternative_base, raw_obs, obs, ckind in page2_prose_specs:
    ev = [P[2]["t"], P[2]["i"]]
    r = reaction(rid, ev, categories=("optimization",))
    set_location(r, "Page 2 optimization discussion following Figure 1", 2, "prose after Figure 1")
    set_fact_name(r, "Photochemical optimization of 1a with styrene", "optimization experiment")
    parts = [participant("p_1a", "m_1a", "reactant", ev),
             participant("p_2a", "m_2a", "reactant", ev)]
    if alternative_base:
        parts.append(participant("p_alternative_base", alternative_base, "additive", ev))
    else:
        parts.append(participant("p_meoh", "m_meoh", "solvent", ev))
        if rid == "optimization_insufficient_meoh":
            parts.append(participant("p_mecn", "m_me_c_n", "solvent", ev))
    if rid == "optimization_insufficient_meoh":
        parts.append(participant("p_3a", "m_3a", "product", ev, product_role="desired product"))
    r["participants"] = parts
    add_step(r, "optimization_run", ev, [p["participant_id"] for p in parts if p["role"] != "product"],
             operation="optimization experiment; exact quantities, irradiation time, and solvent ratio are in unavailable SI")
    targets = ["m_3a"] if rid == "optimization_insufficient_meoh" else []
    add_outcome(r, "qualitative_result", "other" if ckind == "unresolved" else "no_reaction", targets, ev,
                raw_obs, obs, scope={"desired_transformation": "1a + styrene to 3a"})
    classification(r, ckind)
    r["notes"] = [f(raw_obs, obs, ev)]
    r["review"]["unresolved"] = [
        "The main text supplies only a qualitative optimization result; exact amounts, time, yield, and detection method are unavailable without SI.",
        "For 'not effective', no numerical yield or strict no-reaction claim was inferred.",
    ]
    bundle["reactions"].append(r)

ensure_molecule("m_vinyl_ether_scope", "vinyl ether substrate (identity not specified in prose)", [P[2]["t"]],
                template=generic_template("Vinyl ether substrate; substituents are not specified in the page-2 prose", ["R"], [P[2]["t"]]),
                chirality="unresolved", geometry="unresolved", representation="not_encoded")
ensure_molecule("m_vinyl_ether_product", "desired vinyl-ether difunctionalization product", [P[2]["t"]],
                chirality="unresolved", geometry="unresolved", representation="not_encoded",
                identity_note="The product is reported only qualitatively and is not assigned a compound number on pages 1–5.")
r = reaction("scope_vinyl_ether_low_yield", [P[2]["t"], P[2]["i"]], categories=("substrate_scope",))
set_location(r, "Page 2 scope discussion following Figure 2", 2, "prose adjacent to Figure 2")
set_fact_name(r, "1,2-alkynylalkenylation of vinyl ether")
r["general_scheme_ids"] = [GENERIC_ID]
set_common_photo(r, r["evidence_ids"], "m_1a", "m_vinyl_ether_scope", "m_vinyl_ether_product",
                 substrate_eq=1.5, catalyst_pct=1, kh_eq=2, solvents=("MeCN", "MeOH"), volumes=(0.5, 1.5), power=5, time_h=None)
add_outcome(r, "yield", "yield", ["m_vinyl_ether_product"], r["evidence_ids"], "<10%",
            {"number": 10, "unit": "%", "qualifier": "lt", "raw_text": "<10%"}, product_role="desired product")
classification(r, "alkynylalkenylation")
r["review"]["unresolved"] = ["Specific vinyl ether and desired-product identities are not stated in the supplied prose.",
                                 "Reaction time and yield measurement basis are unreported on pages 1–5.",
                                 "Complete procedure and characterization require unavailable SI."]
bundle["reactions"].append(r)


# Figure 3 reagent assessment.
fig3 = [
    ("5a", "1b", 60, 30, None), ("5b", "1c", 57, 30, None),
    ("5c", "1d", 53, 36, None), ("5d", "1e", 71, 36, "dr = 1:1"),
    ("5e", "1f", 33, 48, "dr = 1:1"), ("5f", "1g", 34, 48, None),
    ("5g", "1h", 65, 37, None), ("5h", "1i", 78, 36, None),
    ("5i", "1j", 72, 36, None), ("5j", "1k", 63, 34, None),
    ("5k", "1l", 55, 28, None), ("5l", "1m", 23, 18, "E"),
    ("5m", "1n", 47, 72, "E,E"),
]
for product_label, reagent_label, yld, hours, sel in fig3:
    ev = [P[3]["i"], P[3]["t"]]
    reagent_mid = f"m_{reagent_label}"
    ensure_molecule(reagent_mid, reagent_label, ev, original_name=f"dual-function reagent {reagent_label}",
                    chirality="unresolved", geometry="unspecified", representation="not_encoded",
                    axial={"allene_axis": "unassigned; compact drawing inspected but no source descriptor"},
                    identity_note="Complete candidate SMILES not asserted from the compact Figure 3 drawing.")
    product_mid = f"m_product_{product_label}"
    if sel in ("E", "E,E"):
        skw = dict(chirality="unspecified", geometry="specified",
                   double_bonds={"reported_geometry": sel},
                   assignment={"geometry": sel, "basis": "Figure 3 label/prose", "absolute_chirality": "unassigned"})
    else:
        skw = stereo_for_product(product_label, sel)
    ensure_molecule(product_mid, product_label, ev, original_name=f"product {product_label}",
                    representation="not_encoded", identity_note="Source drawing retained; complete candidate SMILES not asserted.", **skw)
    r = reaction(f"scope_{product_label}", ev, categories=("substrate_scope",))
    set_location(r, f"Figure 3, {reagent_label} to {product_label}", 3, "Figure 3")
    set_fact_name(r, "1,2-alkynylalkenylation: assessment of dual-function reagents")
    r["general_scheme_ids"] = [GENERIC_ID]
    set_common_photo(r, ev, reagent_mid, "m_2a", product_mid, substrate_eq=1.5, catalyst_pct=1,
                     kh_eq=2, solvents=("MeCN", "MeOH"), volumes=(0.5, 1.5), power=5, time_h=hours)
    add_outcome(r, "yield", "yield", [product_mid], ev, f"{yld}%", {"number": yld, "unit": "%"}, product_role="main product")
    if sel:
        metric = "double-bond geometry" if sel in ("E", "E,E") else "diastereomeric ratio"
        add_outcome(r, "selectivity", "selectivity", [product_mid], ev, sel if sel in ("E", "E,E") else sel, sel,
                    scope={"metric": metric, "order": "source-reported", "major_configuration": sel if sel in ("E", "E,E") else "unassigned"})
    classification(r, "alkynylalkenylation")
    r["notes"] = [f("Yield measurement basis is not stated in Figure 3.", "Yield measurement basis is not stated in Figure 3.", ev)]
    r["review"]["unresolved"] = ["Yield/selectivity measurement methods are unreported on pages 1–5.",
                                     "Complete structures and procedures remain dependent on unavailable SI.",
                                     "Exact atom mapping is unresolved."]
    bundle["reactions"].append(r)


# The prose reports plural vinylcyclopropane radical-clock experiments but
# defers structures and individual results to SI.  Retain one unresolved family
# record rather than inventing a count or assigning products to unknown VCPs.
ensure_molecule("m_vcp_family", "vinylcyclopropanes (VCPs)", [P[3]["t"]],
                template=generic_template("Vinylcyclopropane radical-clock family", ["electronic substituents"], [P[3]["t"]]),
                chirality="unresolved", geometry="unresolved", representation="not_encoded")
r = reaction("vcp_radical_clock_family", [P[3]["t"], P[3]["i"]], record_type="unresolved",
             categories=("mechanistic_study",))
set_location(r, "Page 3 radical-clock prose; details referred to SI", 3, "prose adjacent to Figure 2 continuation")
set_fact_name(r, "Vinylcyclopropane radical-clock experiments", "radical-clock study")
r["participants"] = [participant("p_1a", "m_1a", "reactant", r["evidence_ids"]),
                     participant("p_vcp", "m_vcp_family", "reactant", r["evidence_ids"])]
add_step(r, "reported_family", r["evidence_ids"], ["p_1a", "p_vcp"],
         operation="VCP radical-clock experiments; individual structures and conditions are in unavailable SI")
add_outcome(r, "reported_outcome_family", "other", [], r["evidence_ids"],
            "either the ring-opening product or the difunctionalized product depending on the electronic property of VCPs",
            "ring-opening versus difunctionalization depended on VCP electronics",
            scope={"individual_experiment_count": "unreported on pages 1–5", "product_identities": "unreported"})
classification(r, "unresolved")
r["review"]["unresolved"] = ["Plural VCP experiments are reported, but their number, structures, conditions, and individual outcomes are unavailable without SI.",
                                 "The family was not split into fabricated per-substrate records."]
bundle["reactions"].append(r)


# Figure 4 enynylalkenylation scope.
ensure_molecule("m_6a", "6a", [P[3]["i"]], original_name="product 6a",
                smiles="C=C(C)C#CC(c1ccccc1)CC=Cc1ccccc1",
                chirality="unspecified", geometry="unspecified", representation="partial",
                identity_note="Candidate graph contains unassigned tetrahedral and alkene stereochemistry.")
ensure_molecule("m_6b", "6b", [P[3]["i"], P[4]["i"]], original_name="product 6b",
                smiles="C=C(C)C#CC(c1ccc(OC)cc1)CC=Cc1ccccc1",
                chirality="unspecified", geometry="unspecified", representation="partial",
                identity_note="Candidate graph transcribed from Figure 4/Scheme 2; configuration is unassigned.")
co_substrates = {
    "6a": ("m_2a", None),
    "6b": ("m_p_methoxystyrene", "C=Cc1ccc(OC)cc1"),
    "6c": ("m_p_cyanostyrene", "C=Cc1ccc(C#N)cc1"),
    "6d": ("m_p_fluorostyrene", "C=Cc1ccc(F)cc1"),
    "6e": ("m_enyne_for_6e", None),
    "6f": ("m_diene_for_6f", None),
    "6g": ("m_ethyl_acrylate", "C=CC(=O)OCC"),
    "6h": ("m_2a", None), "6i": ("m_2a", None), "6j": ("m_2a", None),
    "6k": ("m_1_hexene", "C=CCCCC"),
}
reagents6 = {"6h": "m_1b", "6i": "m_1c", "6j": "m_1e", "6k": "m_1n"}
fig4 = [
    ("6a", 61, 9, None, None), ("6b", 53, 22, None, None),
    ("6c", 52, 20, None, None), ("6d", 54, 19, None, None),
    ("6e", 35, 20, None, None), ("6f", 38, 20, None, "a"),
    ("6g", 31, 9, None, None), ("6h", 68, 20, None, None),
    ("6i", 69, 20, None, None), ("6j", 50, 20, "internal:terminal = 2.7:1", None),
    ("6k", 35, 17, None, None),
]
for label, yld, hours, sel, foot in fig4:
    ev = [P[3]["i"], P[3]["t"], P[4]["t"]]
    substrate_mid, smi = co_substrates[label]
    if substrate_mid != "m_2a":
        ensure_molecule(substrate_mid, f"alkene substrate for {label}", [P[3]["i"]],
                        original_name=f"alkene substrate drawn for {label}", smiles=smi,
                        chirality="achiral" if smi else "unresolved",
                        geometry="not_applicable" if smi and smi.startswith("C=C") else "unresolved",
                        representation="encoded" if smi else "not_encoded",
                        identity_note="Candidate structure is source-derived; compact drawing remains machine-unverified." if smi else "Complete candidate SMILES not asserted.")
    reagent_mid = reagents6.get(label, "m_1a")
    product_mid = "m_6a" if label == "6a" else "m_6b" if label == "6b" else f"m_product_{label}"
    if label not in ("6a", "6b"):
        ensure_molecule(product_mid, label, [P[3]["i"]], original_name=f"product {label}",
                        chirality="unspecified", geometry="mixture" if label == "6j" else "unspecified",
                        mixture={"type": "internal/terminal constitutional-isomer mixture", "ratio": "2.7:1"} if label == "6j" else None,
                        representation="not_encoded", identity_note="Complete candidate SMILES not asserted from compact Figure 4 drawing.")
    r = reaction(f"scope_{label}", ev, categories=("substrate_scope",))
    set_location(r, f"Figure 4, product {label}", 3, "Figure 4")
    set_fact_name(r, "1,2-enynylalkenylation of alkenes")
    r["general_scheme_ids"] = [GENERIC_ID]
    parts = [
        participant("p_dfr", reagent_mid, "reactant", ev, amount=("0.2 mmol", {"number": 0.2, "unit": "mmol"})),
        participant("p_alkene", substrate_mid, "reactant", ev,
                    equivalents=("1.5 equiv" if foot == "a" else "3 equiv",
                                 {"number": 1.5 if foot == "a" else 3, "unit": "equiv"}), reference=reagent_mid),
        participant("p_photocatalyst", "fac_irppy3", "catalyst", ev,
                    mol_percent=("2 mol%", {"number": 2, "unit": "mol%"}), reference=reagent_mid),
        participant("p_base", "m_kh2po4", "additive", ev,
                    equivalents=("1 equiv", {"number": 1, "unit": "equiv"}), reference=reagent_mid),
        participant("p_solvent", "m_dmf", "solvent", ev, volume=("2 mL", {"number": 2, "unit": "mL"})),
        participant("p_product", product_mid, "product", ev, product_role="main product"),
    ]
    r["participants"] = parts
    r["reaction_scale"] = f("1 (0.2 mmol)", {"reference_molecule_id": reagent_mid, "amount": {"number": 0.2, "unit": "mmol"}}, ev)
    r["solvent_conditions"]["composition"] = f("DMF", ["DMF"], ev)
    r["solvent_conditions"]["total_volume"] = f("2 mL", {"number": 2, "unit": "mL"}, ev)
    add_step(r, "reaction", ev, [p["participant_id"] for p in parts if p["role"] != "product"],
             operation="irradiated with 30 W blue LEDs", temperature=("rt", "rt"),
             time=(f"{hours} h", {"number": hours, "unit": "h"}), atmosphere="N2")
    r["special_conditions"]["light_source"] = f("blue LEDs", "blue LEDs", ev)
    r["special_conditions"]["power"] = f("30 W", {"number": 30, "unit": "W"}, ev)
    r["procedure"]["original_text"] = f("Figure 4 reaction conditions", "Figure 4 reaction conditions", ev)
    add_outcome(r, "yield", "yield", [product_mid], ev, f"{yld}%", {"number": yld, "unit": "%"}, product_role="main product")
    if sel:
        add_outcome(r, "selectivity", "selectivity", [product_mid], ev, sel, sel,
                    scope={"metric": "internal:terminal product ratio", "order": ["internal", "terminal"], "major": "internal"})
    classification(r, "enynylalkenylation")
    r["notes"] = [f("Yield/selectivity measurement basis is not stated in Figure 4.",
                    "Yield/selectivity measurement basis is not stated in Figure 4.", ev),
                  f("Adjacent page-4 prose reports over-functionalization for terminal-enynyl products as a reason for lower yields.",
                    "Adjacent page-4 prose reports over-functionalization for terminal-enynyl products as a reason for lower yields.", [P[4]["t"]])]
    r["review"]["unresolved"] = ["Yield/selectivity measurement methods are unreported.",
                                     "Complete structures/procedures remain unavailable without SI.",
                                     "Exact atom mapping is unresolved."]
    bundle["reactions"].append(r)


# Scheme 2 product transformations.
for label in ("7", "8", "9", "10", "11", "12"):
    kwargs = dict(chirality="unspecified", geometry="unspecified", representation="not_encoded",
                  identity_note="Product structure retained by Scheme 2 drawing; complete candidate SMILES not asserted.")
    if label == "8":
        kwargs.update(geometry="specified", double_bonds={"reported_geometry": "Z,E"},
                      assignment={"double_bond_geometry": "Z,E"})
    elif label == "10":
        kwargs.update(geometry="specified", double_bonds={"reported_geometry": "Z,E"},
                      assignment={"double_bond_geometry": "Z,E"})
    elif label == "9":
        kwargs.update(smiles="CC(C)C#CC(c1ccccc1)CC=Cc1ccccc1", representation="partial",
                      atom_centers=[{"candidate_graph": "tetrahedral center adjacent to the alkyne", "configuration": "unassigned"}],
                      identity_note="Candidate graph transcribed from Scheme 2; the tetrahedral center and surviving alkene geometry are unassigned and the structure remains machine-unverified.")
    elif label == "11":
        kwargs.update(axial={"allene_axis": "unassigned"}, chirality="unspecified")
    elif label == "12":
        kwargs.update(chirality="mixture", mixture={"type": "diastereomeric mixture", "ratio": "1:1"})
    ensure_molecule(f"m_product_{label}", label, [P[4]["i"]], original_name=f"product {label}", **kwargs)

scheme2_specs = [
    ("scheme2_epoxide_7", "m_3a", "m_product_7", 73, "epoxidation", "a"),
    ("scheme2_lindlar_8", "m_3a", "m_product_8", 68, "hydrogenation", "b"),
    ("scheme2_deoxygenation_9", "m_3a", "m_product_9", 60, "deoxygenation", "c"),
    ("scheme2_reduction_10", "m_product_9", "m_product_10", 78, "hydrogenation", "d"),
    ("scheme2_allene_11", "m_3a", "m_product_11", 77, "allene_formation", "e"),
    ("scheme2_cyanoalkylation_12", "m_6b", "m_product_12", 54, "cyanoalkylation", "f"),
]
conditions = {
    "a": dict(reagents=[("m_mcpba", "m-CPBA", "3 equiv", 3, "reagent"),
                         ("m_nahco3", "NaHCO3", "3 equiv", 3, "additive")],
              solvents=[("m_dcm", "DCM", None), ("m_water", "H2O", None)], temp=None, time=("20 h", {"number": 20, "unit": "h"}), pressure=None),
    "b": dict(reagents=[("m_lindlar", "Lindlar's catalyst", "20 mol%", 20, "catalyst"),
                         ("m_quinoline", "quinoline", "0.5 equiv", 0.5, "additive"),
                         ("m_h2", "H2", None, None, "reagent")],
              solvents=[("m_etoac", "dry EtOAc", None)], temp=("rt", "rt"), time=("15 min", {"number": 15, "unit": "min"}),
              pressure=("1 atm H2", "1 atm H2"), extra="Lindlar's catalyst also reported as 5 wt %"),
    "c": dict(reagents=[("m_cp2zrcl2", "Cp2ZrCl2", "1.6 equiv", 1.6, "reagent"),
                         ("m_n_buli", "n-BuLi (2.5 M)", "3.2 equiv", 3.2, "reagent")],
              solvents=[("m_thf", "THF", None)], temp=None, time=None, pressure=None),
    "d": dict(reagents=[("m_dibal", "DIBAL-H", "2 equiv", 2, "reagent")], solvents=[("m_toluene", "toluene", None)],
              temp=("50 °C", {"number": 50, "unit": "°C"}), time=("12 h", {"number": 12, "unit": "h"}), pressure=None),
    "e": dict(reagents=[("m_dibal", "DIBAL-H", "4 equiv", 4, "reagent")], solvents=[("m_toluene", "toluene", None)],
              temp=("80 °C", {"number": 80, "unit": "°C"}), time=("21 h", {"number": 21, "unit": "h"}), pressure=None),
    "f": dict(reagents=[("m_cu_oac2", "Cu(OAc)2", "5 mol%", 5, "catalyst"),
                         ("m_phen", "1,10-phenanthroline", "7 mol%", 7, "ligand"),
                         ("m_lpo", "LPO", "3 equiv", 3, "reagent"),
                         ("m_tmscn", "TMSCN", "4 equiv", 4, "reagent")],
              solvents=[("m_dcm", "dry DCM", None)], temp=("50 °C", {"number": 50, "unit": "°C"}),
              time=("42 h", {"number": 42, "unit": "h"}), pressure=None),
}
for rid, start, prod, yld, ckind, key in scheme2_specs:
    ev = [P[4]["i"], P[4]["t"]]
    c = conditions[key]
    r = reaction(rid, ev, categories=("derivatization",))
    set_location(r, f"Scheme 2, condition {key}", 4, "Scheme 2")
    names = {"epoxidation": "epoxidation of 3a", "hydrogenation": "selective reduction",
             "deoxygenation": "removal of propargylic methoxy group", "allene_formation": "formal allenylalkenylation",
             "cyanoalkylation": "copper-catalyzed radical cyanoalkylation"}
    set_fact_name(r, names[ckind], ckind, inferred=ckind in ("deoxygenation", "allene_formation"))
    ensure_molecule(start, start.replace("m_", ""), ev) if start not in molecules else None
    parts = [participant("p_start", start, "reactant", ev), participant("p_product", prod, "product", ev, product_role="main product")]
    for idx, (mid, label, rawq, number, role) in enumerate(c["reagents"], 1):
        ensure_molecule(mid, label, ev, original_name=label)
        kw = {}
        if rawq:
            if "mol%" in rawq:
                kw = {"mol_percent": (rawq, {"number": number, "unit": "mol%"}), "reference": start}
            elif "equiv" in rawq:
                kw = {"equivalents": (rawq, {"number": number, "unit": "equiv"}), "reference": start}
        parts.append(participant(f"p_reagent_{idx}", mid, role, ev, **kw))
    for idx, (mid, label, _) in enumerate(c["solvents"], 1):
        ensure_molecule(mid, label, ev, original_name=label)
        parts.append(participant(f"p_solvent_{idx}", mid, "solvent", ev))
    r["participants"] = parts
    r["solvent_conditions"]["composition"] = f("/".join(x[1] for x in c["solvents"]), [x[1] for x in c["solvents"]], ev)
    if key == "c":
        add_step(r, "low_temperature_hold", ev, [p["participant_id"] for p in parts if p["role"] != "product"],
                 operation="reaction at -78 °C", temperature=("-78 °C", {"number": -78, "unit": "°C"}),
                 time=("1 h", {"number": 1, "unit": "h"}))
        add_step(r, "room_temperature_hold", ev, [p["participant_id"] for p in parts if p["role"] != "product"],
                 operation="then held at rt; warming duration not stated", temperature=("rt", "rt"),
                 time=("3 h", {"number": 3, "unit": "h"}))
    else:
        add_step(r, "reaction", ev, [p["participant_id"] for p in parts if p["role"] != "product"],
                 operation="Scheme 2 condition " + key, temperature=c["temp"], time=c["time"], pressure=c["pressure"])
    if key == "b":
        # Preserve both catalyst specifications as reported without pretending wt% and mol% are interchangeable.
        r["notes"].append(f("Lindlar's catalyst is reported as 5 wt % and 20 mol %; 5 wt % describes the catalyst preparation/grade rather than a mass charge.",
                            "Lindlar's catalyst: 5 wt %, 20 mol %", ev))
    if key == "c":
        r["notes"].append(f("n-BuLi concentration 2.5 M is reported, but solution volume is not stated.",
                            "n-BuLi concentration 2.5 M is reported, but solution volume is not stated.", ev))
    r["procedure"]["original_text"] = f("Scheme 2 condition " + key, "Scheme 2 condition " + key, ev)
    add_outcome(r, "yield", "yield", [prod], ev, f"{yld}%", {"number": yld, "unit": "%"}, product_role="main product")
    if prod == "m_product_12":
        add_outcome(r, "selectivity", "selectivity", [prod], ev, "dr = 1:1", "dr = 1:1",
                    scope={"metric": "diastereomeric ratio", "order": "unreported", "major_configuration": "none assigned"})
    classification(r, ckind)
    r["review"]["unresolved"] = ["Yield/selectivity measurement method is unreported in Scheme 2.",
                                     "Reaction scale, workup, and purification are unreported on pages 1–5.",
                                     "Exact atom mapping and complete product stereochemistry are unresolved."]
    bundle["reactions"].append(r)


# Main-text mechanistic experiments explicitly reported; SI-only details remain missing.
ensure_molecule("m_tempo", "TEMPO", [P[4]["t"]], original_name="TEMPO")
mech_specs = [
    ("stern_volmer_quenching", "Stern–Volmer quenching study", ["m_1a", "fac_irppy3"],
     "excited-state Ir(ppy)3 was readily oxidatively quenched by 1a", "Figure S1"),
    ("reduction_potential_1a", "reduction-potential measurement", ["m_1a"],
     "Ep/2 = -1.55 V vs. SCE", "Figure S2"),
    ("quantum_yield_1a_styrene", "reaction quantum-yield measurement", ["m_1a", "m_2a", "m_3a"],
     "quantum yield = 0.14", "Figure S3"),
    ("light_on_off", "light on-off experiment", ["m_1a", "m_2a", "m_3a"],
     "reaction outcome was not increased without light irradiation", "Figure S4"),
    ("tempo_trapping", "TEMPO radical-trapping experiment", ["m_1a", "m_tempo"],
     "TEMPO trapping supports formation of allenyl radical I; details are in unavailable SI", "SI details not supplied"),
]
for rid, name, mids, observation, source_ref in mech_specs:
    ev = [P[4]["t"]]
    r = reaction(rid, ev, categories=("mechanistic_study",))
    set_location(r, f"Page 4 prose; {source_ref}", 4, source_ref)
    set_fact_name(r, name, name)
    parts = []
    for idx, mid in enumerate(mids, 1):
        role = "product" if mid == "m_3a" else "reactant"
        parts.append(participant(f"p_{idx}", mid, role, ev, product_role="monitored product" if role == "product" else None))
    r["participants"] = parts
    add_step(r, "measurement", ev, [p["participant_id"] for p in parts if p["role"] != "product"],
             operation=name + "; detailed SI conditions unavailable")
    add_outcome(r, "observation", "other", ["m_3a"] if "m_3a" in mids else [], ev, observation, observation,
                scope={"source_reference": source_ref})
    classification(r, "not_applicable" if rid in ("stern_volmer_quenching", "reduction_potential_1a") else "unresolved")
    r["review"]["unresolved"] = ["The paper cites SI figures/details, but no SI file was supplied.",
                                     "Exact apparatus, quantities, timing, and measurement protocol are unavailable on pages 1–5."]
    bundle["reactions"].append(r)


# Proposed mechanism is one record, not executed intermediates.
for mid, label in [("m_intermediate_i", "intermediate I"), ("m_intermediate_ii", "intermediate II"),
                   ("m_intermediate_iii", "intermediate III"), ("m_intermediate_iv", "intermediate IV"),
                   ("m_intermediate_v", "intermediate V"), ("m_intermediate_vi", "intermediate VI")]:
    ensure_molecule(mid, label, [P[4]["i"]], original_name=label, chirality="unresolved", geometry="unresolved",
                    representation="not_encoded", identity_note="Proposed mechanistic intermediate; not an isolated molecule record.")
r = reaction("figure5_proposed_mechanism", [P[4]["i"], P[4]["t"]], record_type="proposed_mechanism")
set_location(r, "Figure 5", 4, "Figure 5")
set_fact_name(r, "Proposed reaction mechanism", "author-proposed photocatalytic radical/cation pathway")
r["participants"] = [participant("p_1a", "m_1a", "reactant", r["evidence_ids"]),
                     participant("p_alkene", generic_alkene, "reactant", r["evidence_ids"]),
                     participant("p_product", generic_product, "product", r["evidence_ids"], product_role="generic product branch")]
add_step(r, "proposed_sequence", r["evidence_ids"], ["p_1a", "p_alkene"],
         operation="SET, radical addition, 5-exo-trig migration, SO2 extrusion, tautomerization, oxidation, then nucleophile trapping or deprotonation")
classification(r, "unresolved")
r["notes"] = [f("Intermediates I–VI and the two terminal branches are author-proposed and were not expanded into executed experiments.",
                "Intermediates I–VI and the two terminal branches are author-proposed and were not expanded into executed experiments.", r["evidence_ids"])]
r["review"]["unresolved"] = ["Proposed pathway is not experimental proof of each elementary step.",
                                 "Exact atom mapping and intermediate structures are not asserted."]
bundle["reactions"].append(r)


# Non-synthetic observations and negative inventory.
bundle["source_context"] = [
    {"context_id": f"{DOC}:metadata", "document_id": DOC, "kind": "document_metadata",
     "evidence_ids": [P[1]["i"], P[1]["t"]],
     "title": "Alkene Difunctionalization Triggered by a Stabilized Allenyl Radical: Concomitant Installation of Two Unsaturated C–C Bonds",
     "doi": "10.1002/anie.202106145"},
    {"context_id": f"{DOC}:mechanistic_photophysics", "document_id": DOC, "kind": "mechanistic_context",
     "evidence_ids": [P[4]["t"]],
     "observations": ["Ir(ppy)3 excited-state reduction potential quoted as -1.73 V vs SCE from reference 15.",
                      "Authors interpret quantum yield 0.14 and light on-off behavior as supporting a photocatalytic cycle."]},
    {"context_id": f"{DOC}:no_standalone_spectra", "document_id": DOC, "kind": "negative_inventory",
     "evidence_ids": [P[n]["i"] for n in range(1, 6)],
     "finding": "No standalone mass-spectrum peak table, molecular-weight table, or fragment-relative-abundance block occurs on pages 1–5."},
    {"context_id": f"{DOC}:si_unavailable", "document_id": DOC, "kind": "source_limitation",
     "evidence_ids": [P[1]["i"]],
     "finding": "The paper advertises online Supporting Information, but no SI file was explicitly supplied or associated in this workspace."},
]


bundle["molecules"] = list(molecules.values())

# Page-bounded raw inventory and required coverage object.
page_records = {
    1: [f"{DOC}:scheme1a_prior_cross_coupling", f"{DOC}:scheme1b_prior_migration", GENERIC_ID],
    2: [GENERIC_ID] + [f"{DOC}:scope_{x[0]}" for x in scope2]
       + [f"{DOC}:{x[0]}" for x in page2_prose_specs]
       + [f"{DOC}:scope_vinyl_ether_low_yield"],
    3: [GENERIC_ID] + [f"{DOC}:scope_{x[0]}" for x in scope2 if x[0] in {"3ab","3ac","3ad","3ae","3af","3ag","3ah","3ai","3aj","4a","4b","4c","4d","4e","4f"}]
       + [f"{DOC}:scope_{x[0]}" for x in fig3] + [f"{DOC}:vcp_radical_clock_family"]
       + [f"{DOC}:scope_{x[0]}" for x in fig4],
    4: [f"{DOC}:scope_{x[0]}" for x in fig4]
       + [f"{DOC}:{x[0]}" for x in scheme2_specs]
       + [f"{DOC}:{x[0]}" for x in mech_specs]
       + [f"{DOC}:figure5_proposed_mechanism"],
    5: [],
}
regions = {
    1: [
        {"region": "title/abstract and bibliographic header", "finding": "Current paper metadata and broad transformation claim; no concrete current experiment.", "evidence_ids": [P[1]["i"], P[1]["t"]]},
        {"region": "Scheme 1A", "finding": "Two generic cited prior-work reaction families, retained as literature examples rather than current experiments.", "evidence_ids": [P[1]["i"]]},
        {"region": "Scheme 1B", "finding": "Present-work generic dual-function-reagent docking/migration concept; generic scheme only.", "evidence_ids": [P[1]["i"], P[1]["t"]]},
    ],
    2: [
        {"region": "Figure 1", "finding": "Optimized 1a + styrene experiment giving 3a, 72% after 47 h; merged with the same 3a scope identity.", "evidence_ids": [P[2]["i"], P[2]["t"]]},
        {"region": "Figure 2 main grid", "finding": "Product outcomes 3b–3aj with printed yields, times, cis/dr annotations, and footnotes a/b.", "evidence_ids": [P[2]["i"]]},
        {"region": "Figure 2 co-solvent panel and caption", "finding": "Products 4a–4f plus common conditions and footnotes a–c; co-solvent volumes retained.", "evidence_ids": [P[2]["i"], P[2]["t"]]},
        {"region": "optimization/scope prose around Figures 1–2", "finding": "Two methanol-condition observations, four ineffective alternative-base experiments, and an unnumbered vinyl-ether result below 10%; retained independently without inventing SI details.", "evidence_ids": [P[2]["i"], P[2]["t"]]},
    ],
    3: [
        {"region": "left-column continuation of Figure 2 discussion", "finding": "Source discussion for 3ab–3aj and 4a–4f, including cis outcomes, natural-product derivatization, and an SI-only plural VCP radical-clock series retained as unresolved.", "evidence_ids": [P[3]["i"], P[3]["t"]]},
        {"region": "Figure 3", "finding": "Thirteen distinct dual-function reagent scope experiments 5a–5m with yields/times and dr/E/E,E annotations.", "evidence_ids": [P[3]["i"], P[3]["t"]]},
        {"region": "Figure 4", "finding": "Eleven distinct enynylalkenylation experiments 6a–6k with yields/times, one footnote-a loading change, and 6j internal:terminal ratio.", "evidence_ids": [P[3]["i"], P[3]["t"]]},
    ],
    4: [
        {"region": "top-left continuation of Figure 4 discussion", "finding": "Discussion of 6a–6k, including over-functionalization and internal-enyne cases; linked to page-3 scope records.", "evidence_ids": [P[4]["i"], P[4]["t"]]},
        {"region": "Scheme 2", "finding": "Six distinct product-transformation arrows yielding 7–12; each condition letter is a separate experiment.", "evidence_ids": [P[4]["i"], P[4]["t"]]},
        {"region": "mechanistic-study prose", "finding": "Stern–Volmer, reduction potential, quantum yield, light on-off, and TEMPO-trapping experiments; detailed protocols refer to unavailable SI.", "evidence_ids": [P[4]["t"]]},
        {"region": "Figure 5", "finding": "Author-proposed mechanism, retained as one proposed-mechanism record rather than executed intermediate experiments.", "evidence_ids": [P[4]["i"], P[4]["t"]]},
    ],
    5: [
        {"region": "summary/acknowledgements/conflict statement", "finding": "No new synthesis experiment, table row, or reaction arrow.", "evidence_ids": [P[5]["i"], P[5]["t"]]},
        {"region": "references and manuscript dates", "finding": "Bibliography only; cited works were not imported as current-paper experiments.", "evidence_ids": [P[5]["i"], P[5]["t"]]},
    ],
}
coverage = {"document_id": DOC, "semantic_accuracy_verified": False, "pages": []}
for page in range(1, 6):
    coverage["pages"].append({
        "page": page,
        "state": "attempted",
        "regions": regions[page],
        "record_ids": list(dict.fromkeys(page_records[page])),
        "unresolved": (["Supporting-information procedures/characterization are not available and were not inferred."]
                       if page in (1, 2, 3, 4) else
                       ["References were inventoried only; cited-paper experiments were not source-extracted from this PDF."]),
    })

raw_inventory = {
    "document_id": DOC,
    "pages_attempted": [1, 2, 3, 4, 5],
    "source_images": [f"assets/{DOC}/page{n}.jpg" for n in range(1, 6)],
    "counts": {
        "literature_examples": 2,
        "generic_schemes": 1,
        "figure1_2_specific_experiments": len(scope2),
        "page2_prose_optimization_experiments": len(page2_prose_specs),
        "page2_prose_scope_experiments": 1,
        "figure3_specific_experiments": len(fig3),
        "vcp_unresolved_experiment_family": 1,
        "figure4_specific_experiments": len(fig4),
        "scheme2_derivatizations": len(scheme2_specs),
        "mechanistic_experiments": len(mech_specs),
        "proposed_mechanisms": 1,
    },
    "deduplication": ["Figure 1 optimized 3a and the 3a scope identity are represented by one reaction record with merged evidence."],
    "page_regions": regions,
    "limitations": ["No SI supplied.", "Compact drawings not confidently transcribed remain unresolved rather than forced to SMILES.",
                    "Page enumeration and schema validity do not establish semantic or chemical accuracy."],
}

mentions = st.source_mentions(json.loads((RUN / "dataset.json").read_text()), ["1a", "2a", "3a", "6b", "TEMPO"], context=180)

(RUN / "submitted.json").write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")
(RAW / "page_extraction.json").write_text(json.dumps(raw_inventory, ensure_ascii=False, indent=2) + "\n")
(RAW / "source_mentions.json").write_text(json.dumps(mentions, ensure_ascii=False, indent=2) + "\n")
(RAW / "stereo_checks.json").write_text(json.dumps(stereo_diagnostics, ensure_ascii=False, indent=2) + "\n")
covdir = RUN / "coverage"
covdir.mkdir(parents=True, exist_ok=True)
(covdir / "pages_0001_0005.json").write_text(json.dumps(coverage, ensure_ascii=False, indent=2) + "\n")

print(json.dumps({"molecules": len(bundle["molecules"]), "reactions": len(bundle["reactions"]),
                  "source_context": len(bundle["source_context"]),
                  "coverage_pages": [x["page"] for x in coverage["pages"]]}, ensure_ascii=False))
