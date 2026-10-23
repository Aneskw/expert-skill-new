#!/usr/bin/env python3
"""Build the source-derived v3.1 machine submission for the six-page paper.

This file is retained as raw model work. It deliberately keeps complex drawn
structures unresolved and all imported records at needs_review.
"""
from __future__ import annotations

import copy
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

RUN = ROOT / "runs/current"
DOC = "doc_9fe9c91a39ef967d"
E = {n: f"{DOC}:page{n}:image" for n in range(1, 7)}
T = {n: f"{DOC}:page{n}:text" for n in range(1, 7)}
_MISSING = object()


def f(raw=None, normalized=_MISSING, evidence=(), method="original_extraction",
      status=None, derivation=None, alternatives=()):
    if normalized is _MISSING:
        normalized = raw
    return st.fact(raw, normalized, evidence, method=method, status=status,
                   derivation=derivation, alternatives=alternatives)


def inferred(normalized, evidence, derivation, raw=None):
    return f(raw, normalized, evidence, method="inference", derivation=derivation)


def derived(normalized, evidence, derivation, raw=None):
    return f(raw, normalized, evidence, method="structure_derivation", derivation=derivation)


bundle = json.loads((RUN / "submitted.json").read_text())
bundle["molecules"] = []
bundle["reactions"] = []
bundle["source_context"] = []
molecules = {}

SMILES = {
    "1a": "Cc1ccc(N=C2C=CC=CC=C2)cc1",
    "2a": "C=C=CN(C)S(=O)(=O)c1ccc(C)cc1",
    "3a": "CN(C#Cc1ccccc1)S(=O)(=O)c1ccc(C)cc1",
    "DCM": "ClCCl",
    "DCE": "ClCCCl",
    "toluene": "Cc1ccccc1",
    "hexane": "CCCCCC",
    "THF": "C1CCOC1",
    "acetonitrile": "CC#N",
    "triphenylphosphine": "P(c1ccccc1)(c1ccccc1)c1ccccc1",
    "ethyl acetate": "CCOC(C)=O",
}


def mid(label):
    slug = re.sub(r"[^A-Za-z0-9]+", "-", label).strip("-").lower()
    return f"{DOC}:mol:{slug}"


def ensure_molecule(label, evidence, name=None, note=None, source_label=None):
    evidence = list(dict.fromkeys(evidence))
    if label in molecules:
        m = molecules[label]
        occ = m["occurrences"][0]
        occ["evidence_ids"] = list(dict.fromkeys(occ["evidence_ids"] + evidence))
        for key in ("compound_label",):
            occ[key]["evidence_ids"] = list(dict.fromkeys(occ[key]["evidence_ids"] + evidence))
        return m["molecule_id"]
    m = st.new_molecule(mid(label), DOC, evidence)
    shown = source_label if source_label is not None else label
    m["occurrences"][0]["compound_label"] = f(shown, shown, evidence)
    m["names"]["original"] = f(name or shown, name or shown, evidence)
    m["names"]["english"] = f(name or shown, name or shown, evidence)
    m["names"]["abbreviations"] = f([label], [label], evidence)
    m["names"]["aliases"] = f()
    m["cas_numbers"] = f()
    m["external_properties"] = []
    if label in SMILES:
        smi = SMILES[label]
        m["representation_type"] = "smiles"
        m["smiles"] = derived(
            smi, evidence,
            "Candidate connectivity transcribed from the source drawing/name; RDKit parsing and stereo diagnostics do not establish source correctness.",
            raw={"source_label": shown, "candidate_from_drawing": True},
        )
        m["identity_status"] = "source_identified"
        if label in {"1a", "2a", "3a"}:
            statuses = ("unresolved", "unresolved", "partial")
            stereo_reason = "Complex substrate drawing inspected; absolute/axial/E-Z assignments are not stated and the candidate SMILES is not treated as a verified stereochemical identity."
        else:
            statuses = ("achiral", "not_applicable", "not_applicable")
            stereo_reason = "No potential stereogenic atom/axis or E/Z element is present in the candidate small-molecule graph; this is a structure-derived assessment."
        for key, val in zip(("chirality_status", "geometry_status", "representation_status"), statuses):
            m["stereochemistry"][key] = inferred(val, evidence, stereo_reason, raw="stereo review of source drawing and candidate graph")
        m["identity_note"] = f(
            note or "Concrete candidate structure retained for machine review; no CAS was supplied.",
            note or "Concrete candidate structure retained for machine review; no CAS was supplied.", evidence)
    else:
        m["representation_type"] = "unresolved"
        m["identity_status"] = "source_identified"
        why = note or "The source drawing/label was retained, but no complete source-verified SMILES was assigned in this extraction."
        m["identity_note"] = f(why, why, evidence)
        for key in ("chirality_status", "geometry_status", "representation_status"):
            m["stereochemistry"][key] = inferred(
                "unresolved", evidence,
                "The displayed structure was inspected, but stereochemical identity and a complete machine representation remain unresolved.",
                raw="source drawing inspected; no stereochemical assignment claimed")
    molecules[label] = m
    bundle["molecules"].append(m)
    return m["molecule_id"]


def q(number, unit, raw_text=None, qualifier=None):
    d = {"number": number, "unit": unit}
    if raw_text is not None:
        d["raw_text"] = raw_text
    if qualifier is not None:
        d["qualifier"] = qualifier
    return d


def participant(r, key, label, role, evidence, *, stage="reaction", name=None,
                equiv=None, molpct=None, mass=None, amount=None, volume=None,
                reference=None, product_role=None, note=None, source_label=None):
    molecule_id = ensure_molecule(label, evidence, name=name, note=note, source_label=source_label)
    p = st.new_participant(f"p-{key}", molecule_id, role, evidence, stage=stage)
    if product_role is not None:
        p["product_role"] = f(product_role, product_role, evidence)
    values = {"equivalents": equiv, "mol_percent": molpct, "mass": mass,
              "amount": amount, "volume": volume}
    for field, value in values.items():
        if value is not None:
            raw = value.get("raw_text") if isinstance(value, dict) else value
            p["charges"][field] = f(raw, value, evidence)
    if reference is not None:
        ref_mid = ensure_molecule(reference, evidence)
        p["equivalence_reference"] = f(reference, ref_mid, evidence)
    r["participants"].append(p)
    return p


def classification(r, evidence, mode="cycloaddition"):
    c = r["classification"]
    c["model_votes"] = []
    c["mapped_reaction_smiles"] = f()
    if mode == "none":
        c["bond_changes"] = []
        c["labels"] = inferred([], evidence, "No observed product-forming transformation; no formation label is assigned.", raw="0% conversion/yield or no reaction")
        c["coverage"] = inferred("not_applicable", evidence, "No observed transformation is classified.", raw="no reaction")
        c["limitations"] = ["Negative/control outcome; intended-product bonds were not inferred."]
        return
    if mode == "unresolved":
        c["bond_changes"] = []
        c["labels"] = f()
        c["coverage"] = inferred("unresolved", evidence, "Exact atom mapping and net bond-order/hydrogen changes could not be assigned safely from the displayed rearrangement.", raw="source arrow inspected")
        c["limitations"] = ["No atom-mapped structures were supplied or invented.", "Hydrogen migration and double-bond relocation require expert mapping."]
        return
    changes = [
        {"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "N"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
        {"change_id": "b2", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
    ]
    deriv = "The [8+2] cycloaddition drawings show one new C-N and one new C-C sigma connection. Atom identities, bond-order relocations, and the reported 1,5-H shift remain unmapped."
    c["bond_changes"] = [derived(x, evidence, deriv) for x in changes]
    c["labels"] = derived([
        {"category_number": 5, "code": "C-N", "supporting_change_ids": ["b1"]},
        {"category_number": 2, "code": "C-C", "supporting_change_ids": ["b2"]},
    ], evidence, deriv)
    c["coverage"] = inferred("partial", evidence, "Only the two visually explicit ring-closing connections are classified; hydrogen and bond-order changes are unresolved.", raw="[8+2] cycloaddition arrow and product drawing")
    c["limitations"] = ["No complete atom mapping was asserted.", "The 1,5-hydrogen shift is discussed by the authors but its exact H source/target is not classified here.", "Other bond-order changes remain unresolved."]


def base_reaction(suffix, evidence, categories, location, name, transformation,
                  record_type="specific_experiment"):
    ev = list(dict.fromkeys(evidence))
    r = st.new_reaction(f"{DOC}:{suffix}", DOC, ev)
    r["record_type"] = record_type
    r["experiment_categories"] = f(categories, categories, ev)
    r["location"] = f(location, location, ev)
    r["reaction_name"] = f(name, name, ev)
    r["transformation_type"] = inferred(transformation, ev, "Normalized from the source reaction name and displayed reactant/product connectivity.", raw=name)
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                   "checked_fields": [], "unresolved": []}
    r["completeness"] = "partial"
    return r


def set_solvent_and_step(r, solvent, evidence, *, concentration=None, volume=None,
                         temperature=None, time=None, operation="combine and hold"):
    sp = participant(r, "solvent", solvent, "solvent", evidence, volume=volume)
    r["solvent_conditions"]["composition"] = f(solvent, [sp["molecule_id"]], evidence)
    if volume is not None:
        r["solvent_conditions"]["total_volume"] = f(volume.get("raw_text"), volume, evidence)
    if concentration is not None:
        r["solvent_conditions"]["concentration"] = f(concentration.get("raw_text"), concentration, evidence)
        r["solvent_conditions"]["concentration_basis"] = f("relative to the 1 equiv unsaturated partner", "1 equiv unsaturated partner", evidence)
    s = st.new_step("reaction", evidence)
    s["participant_ids"] = [p["participant_id"] for p in r["participants"]]
    s["operation"] = f(operation, operation, evidence)
    if temperature is not None:
        if isinstance(temperature, dict) and temperature.get("conflict"):
            alts = [f(a["raw"], a["value"], a["evidence"]) for a in temperature["alternatives"]]
            s["temperature"] = f(temperature["raw"], None, evidence, status="conflict", alternatives=alts)
        else:
            raw = temperature.get("raw_text") if isinstance(temperature, dict) else temperature
            s["temperature"] = f(raw, temperature, evidence)
    if time is not None:
        raw = time.get("raw_text") if isinstance(time, dict) else time
        s["time"] = f(raw, time, evidence)
    r["steps"].append(s)
    return s


def add_yield(r, target_labels, percent, evidence, basis=None, outcome_id="yield", scope=None):
    targets = [ensure_molecule(x, evidence) for x in target_labels]
    o = st.new_outcome(outcome_id, "yield", targets, evidence)
    o["product_role"] = f("product" if len(targets) == 1 else "combined product pair", "product" if len(targets) == 1 else "combined product pair", evidence)
    o["value"] = f(f"{percent}%", q(percent, "%", f"{percent}%"), evidence)
    if basis:
        o["yield_type"] = f(basis, basis, evidence)
        o["measurement_basis"] = f(basis, basis, evidence)
    if scope:
        o["scope"] = f(scope, scope, evidence)
    r["outcomes"].append(o)
    return o


def add_conversion(r, target_label, percent, evidence, basis=None):
    target = ensure_molecule(target_label, evidence)
    o = st.new_outcome("conversion", "conversion", [target], evidence)
    o["product_role"] = f("starting material", "starting material", evidence)
    o["value"] = f(f"{percent}%", q(percent, "%", f"{percent}%"), evidence)
    if basis:
        o["measurement_basis"] = f(basis, basis, evidence)
    r["outcomes"].append(o)


def add_selectivity(r, labels, ratio, evidence, method, order):
    targets = [ensure_molecule(x, evidence) for x in labels]
    o = st.new_outcome("selectivity", "selectivity", targets, evidence)
    o["product_role"] = f("isomer ratio", "isomer ratio", evidence)
    o["value"] = f(ratio, ratio, evidence)
    o["measurement_basis"] = f(method, method, evidence) if method else f()
    o["scope"] = f({"metric": "isomer ratio", "order": order}, {"metric": "isomer ratio", "order": order}, evidence)
    r["outcomes"].append(o)


def add_no_reaction(r, evidence, decomposition=False, target="1a"):
    o = st.new_outcome("no-reaction", "no_reaction", [], evidence)
    o["value"] = f("0% conversion/yield; no reaction", True, evidence)
    o["scope"] = f("observed transformation", "observed transformation", evidence)
    r["outcomes"].append(o)
    if decomposition:
        d = st.new_outcome("decomposition", "decomposition", [ensure_molecule(target, evidence)], evidence)
        d["value"] = f("8-azaheptafulvene decomposition", True, evidence)
        d["scope"] = f("starting azaheptafulvene", "starting azaheptafulvene", evidence)
        r["outcomes"].append(d)


def finish(r, class_mode="cycloaddition", unresolved=()):
    classification(r, r["evidence_ids"], class_mode)
    r["review"]["unresolved"].extend(unresolved)
    bundle["reactions"].append(r)
    return r


# Source-wide identity notes.
ensure_molecule("1a", [E[2], E[3], E[4], E[5]],
    name="8-aryl-8-azaheptafulvene 1a",
    note="Identity wording conflict retained: prose/Typical Procedure calls 1a '8-phenyl-8-azaheptafulvene', while the displayed structure and formula-implied mass show a para-methyl aryl substituent. Candidate SMILES follows the drawing, not a silent textual correction.")
ensure_molecule("2a", [E[2], E[3], E[5]], name="tosylallenamide 2a")
ensure_molecule("3a", [E[4], E[5]], name="tosylynamide 3a")


# Table 1: seventeen independent optimization experiments.
table1 = [
    (1, "(ArO)3PAuNTf2", "DCM", "rt", q(18, "h", "18 h"), 100, 40, "NMR"),
    (2, "(ArO)3PAuNTf2", "DCE", q(84, "°C", "84 °C"), q(15, "min", "15 min"), 100, 72, "NMR"),
    (3, "(ArO)3PAuNTf2", "DCE", q(60, "°C", "60 °C"), q(1, "h", "1 h"), 100, 82, "isolated"),
    (4, "(ArO)3PAuNTf2", "toluene", q(60, "°C", "60 °C"), q(1, "h", "1 h"), 60, 49, "NMR"),
    (5, "(ArO)3PAuNTf2", "hexane", q(60, "°C", "60 °C"), q(1, "h", "1 h"), 80, 36, "NMR"),
    (6, "(ArO)3PAuNTf2", "THF", q(60, "°C", "60 °C"), q(1, "h", "1 h"), 60, 46, "NMR"),
    (7, "(ArO)3PAuNTf2", "acetonitrile", q(60, "°C", "60 °C"), q(1, "h", "1 h"), 64, 51, "NMR"),
    (8, "JohnPhosAuNTf2", "DCE", q(60, "°C", "60 °C"), q(1, "h", "1 h"), 11, 7, "NMR"),
    (9, "IPrAuNTf2", "DCE", q(60, "°C", "60 °C"), q(1, "h", "1 h"), 0, 0, None),
    (10, "PicAuCl2", "DCE", q(60, "°C", "60 °C"), q(1, "h", "1 h"), 37, 40, "NMR"),
    (11, None, "DCE", q(60, "°C", "60 °C"), q(1, "h", "1 h"), 0, 0, None),
    (12, "AgNTf2", "DCE", q(60, "°C", "60 °C"), q(1, "h", "1 h"), 0, 0, None),
    (13, "ZnCl2", "DCE", q(60, "°C", "60 °C"), q(1, "h", "1 h"), 0, 0, None),
    (14, "In(OTf)3", "DCE", q(60, "°C", "60 °C"), q(1, "h", "1 h"), 0, 0, None),
    (15, "HBF4·Et2O", "DCE", q(60, "°C", "60 °C"), q(1, "h", "1 h"), 0, 0, None),
    (16, "(ArO)3PAuNTf2", "DCE", q(60, "°C", "60 °C"), q(1.5, "h", "1.5 h"), 100, 95, "isolated"),
    (17, "(ArO)3PAuNTf2", "DCE", q(60, "°C", "60 °C"), q(4, "h", "4 h"), 100, 96, "isolated"),
]
for row, catalyst, solvent, temp, tm, conv, yld, basis in table1:
    ev = [E[2]]
    if row == 16:
        ev += [E[3], E[5]]
    r = base_reaction(f"table1-entry-{row:02d}", ev,
        ["optimization"] + (["substrate_scope"] if row == 16 else []) + (["scale_up"] if row == 17 else []),
        {"page": 2, "region": f"Table 1 entry {row}"},
        "Gold(I)-catalyzed [8+2] cycloaddition of allenamide 2a and azaheptafulvene 1a",
        "[8+2] cycloaddition")
    eq1 = 1.3 if row in (16, 17) else 1.1
    amount1 = q(0.26, "mmol", "0.26 mmol") if row == 16 else q(1, "mmol", "1 mmol") if row == 17 else q(0.2, "mmol", "0.2 mmol")
    p1 = participant(r, "1a", "1a", "reactant", ev, equiv=q(eq1, "equiv", f"{eq1} equiv"), amount=amount1, reference="2a",
                     mass=q(47.1, "mg", "47.1 mg") if row == 16 else None)
    p2 = participant(r, "2a", "2a", "reactant", ev, equiv=q(1, "equiv", "1 equiv"), reference="2a",
                     amount=q(0.2, "mmol", "0.2 mmol") if row == 16 else None,
                     mass=q(44.7, "mg", "44.7 mg") if row == 16 else None)
    if conv > 0 or yld > 0:
        participant(r, "4a", "4a", "product", ev, product_role="main product",
                    name="1,8-dihydrocyclohepta[b]pyrrole 4a",
                    note="Product drawing and characterization are retained, but no complete source-verified SMILES was assigned.")
    loading = 2.5 if row == 17 else 10 if row in (12, 13, 14, 15) else 5
    if catalyst:
        participant(r, "catalyst", catalyst, "catalyst", ev,
                    molpct=q(loading, "mol%", f"{loading} mol%"), reference="2a",
                    mass=q(11.0, "mg", "11.0 mg") if row == 16 else None)
    conc = q(0.05 if row == 16 else 0.1, "M", "0.05 M" if row == 16 else "0.1 M")
    vol = q(4, "mL", "4 mL") if row == 16 else None
    set_solvent_and_step(r, solvent, ev, concentration=conc, volume=vol, temperature=temp, time=tm)
    r["reaction_scale"] = f(
        "2a 0.2 mmol" if row == 16 else "Table 1 footnote [k]: 1a (1 mmol, 1.3 equiv)" if row == 17 else "1a 0.2 mmol; 2a 1 equiv",
        {"reference_molecule_id": p2["molecule_id"], "amount": q(0.2, "mmol", "0.2 mmol")} if row == 16 else {"reference_molecule_id": p2["molecule_id"], "amount": "unresolved scale-up basis"} if row == 17 else {"reference_molecule_id": p1["molecule_id"], "amount": q(0.2, "mmol", "0.2 mmol")}, ev)
    if conv:
        add_conversion(r, "2a", conv, ev, "NMR with 4-bromobenzaldehyde as internal standard")
    if yld:
        add_yield(r, ["4a"], yld, ev, basis)
    else:
        add_no_reaction(r, ev, decomposition=row in (12, 13, 14, 15))
    r["procedure"]["original_text"] = f({"table": 1, "entry": row, "catalyst": catalyst or "none", "solvent": solvent,
                                               "temperature": temp, "time": tm, "conversion_percent": conv, "yield_percent": yld},
                                              {"table": 1, "entry": row}, ev)
    unresolved = []
    if row == 10:
        unresolved.append("Table 1 reports 37% conversion and 40% NMR yield; both values are retained without correction.")
    if row == 17:
        unresolved.append("Scale-up wording is ambiguous: footnote [k] says 1a (1 mmol, 1.3 equiv) while the prose says the reaction was scaled to 1 mmol; denominator amount is not silently inferred.")
    if row == 16:
        participant(r, "pPh3", "triphenylphosphine", "additive", [E[5]], stage="workup", mass=q(5, "mg", "5 mg"), molpct=q(10, "mol%", "10 mol%"), reference="2a")
        participant(r, "alumina", "neutral alumina", "purification_agent", [E[5]], stage="purification")
        participant(r, "hexane-purif", "hexane", "purification_agent", [E[5]], stage="purification")
        participant(r, "etoac", "ethyl acetate", "purification_agent", [E[5]], stage="purification")
        work = st.new_step("quench-and-evaporate", [E[5]], stage="workup")
        work["participant_ids"] = ["p-pPh3"]
        work["operation"] = f("triphenylphosphine was added, then the solvent was removed", "add triphenylphosphine then remove solvent", [E[5]])
        pur = st.new_step("purification", [E[5]], stage="purification")
        pur["participant_ids"] = ["p-alumina", "p-hexane-purif", "p-etoac"]
        pur["operation"] = f("neutral alumina column chromatography (hexane/EtOAc 5:1)", "neutral alumina column chromatography; hexane/ethyl acetate 5:1", [E[5]])
        r["steps"] += [work, pur]
        r["procedure"]["original_text"] = f("Typical procedure for 4a: 1a (47.1 mg, 0.26 mmol, 1.3 equiv), 2a (44.7 mg, 0.2 mmol, 1 equiv), catalyst (11.0 mg, 5 mol%) in DCE (4 mL, 0.05 M), 60 °C, 1.5 h; PPh3 (5 mg, 10 mol%), evaporation, neutral alumina (hexane/EtOAc 5:1).", "Typical procedure for 4a", [E[5]])
        # Isolated product amount is a separate outcome.
        o = st.new_outcome("product-amount", "product_amount", [ensure_molecule("4a", ev)], [E[5]])
        o["product_role"] = f("main product", "main product", [E[5]])
        o["value"] = f("79.5 mg", q(79.5, "mg", "79.5 mg"), [E[5]])
        o["measurement_basis"] = f("isolated yellow solid", "isolated", [E[5]])
        r["outcomes"].append(o)
    finish(r, "cycloaddition" if (conv or yld) else "none", unresolved)


# Scheme 1 scope entries other than 4a (already represented by Table 1 entry 16).
scope1 = [
    ("4b", "1b", "2a", q(1.5, "h", "1.5 h"), q(60, "°C", "60 °C"), 80, None, None),
    ("4c", "1c", "2a", q(3, "h", "3 h"), q(60, "°C", "60 °C"), 79, None, None),
    ("4d", "1d", "2a", q(3, "h", "3 h"), q(60, "°C", "60 °C"), 84, None, None),
    ("4e", "1a", "2-for-4e", q(14, "h", "14 h"), q(60, "°C", "60 °C"), 80, None, None),
    ("4f", "1a", "2-for-4f", q(12, "h", "12 h"), q(60, "°C", "60 °C"), 73, None, None),
    ("4g", "1a", "2-for-4g", q(12, "h", "12 h"), q(60, "°C", "60 °C"), 85, None, None),
    ("4h", "1a", "2-for-4h", q(12, "h", "12 h"), q(60, "°C", "60 °C"), 80, None, None),
    ("4i", "1a", "2-for-4i", q(12, "h", "12 h"), q(60, "°C", "60 °C"), 67, "5i", "10:1"),
    ("4j", "1a", "2-for-4j", q(2, "h", "2 h"), q(60, "°C", "60 °C"), 71, "5j", "10:1"),
    ("4k", "1a", "2-for-4k", q(12, "h", "12 h"), q(60, "°C", "60 °C"), 72, "5k", "10:1"),
    ("4l", "1a", "2-for-4l", q(12, "h", "12 h"), q(60, "°C", "60 °C"), 70, "5l", "4:1"),
    ("4m", "1a", "2-for-4m", q(3.5, "h", "3.5 h"), q(40, "°C", "40 °C"), 38, None, None),
    ("4n", "1a", "2-for-4n", q(12, "h", "12 h"), q(60, "°C", "60 °C"), 71, None, None),
    ("4o", "1a", "2-for-4o", q(12, "h", "12 h"), q(60, "°C", "60 °C"), 72, None, None),
    ("4p", "1a", "2p", q(6, "h", "6 h"), q(40, "°C", "40 °C"), 41, None, None),
]
for product, aza, allene, tm, temp, yld, minor, ratio in scope1:
    ev = [E[3], E[2]]
    r = base_reaction(f"scheme1-{product}", ev, ["substrate_scope"],
                      {"page": 3, "region": f"Scheme 1 product {product}"},
                      "Cycloaddition reaction of 8-aryl-8-azaheptafulvenes 1 and allenamides 2", "[8+2] cycloaddition")
    participant(r, "aza", aza, "reactant", ev, equiv=q(1.3, "equiv", "1.3 equiv"), reference=allene,
                source_label=f"1 precursor for {product}")
    participant(r, "allene", allene, "reactant", ev, equiv=q(1, "equiv", "1 equiv"), reference=allene,
                source_label=f"2 precursor for {product}")
    participant(r, "major", product, "product", [E[3]], product_role="major product")
    targets = [product]
    if minor:
        participant(r, "minor", minor, "product", [E[3]], product_role="minor 1,4-dihydro isomer")
        targets.append(minor)
    participant(r, "catalyst", "(ArO)3PAuNTf2", "catalyst", ev, molpct=q(5, "mol%", "5 mol%"), reference=allene)
    set_solvent_and_step(r, "DCE", ev, concentration=q(0.05, "M", "0.05 M"), temperature=temp, time=tm)
    add_yield(r, targets, yld, [E[3]], basis=None,
              scope="combined yield of 4/5 pair" if minor else "yield printed beside product in Scheme 1")
    if ratio:
        add_selectivity(r, [product, minor], ratio, [E[3]], None, f"{product}:{minor}")
    finish(r, unresolved=["Scheme 1 does not state the measurement basis for this yield." if product != "4p" else "Scheme 1 gives the yield but not its measurement basis; the substrate structure is retained from the drawing without a full SMILES."])


# Scheme 2 mechanistic time-course experiments.
for suffix, aza, allene, temp, tm, products, ratio in [
    ("1c-2a-0p5h", "1c", "2a", q(60, "°C", "60 °C"), q(0.5, "h", "0.5 h"), ["6c", "4c"], "1.2:1"),
    ("1c-2a-1p5h", "1c", "2a", q(60, "°C", "60 °C"), q(1.5, "h", "1.5 h"), ["6c", "4c"], "1:6"),
    ("1c-2a-3h", "1c", "2a", q(60, "°C", "60 °C"), q(3, "h", "3 h"), ["4c"], "6c not observed; 4c only"),
    ("1a-2p-rt-27h", "1a", "2p", "rt", q(27, "h", "27 h"), ["6p", "4p"], "2.6:1"),
    ("1a-2p-40c-3h", "1a", "2p", q(40, "°C", "40 °C"), q(3, "h", "3 h"), ["6p", "4p"], "9.2:1"),
    ("1a-2p-40c-6h", "1a", "2p", q(40, "°C", "40 °C"), q(6, "h", "6 h"), ["4p"], "6p not observed; 4p only"),
]:
    ev = [E[3]] + ([E[6]] if suffix == "1a-2p-40c-3h" else [])
    r = base_reaction(f"scheme2-{suffix}", ev, ["mechanistic_study"],
                      {"page": 3, "region": f"Scheme 2 row {suffix}"},
                      "Formation and transformation of intermediate cycloadducts 6", "[8+2] cycloaddition/time-course isomerization")
    participant(r, "aza", aza, "reactant", ev, source_label=aza)
    participant(r, "unsaturated", allene, "reactant", ev, source_label=allene)
    for i, prod in enumerate(products):
        participant(r, f"product-{i+1}", prod, "product", ev, product_role="observed isomer")
    participant(r, "catalyst", "(ArO)3PAuNTf2", "catalyst", ev, molpct=q(5, "mol%", "5 mol%"), reference=allene)
    set_solvent_and_step(r, "DCE", ev, concentration=q(0.05, "M", "0.05 M"), temperature=temp, time=tm)
    if len(products) == 2:
        add_selectivity(r, products, ratio, [E[3]], "1H NMR of crude reaction mixture", f"{products[0]}:{products[1]}")
    else:
        add_selectivity(r, products, ratio, [E[3]], "1H NMR of crude reaction mixture", products[0])
    if suffix == "1a-2p-40c-3h":
        add_yield(r, ["6p"], 52, [E[6]], "isolated", scope="isolated intermediate 6p")
    finish(r, unresolved=["Time-course ratios are crude-mixture 1H NMR results; absolute stereochemistry is not assigned."])

# Independent thermal conversion of isolated 6p to 4p without gold.
r = base_reaction("scheme2-isolated-6p-to-4p", [E[3]], ["mechanistic_study"],
                  {"page": 3, "region": "text following Scheme 2"},
                  "Transformation of isolated 6p to 4p without gold catalyst", "1,5-hydrogen shift/isomerization")
participant(r, "6p", "6p", "reactant", [E[3]])
participant(r, "4p", "4p", "product", [E[3]], product_role="isomerization product")
set_solvent_and_step(r, "DCE", [E[3]], temperature=q(40, "°C", "40 °C"), time=q(15, "h", "15 h"), operation="heat isolated 6p without gold catalyst")
o = st.new_outcome("conversion", "conversion", [ensure_molecule("6p", [E[3]])], [E[3]])
o["value"] = f("quantitative transformation", "quantitative", [E[3]])
o["product_role"] = f("starting intermediate", "starting intermediate", [E[3]])
r["outcomes"].append(o)
finish(r, "unresolved", ["Quantitative is qualitative in the source; it is not converted to an invented numeric percentage.", "Exact hydrogen/bond-order mapping is unresolved."])


# Table 2 optimization experiments; entry 3 is also Scheme 4 product 7a and the page-5 Typical Procedure.
table2 = [
    (1, "(ArO)3PAuNTf2", "DCE", q(60, "°C", "60 °C"), q(21, "h", "21 h"), 79, 90, "isolated"),
    (2, "JohnPhosAuNTf2", "DCE", q(60, "°C", "60 °C"), q(21, "h", "21 h"), 52, 56, "NMR with 4-bromobenzaldehyde as internal standard"),
    (3, "JohnPhosAuNTf2", "toluene", q(110, "°C", "110 °C"), q(4.5, "h", "4.5 h"), 90, 100, "isolated"),
    (4, None, "toluene", q(110, "°C", "110 °C"), q(4.5, "h", "4.5 h"), 0, None, None),
    (5, "HBF4·Et2O", "toluene", q(110, "°C", "110 °C"), q(4.5, "h", "4.5 h"), 0, None, None),
    (6, "In(OTf)3", "toluene", q(110, "°C", "110 °C"), q(4.5, "h", "4.5 h"), 0, None, None),
]
for row, catalyst, solvent, temp, tm, yld, conv, basis in table2:
    ev = [E[4]] + ([E[5]] if row == 3 else [])
    r = base_reaction(f"table2-entry-{row:02d}", ev,
        ["optimization"] + (["substrate_scope"] if row == 3 else []),
        {"page": 4, "region": f"Table 2 entry {row}"},
        "Cycloaddition reaction between 8-aryl-8-azaheptafulvene 1a and ynamide 3a", "[8+2] cycloaddition")
    p1 = participant(r, "1a", "1a", "reactant", ev, amount=q(0.3, "mmol", "0.3 mmol") if row == 3 else q(0.2, "mmol", "0.2 mmol"), equiv=q(1.5, "equiv", "1.5 equiv"), reference="3a", mass=q(54.3, "mg", "54.3 mg") if row == 3 else None)
    p3 = participant(r, "3a", "3a", "reactant", ev, amount=q(0.2, "mmol", "0.2 mmol") if row == 3 else None, equiv=q(1, "equiv", "1 equiv"), reference="3a", mass=q(57.0, "mg", "57.0 mg") if row == 3 else None)
    if yld:
        participant(r, "7a", "7a", "product", ev, product_role="main product", name="dihydrocyclohepta[b]pyrrole 7a",
                    note="Product-name locants are inconsistent across the supplied pages (abstract/conclusion describe a 1,4-dihydro product, while the Typical Procedure heading says 2-amido-1,6-dihydro); the drawing/label 7a is retained without choosing a corrected locant name.")
    if catalyst:
        participant(r, "catalyst", catalyst, "catalyst", ev, molpct=q(5, "mol%", "5 mol%"), reference="3a", mass=q(8.0, "mg", "8.0 mg") if row == 3 else None)
    actual_temp = temp
    if row == 3:
        actual_temp = {"conflict": True, "raw": ["110 °C (Table 2/Scheme 4)", "100 °C (Typical Procedure)"],
                       "alternatives": [
                           {"raw": "110 °C", "value": q(110, "°C", "110 °C"), "evidence": [E[4]]},
                           {"raw": "100 °C", "value": q(100, "°C", "100 °C"), "evidence": [E[5]]},
                       ]}
    set_solvent_and_step(r, solvent, ev, concentration=q(0.1, "M", "0.1 M"), volume=q(2, "mL", "2 mL") if row == 3 else None, temperature=actual_temp, time=tm)
    r["reaction_scale"] = f("3a (0.2 mmol, 1 equiv)", {"reference_molecule_id": p3["molecule_id"], "amount": q(0.2, "mmol", "0.2 mmol")}, ev)
    if yld:
        add_yield(r, ["7a"], yld, ev, basis)
        if conv is not None:
            add_conversion(r, "3a", conv, ev, "conversion stated in Table 2 footnote")
    else:
        add_no_reaction(r, ev)
    r["procedure"]["original_text"] = f({"table": 2, "entry": row, "catalyst": catalyst or "none", "solvent": solvent,
                                               "temperature": temp, "time": tm, "yield_percent": yld},
                                              {"table": 2, "entry": row}, ev)
    unresolved = []
    if row == 3:
        unresolved.append("Source temperature conflict retained: Table 2 and Scheme 4 state 110 °C; the Typical Procedure states 100 °C.")
        unresolved.append("7a nomenclature conflict retained: main-text overview uses 1,4-dihydro wording, while the Typical Procedure heading uses 2-amido-1,6-dihydro wording.")
        participant(r, "pPh3", "triphenylphosphine", "additive", [E[5]], stage="workup", mass=q(5, "mg", "5 mg"), molpct=q(10, "mol%", "10 mol%"), reference="3a")
        participant(r, "silica", "silica gel", "purification_agent", [E[5]], stage="purification")
        participant(r, "hexane-purif", "hexane", "purification_agent", [E[5]], stage="purification")
        participant(r, "etoac", "ethyl acetate", "purification_agent", [E[5]], stage="purification")
        work = st.new_step("quench-and-evaporate", [E[5]], stage="workup")
        work["participant_ids"] = ["p-pPh3"]
        work["operation"] = f("triphenylphosphine was added and the solvent was removed", "add triphenylphosphine then remove solvent", [E[5]])
        pur = st.new_step("purification", [E[5]], stage="purification")
        pur["participant_ids"] = ["p-silica", "p-hexane-purif", "p-etoac"]
        pur["operation"] = f("silica column chromatography (hexane/EtOAc 20:1)", "silica column chromatography; hexane/ethyl acetate 20:1", [E[5]])
        r["steps"] += [work, pur]
        r["procedure"]["original_text"] = f("Typical procedure for 7a: 1a (54.3 mg, 0.3 mmol, 1.5 equiv), 3a (57.0 mg, 0.2 mmol, 1 equiv), JohnPhosAuNTf2 (8.0 mg, 5 mol%) in toluene (2 mL, 0.1 M), 100 °C, 4.5 h; PPh3 (5 mg, 10%), evaporation, silica (hexane/EtOAc 20:1).", "Typical procedure for 7a", [E[5]])
        o = st.new_outcome("product-amount", "product_amount", [ensure_molecule("7a", ev)], [E[5]])
        o["product_role"] = f("main product", "main product", [E[5]])
        o["value"] = f("86.4 mg", q(86.4, "mg", "86.4 mg"), [E[5]])
        o["measurement_basis"] = f("isolated yellow solid", "isolated", [E[5]])
        r["outcomes"].append(o)
    finish(r, "cycloaddition" if yld else "none", unresolved)


# Scheme 4 scope, excluding 7a which is Table 2 entry 3.
scope4 = [
    ("7b", "1b", "3a", q(6, "h", "6 h"), 62, False, "tosyl ynamide"),
    ("7c", "1c", "3a", q(5.5, "h", "5.5 h"), 79, False, "tosyl ynamide"),
    ("7d", "1a", "3-for-7d", None, 0, False, "R4=methyl; complex mixture described in text"),
    ("7e", "1a", "3-for-7e", None, 0, False, "R4=methylethenyl; complex mixture described in text"),
    ("7f", "1a", "3f", q(45, "min", "45 min"), 22, False, "terminal ynamide, R4=H"),
    ("7f-route-3m", "1a", "3m", q(40, "min", "40 min"), 59, True, "R4=TMS; unstable adduct followed by desilylation"),
    ("7g", "1a", "3-for-7g", q(40, "min", "40 min"), 76, False, "carboxylate-substituted ynamide"),
    ("7h", "1a", "3-for-7h", q(1, "h", "1 h"), 76, False, "oxazolidinone ynamide, R4=Ph"),
    ("7i", "1a", "3-for-7i", q(2.5, "h", "2.5 h"), 74, False, "oxazolidinone ynamide, R4=p-tolyl"),
    ("7j", "1a", "3-for-7j", q(1, "h", "1 h"), 86, False, "oxazolidinone ynamide, R4=1-naphthyl"),
    ("7k", "1a", "3-for-7k", q(1.5, "h", "1.5 h"), 74, False, "oxazolidinone ynamide, R4=p-bromophenyl"),
    ("7l", "1a", "3-for-7l", q(2.5, "h", "2.5 h"), 81, False, "oxazolidinone ynamide, R4=p-methoxyphenyl"),
]
for product_tag, aza, ynamide, tm, yld, desilylation, desc in scope4:
    product = "7f" if product_tag == "7f-route-3m" else product_tag
    ev = [E[4]]
    r = base_reaction(f"scheme4-{product_tag}", ev, ["substrate_scope"],
                      {"page": 4, "region": f"Scheme 4 product/route {product_tag}"},
                      "Cycloadducts from 8-aryl-8-azaheptafulvenes 1 and ynamides 3",
                      "[8+2] cycloaddition followed by desilylation" if desilylation else "[8+2] cycloaddition")
    participant(r, "aza", aza, "reactant", ev, equiv=q(1.5, "equiv", "1.5 equiv"), reference=ynamide,
                source_label=f"1 precursor for {product_tag}")
    participant(r, "ynamide", ynamide, "reactant", ev, equiv=q(1, "equiv", "1 equiv"), reference=ynamide,
                source_label=f"3 precursor for {product_tag}", note=desc)
    if yld:
        participant(r, "product", product, "product", ev, product_role="main product")
    participant(r, "catalyst", "JohnPhosAuNTf2", "catalyst", ev, molpct=q(5, "mol%", "5 mol%"), reference=ynamide)
    set_solvent_and_step(r, "toluene", ev, concentration=q(0.1, "M", "0.1 M"), temperature=q(110, "°C", "110 °C"), time=tm)
    if yld:
        add_yield(r, [product], yld, ev, None, scope=desc)
        unresolved = ["Scheme 4 gives the yield but does not state its measurement basis."]
        if desilylation:
            unresolved.append("The 59% route is described as desilylation of an unstable 3-TMS adduct; desilylation reagents/operations are not reported in the supplied main text.")
            class_mode = "cycloaddition"
        else:
            class_mode = "cycloaddition"
    else:
        add_no_reaction(r, ev)
        # Text specifies complex mixtures rather than a clean product.
        mix = st.new_outcome("complex-mixture", "other", [ensure_molecule(ynamide, ev)], ev)
        mix["value"] = f("complex mixture", "complex mixture", ev)
        mix["scope"] = f(desc, desc, ev)
        r["outcomes"].append(mix)
        unresolved = ["Reaction time and product identities within the complex mixture are not reported."]
        class_mode = "none"
    finish(r, class_mode, unresolved)


# Generic schemes and proposed mechanisms are separate from specific experiments.
for suffix, page, reactant_labels, product_labels, title in [
    ("generic-allene-cycloaddition", 3, ["class-1", "class-2"], ["class-4", "class-5"], "Cycloaddition reaction of 8-aryl-8-azaheptafulvenes 1 and allenamides 2"),
    ("generic-ynamide-cycloaddition", 4, ["class-1", "class-3"], ["class-7"], "Cycloadducts from 8-aryl-8-azaheptafulvenes 1 and ynamides 3"),
]:
    generic_ev = [E[1], E[page]]
    r = base_reaction(suffix, generic_ev, [], {"pages": [1, page], "region": "Figure 1c overview plus generic scope scheme header and arrow"}, title, "generic [8+2] cycloaddition", "generic_scheme")
    for i, label in enumerate(reactant_labels):
        participant(r, f"reactant-{i+1}", label, "reactant", generic_ev, source_label=label.replace("class-", ""), note="Generic R-group class; retained unresolved rather than expanded into a false concrete molecule.")
    for i, label in enumerate(product_labels):
        participant(r, f"product-{i+1}", label, "product", generic_ev, source_label=label.replace("class-", ""), product_role="generic product class")
    finish(r, "cycloaddition", ["Generic R-group scheme; not counted as a concrete experiment and no yield is assigned."])

for suffix, page, labels, title in [
    ("proposed-mechanism-allene", 3, ["class-1", "class-2", "intermediate-I", "intermediate-III", "class-6", "class-4"], "Proposed mechanism for azaheptafulvene/allenamide cycloaddition"),
    ("proposed-mechanism-ynamide", 5, ["class-1", "class-3", "intermediate-II", "intermediate-IV", "intermediate-V", "class-7"], "Proposed mechanism for azaheptafulvene/ynamide cycloaddition"),
]:
    r = base_reaction(suffix, [E[page]], [], {"page": page, "region": "proposed mechanism scheme"}, title, "author-proposed stepwise mechanism", "proposed_mechanism")
    for i, label in enumerate(labels):
        role = "reactant" if i < 2 else "product" if i == len(labels) - 1 else "unknown"
        participant(r, f"species-{i+1}", label, role, [E[page]], source_label=label, note="Mechanistic/generic species; structure retained visually but no concrete registry SMILES asserted.")
    finish(r, "cycloaddition", ["Author-proposed mechanism, not an independently verified elementary-step experiment."])


# Characterization is source context, not additional synthesis records.
bundle["source_context"] = [
    {"context_id": f"{DOC}:context:4a-characterization", "document_id": DOC,
     "kind": "characterization", "target_molecule_ids": [ensure_molecule("4a", [E[5]])],
     "evidence_ids": [E[5]], "summary": "4a: yellow solid, Rf, mp, 1H/13C NMR and HRMS (EI); retained as characterization of Table 1 entry 16, not a separate reaction."},
    {"context_id": f"{DOC}:context:7a-characterization", "document_id": DOC,
     "kind": "characterization", "target_molecule_ids": [ensure_molecule("7a", [E[5]])],
     "evidence_ids": [E[5]], "summary": "7a: yellow solid, Rf, mp, 1H/13C NMR and HRMS (ESI); retained as characterization of Table 2 entry 3, not a separate reaction."},
    {"context_id": f"{DOC}:context:page6-references", "document_id": DOC,
     "kind": "references_and_footnotes", "evidence_ids": [E[6]],
     "summary": "Page 6 is predominantly references. Footnotes 20-21 add instability/isolation context for 6c/6p; no standalone spectral or literature-reference reaction was imported as this paper's experiment."},
]


# Raw visual inventory retained separately from formal records.
inventory = {
    "document_id": DOC,
    "pages_visually_inspected": [1, 2, 3, 4, 5, 6],
    "table1_rows": 17,
    "table2_rows": 6,
    "scheme1_product_rows": ["4a", "4b", "4c", "4d", "4e", "4f", "4g", "4h", "4i/5i", "4j/5j", "4k/5k", "4l/5l", "4m", "4n", "4o", "4p"],
    "scheme2_timecourse_rows": 6,
    "scheme2_additional_text_experiments": 1,
    "scheme4_routes": ["7a", "7b", "7c", "7d", "7e", "7f from 3f", "7f from 3m", "7g", "7h", "7i", "7j", "7k", "7l"],
    "generic_schemes": 2,
    "proposed_mechanisms": 2,
    "temperature_conflict": {"reaction_id": f"{DOC}:table2-entry-03", "values": ["110 °C", "100 °C"]},
    "identity_conflict": {"molecule_id": mid("1a"), "issue": "text says 8-phenyl; drawing/mass indicate para-methyl aryl"},
}
(RUN / "raw_model/visual_inventory.json").write_text(json.dumps(inventory, ensure_ascii=False, indent=2) + "\n")
(RUN / "submitted.json").write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")
print(json.dumps({"molecules": len(bundle["molecules"]), "reactions": len(bundle["reactions"]),
                  "specific_experiments": sum(r["record_type"] == "specific_experiment" for r in bundle["reactions"])}, ensure_ascii=False))
