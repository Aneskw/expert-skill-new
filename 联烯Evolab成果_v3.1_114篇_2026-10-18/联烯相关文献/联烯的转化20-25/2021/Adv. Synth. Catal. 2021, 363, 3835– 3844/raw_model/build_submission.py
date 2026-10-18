#!/usr/bin/env python3
"""Build the source-specific v3.1 submission, including this session's pages 7-10."""
from __future__ import annotations

import json
import sys
from collections import OrderedDict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

RUN = ROOT / "runs/current"
DOC = "doc_8cf88f502ad0d68c"
E = {i: f"{DOC}:page{i}:text" for i in range(1, 11)}
I = {i: f"{DOC}:page{i}:image" for i in range(1, 11)}


def dump(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n")


def f(raw=None, normalized=None, evidence=(), method="original_extraction", derivation=None,
      status=None, verification=None, alternatives=()):
    return st.fact(raw, normalized, evidence, method=method, derivation=derivation,
                   status=status, verification=verification, alternatives=alternatives)


def inferred(normalized, evidence, derivation, raw=None):
    return f(raw, normalized, evidence, method="inference", derivation=derivation)


def derived(normalized, evidence, derivation, raw=None):
    return f(raw, normalized, evidence, method="structure_derivation", derivation=derivation)


def q(number, unit, **extra):
    return {"number": number, "unit": unit, **extra}


def raw_text(value):
    if isinstance(value, dict) and "number" in value and "unit" in value:
        return value.get("raw_text") or f"{value['number']:g} {value['unit']}"
    return str(value)


def make_molecule(mid, label, name, ev, smiles=None, allene=False, assignment=None,
                  chirality=None, geometry=None, representation=None, mixture=None):
    m = st.new_molecule(mid, DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f(label, label, ev)
    m["names"]["original"] = f(name, name, ev)
    m["names"]["english"] = f(name, name, ev)
    m["names"]["abbreviations"] = f()
    m["names"]["aliases"] = f()
    m["cas_numbers"] = f()
    m["identity_note"] = f(f"Source identity retained as {label}; no CAS number was reported.",
                             f"Source identity retained as {label}; no CAS number was reported.", ev)
    m["external_properties"] = []
    if smiles:
        diag = st.stereo_inspect(smiles)
        m["representation_type"] = "smiles"
        m["smiles"] = derived(smiles, ev,
            "Candidate connectivity transcribed from the source drawing/name; RDKit parsing is a graph check, not source verification.",
            raw=f"structure drawing/name for {label}")
        m["generic_template"] = None
        m["stereochemistry"]["atom_centers"] = derived(diag["atom_centers"], ev,
            "RDKit graph diagnostic on the candidate SMILES; it cannot infer unreported configuration.")
        m["stereochemistry"]["double_bonds"] = derived(diag["double_bonds"], ev,
            "RDKit graph diagnostic on the candidate SMILES; cumulene axial stereochemistry requires separate review.")
        if allene:
            m["stereochemistry"]["axial"] = derived(
                [{"fragment": "each C=C=C arm", "assignment": "unassigned",
                  "structural_note": "terminal CH2 end; no source axial descriptor"}], ev,
                "The candidate graph has terminal CH2 allene ends; the source supplies no axial descriptor.")
            ch = chirality or "unspecified"
            geo = geometry or "not_applicable"
            rep = representation or "partial"
        else:
            m["stereochemistry"]["axial"] = f()
            ch = chirality or ("specified" if assignment else "achiral")
            geo = geometry or "not_applicable"
            rep = representation or "encoded"
        m["stereochemistry"]["mixture"] = f(mixture, mixture, ev) if mixture else f()
        m["stereochemistry"]["assignment"] = f(assignment, assignment, ev) if assignment else f()
        m["stereochemistry"]["chirality_status"] = inferred(ch, ev,
            "Status conservatively combines the source drawing/name with candidate-graph diagnostics.")
        m["stereochemistry"]["geometry_status"] = inferred(geo, ev,
            "Status conservatively combines the source drawing/name with candidate-graph diagnostics.")
        m["stereochemistry"]["representation_status"] = inferred(rep, ev,
            "Cumulene/axial and unreported stereochemistry are not declared fully encoded." if allene else
            "Candidate SMILES has no unencoded source-specified stereochemistry.")
        m["stereochemistry"]["isomeric_smiles"] = derived(diag["canonical_isomeric_smiles"], ev,
            "RDKit canonical isomeric SMILES from the candidate graph; no new stereochemical assignment is introduced.")
        m["stereochemistry"]["identity_key"] = derived(diag["stereo_graph_key"], ev,
            "RDKit stereo graph key; this is not a sample identity or source verification.")
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = f()
        m["generic_template"] = None
        m["stereochemistry"]["atom_centers"] = f()
        m["stereochemistry"]["double_bonds"] = f()
        m["stereochemistry"]["axial"] = f()
        m["stereochemistry"]["mixture"] = f(mixture, mixture, ev) if mixture else f()
        m["stereochemistry"]["assignment"] = f(assignment, assignment, ev) if assignment else f()
        m["stereochemistry"]["chirality_status"] = inferred(chirality or "unresolved", ev,
            "No complete candidate graph was asserted; absent source stereomarks are not treated as achiral or racemic.")
        m["stereochemistry"]["geometry_status"] = inferred(geometry or "unresolved", ev,
            "No complete candidate graph was asserted; unreported geometry remains unresolved.")
        m["stereochemistry"]["representation_status"] = inferred(representation or "unresolved", ev,
            "The compact or complex source drawing was not encoded as a complete SMILES.")
        m["stereochemistry"]["isomeric_smiles"] = f()
        m["stereochemistry"]["identity_key"] = f()
    return m


def make_template(mid, label, name, ev, fixed, variables):
    m = make_molecule(mid, label, name, ev, chirality="unresolved", geometry="unresolved",
                      representation="not_encoded")
    m["representation_type"] = "generic_template"
    m["generic_template"] = {
        "fixed_structure": f(fixed, fixed, ev),
        "variables": [{
            "label": v,
            "attachment_points": f("source-marked attachment", ["source-marked attachment"], ev),
            "definition": f(desc, desc, ev),
            "allowed_values": f(),
            "constraints": f(),
        } for v, desc in variables],
        "instance_mapping": [],
        "reason": f("R-group drawing in Scheme 1", "The source presents an R-group family rather than one concrete molecule.", ev),
    }
    return m


def participant(pid, mid, role, ev, stage="reaction", mass=None, amount=None, volume=None,
                equivalents=None, mol_percent=None, reference=None, product_role=None):
    p = st.new_participant(pid, mid, role, ev, stage)
    p["product_role"] = f(product_role, product_role, ev) if product_role else f()
    values = {"mass": mass, "amount": amount, "volume": volume,
              "equivalents": equivalents, "mol_percent": mol_percent}
    for key, item in values.items():
        p["charges"][key] = item if isinstance(item, dict) and "normalized_value" in item else (
            f(str(item), item, ev) if item is not None else f())
    p["equivalence_reference"] = f(reference, reference, ev) if reference else f()
    return p


def step(sid, ev, participant_ids=(), stage="reaction", operation=None, temperature=None,
         time=None, atmosphere=None, pressure=None, inherited=False):
    s = st.new_step(sid, ev, stage)
    s["participant_ids"] = list(participant_ids)
    method = "inherited_general_condition" if inherited else "original_extraction"
    s["operation"] = f(operation, operation, ev, method=method) if operation else f()
    s["temperature"] = f(raw_text(temperature), temperature, ev, method=method) if temperature is not None else f()
    s["time"] = f(raw_text(time), time, ev, method=method) if time is not None else f()
    s["atmosphere"] = f(atmosphere, atmosphere, ev, method=method) if atmosphere else f()
    s["pressure"] = f(pressure, pressure, ev, method=method) if pressure else f()
    return s


def outcome(oid, kind, targets, ev, value=None, role=None, yield_type=None, basis=None, scope=None,
            method="original_extraction"):
    o = st.new_outcome(oid, kind, targets, ev)
    o["product_role"] = f(role, role, ev, method=method) if role else f()
    o["value"] = f(raw_text(value), value, ev, method=method) if value is not None else f()
    o["yield_type"] = f(yield_type, yield_type, ev, method=method) if yield_type else f()
    o["measurement_basis"] = f(basis, basis, ev, method=method) if basis else f()
    o["scope"] = f(scope, scope, ev, method=method) if scope else f()
    return o


def classify(r, mode="cascade"):
    c = r["classification"]
    c["mapped_reaction_smiles"] = f()
    c["model_votes"] = []
    if mode == "not_applicable":
        c["bond_changes"] = []
        c["labels"] = inferred([], r["evidence_ids"], "No observed target transformation is classified.")
        c["coverage"] = inferred("not_applicable", r["evidence_ids"], "No observed target transformation is classified.")
        c["limitations"] = ["No desired-product formation was observed; expected bonds were not assigned."]
        return
    if mode == "unresolved":
        c["bond_changes"] = []
        c["labels"] = inferred([], r["evidence_ids"], "Observed products are structurally unresolved.")
        c["coverage"] = inferred("unresolved", r["evidence_ids"], "Observed products are structurally unresolved.")
        c["limitations"] = ["Product identity is unresolved; no bond-formation label is asserted."]
        return
    changes = []
    if mode in ("cascade", "prior_cycloaddition"):
        count = 3 if mode == "cascade" else 2
        for i in range(count):
            changes.append({"change_id": f"b{i+1}", "atoms": [{"atom_map": None, "element": "C"},
                            {"atom_map": None, "element": "C"}], "before_order": 0,
                            "after_order": 1, "change_type": "formed"})
    elif mode == "hydrogenation":
        changes.append({"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"},
                        {"atom_map": None, "element": "C"}], "before_order": 2,
                        "after_order": 1, "change_type": "order_decreased"})
        for i in (2, 3):
            changes.append({"change_id": f"b{i}", "atoms": [{"atom_map": None, "element": "C"},
                            {"atom_map": None, "element": "H"}], "before_order": 0,
                            "after_order": 1, "change_type": "hydrogen_gain", "hydrogen_count": 1})
    labels = OrderedDict()
    for ch in changes:
        if ch["change_type"] in ("formed", "order_increased", "hydrogen_gain"):
            lab = st.bond_label(ch["atoms"][0]["element"], ch["atoms"][1]["element"], ch["after_order"])
            labels.setdefault(lab["code"], {**lab, "supporting_change_ids": []})["supporting_change_ids"].append(ch["change_id"])
    c["bond_changes"] = [derived(ch, r["evidence_ids"],
        "Class-level net bond change inferred from the source reaction drawing; exact atom mapping is unresolved.") for ch in changes]
    c["labels"] = derived(list(labels.values()), r["evidence_ids"],
        "Taxonomy labels are calculated from the retained partial bond-change list.")
    c["coverage"] = inferred("partial", r["evidence_ids"],
        "Only source-supported class-level changes are retained; exact atoms, hydrogen balance, and coproducts are unresolved.")
    c["limitations"] = [
        "No complete atom-mapped reaction was asserted.",
        "Exact atom identities, all bond-order changes, proton transfers, coproducts, and fullerene cage mapping remain unresolved.",
    ]


def reaction(suffix, ev, record_type="specific_experiment", categories=(), location=None,
             name=None, transformation=None, schemes=()):
    rid = f"{DOC}:{suffix}"
    r = st.new_reaction(rid, DOC, ev)
    r["record_type"] = record_type
    r["experiment_categories"] = f(list(categories), list(categories), ev) if categories else f()
    r["location"] = f(location, location, ev) if location else f()
    r["general_scheme_ids"] = list(schemes)
    r["reaction_name"] = f(name, name, ev) if name else f()
    r["transformation_type"] = inferred(transformation, ev,
        "Normalized from the source reaction name/drawing; this is not an author-supplied taxonomy label.") if transformation else f()
    r["reaction_scale"] = f()
    r["solvent_conditions"] = {k: f() for k in ("composition", "mixture_ratio", "total_volume", "concentration", "concentration_basis")}
    r["special_conditions"] = {k: f() for k in ("light_source", "wavelength", "power", "current", "voltage", "electrodes", "flow_parameters", "reactor")}
    r["procedure"] = {k: f() for k in ("original_text", "general_procedure_reference", "modifications", "addition_order", "pretreatment")}
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                   "checked_fields": [], "unresolved": [
                       "Machine extraction requires expert source review; schema validation does not establish chemistry accuracy.",
                       "No explicitly associated Supporting Information was supplied."]}
    r["completeness"] = "partial"
    return r


# Rebuild from the clean v3.1 init template, never from an earlier submission.
template = json.loads((RUN / "raw_model/standard_init_pages_0007_0010.json").read_text())
M = {}
def add(m): M[m["molecule_id"]] = m

# Generic families and recurrent species.
add(make_template(f"{DOC}:mol:generic_1", "1a-l", "N-tethered 1,5-bisallene family", [E[2], I[2]],
    {"drawing": "RN(CH2-CH=C=CH2)2", "source": "Scheme 1"}, [("R", "nitrogen substituent/tether group")]))
add(make_template(f"{DOC}:mol:generic_3", "3a-l", "C60 fused dihydroazepine adduct family", [E[2], I[2]],
    {"drawing": "C60 fused to a 6,7-membered addend", "source": "Scheme 1"}, [("R", "nitrogen substituent/tether group")]))

bis_smiles = {
    "1a": "Cc1ccc(S(=O)(=O)N(CC=C=C)CC=C=C)cc1",
    "1b": "COc1ccc(S(=O)(=O)N(CC=C=C)CC=C=C)cc1",
    "1c": "Brc1ccc(S(=O)(=O)N(CC=C=C)CC=C=C)cc1",
    "1d": "Fc1ccc(S(=O)(=O)N(CC=C=C)CC=C=C)cc1",
    "1e": "FC(F)(F)c1ccc(S(=O)(=O)N(CC=C=C)CC=C=C)cc1",
    "1f": "FC(F)(F)c1ccccc1S(=O)(=O)N(CC=C=C)CC=C=C",
    "1g": "COC(=O)c1ccc(S(=O)(=O)N(CC=C=C)CC=C=C)cc1",
    "1h": "O=[N+]([O-])c1ccc(S(=O)(=O)N(CC=C=C)CC=C=C)cc1",
    "1j": "c1cc(sc1)S(=O)(=O)N(CC=C=C)CC=C=C",
    "1k": "C[Si](C)(C)CCS(=O)(=O)N(CC=C=C)CC=C=C",
    "1l": "CC(C)(C)S(=O)(=O)N(CC=C=C)CC=C=C",
    "1m": "C(C(=O)OCC)(C(=O)OCC)(CC=C=C)CC=C=C",
    "1n": "C(C(=O)OC)(c1ccccc1)(CC=C=C)CC=C=C",
    "1o": "CS(=O)(=O)C(c1ccccc1)(CC=C=C)CC=C=C",
    "1p": "O=S(=O)(c1ccccc1)C(c1ccccc1)(CC=C=C)CC=C=C",
    "1O": "O(CC=C=C)CC=C=C",
}
for lab, smi in bis_smiles.items():
    page = 2 if lab == "1a" else 4
    nm = "oxygen-tethered 1,5-bisallene" if lab == "1O" else f"bisallene {lab}"
    add(make_molecule(f"{DOC}:mol:{lab}", lab, nm, [E[page], I[page]], smi, allene=True))
add(make_molecule(f"{DOC}:mol:1i", "1i", "heteroaromatic sulfonamide-tethered bisallene 1i", [E[4], I[4]],
                  chirality="unresolved", geometry="unresolved"))

add(make_molecule(f"{DOC}:mol:C60", "C60", "fullerene C60", [E[2], E[4], E[6]],
                  chirality="unresolved", geometry="not_applicable"))
add(make_molecule(f"{DOC}:mol:C70", "C70", "fullerene C70", [E[4], E[5], E[8], I[8]],
                  chirality="unresolved", geometry="not_applicable"))

for lab in ["3a"] + [f"3{x}" for x in "bcdefghijklmnop"]:
    if lab == "3a":
        ev = [E[2], I[2], E[3], I[3], E[4], I[4], E[6], I[6], E[9], I[9]]
    elif lab in ("3b", "3c"):
        ev = [E[4], I[4], E[6], I[6]] + ([E[7], I[7]] if lab == "3c" else [])
    elif lab in [f"3{x}" for x in "defghij"]:
        ev = [E[4], I[4], E[7], I[7]]
    elif lab == "3k":
        ev = [E[4], I[4], E[7], I[7], E[8], I[8]]
    else:
        ev = [E[4], I[4], E[8], I[8]]
    add(make_molecule(f"{DOC}:mol:{lab}", lab, f"C60 fused cycloadduct {lab}", ev,
                      chirality="unspecified", geometry="unspecified", representation="not_encoded"))

for lab, name, assign in [
    ("alpha4a", "alpha-C70 cycloadduct alpha-4a", "alpha C70 regioisomer assignment"),
    ("beta4a", "beta-C70 cycloadduct beta-4a", "beta C70 site-isomer assignment"),
    ("alpha6a", "reduced alpha-C70 cycloadduct alpha-6a", "reduced alpha C70 regioisomer assignment"),
    ("beta6a", "reduced beta-C70 cycloadduct beta-6a", "reduced beta C70 site-isomer assignment"),
]:
    product_ev = [E[5], I[5], E[8], I[8], E[9], I[9]] if "4a" in lab else [E[5], I[5], E[9], I[9]]
    add(make_molecule(f"{DOC}:mol:{lab}", lab.replace("alpha", "α-").replace("beta", "β-"), name, product_ev,
                      assignment=assign, chirality="unresolved", geometry="unspecified", representation="not_encoded"))
add(make_molecule(f"{DOC}:mol:5a", "5a", "selectively reduced C60 cycloadduct 5a", [E[5], I[5], E[9], I[9]],
                  chirality="unspecified", geometry="unspecified", representation="not_encoded"))

# Reagents/solvents/ligands.
simple = [
    ("H2", "H2", "hydrogen", "[H][H]"), ("CH2Cl2", "CH2Cl2", "dichloromethane", "ClCCl"),
    ("oDCB", "o-DCB", "1,2-dichlorobenzene", "Clc1ccccc1Cl"),
    ("toluene", "toluene", "toluene", "Cc1ccccc1"), ("CS2", "CS2", "carbon disulfide", "S=C=S"),
    ("CHCl3", "CHCl3", "chloroform", "ClC(Cl)Cl"), ("O2", "O2", "oxygen", "O=O"),
    ("Et3SiH", "Et3SiH", "triethylsilane", "CC[SiH](CC)CC"),
]
for key, lab, name, smi in simple:
    ev = [E[2], E[6], E[8], E[9]] if key in ("H2", "CH2Cl2", "oDCB") else ([E[5], E[9]] if key == "Et3SiH" else [E[3], E[6], E[8], E[9]])
    add(make_molecule(f"{DOC}:mol:{key}", lab, name, ev, smi))
for key, lab, name, ev, assignment in [
    ("Rh", "[Rh(cod)2]BF4", "bis(cyclooctadiene)rhodium(I) tetrafluoroborate", [E[2], E[6], E[8]], None),
    ("DTBM", "(R)-DTBM-Segphos", "(R)-DTBM-Segphos", [E[2], E[6], E[8]], "(R) ligand designation"),
    ("BINAP", "(R)-BINAP", "(R)-BINAP", [E[2]], "(R) ligand designation"),
    ("H8BINAP", "(R)-H8-BINAP", "(R)-H8-BINAP", [E[2]], "(R) ligand designation"),
    ("TolBINAP", "(R)-Tol-BINAP", "(R)-Tol-BINAP", [E[2]], "(R) ligand designation"),
    ("BF3Et2O", "BF3·Et2O", "boron trifluoride diethyl etherate", [E[5], E[9]], None),
    ("silica", "silica gel", "silica gel", [E[6], E[8], E[9]], None),
]:
    add(make_molecule(f"{DOC}:mol:{key}", lab, name, ev, assignment=assignment,
                      chirality="specified" if assignment else "unresolved", geometry="not_applicable",
                      representation="not_encoded" if assignment else "unresolved"))

# Unresolved entities for cited prior work.
for key, label, name in [
    ("prior_itoh_sub", "Itoh bisallene", "Itoh and co-workers bisallene"),
    ("prior_itoh_prod", "Itoh C60 adduct", "Itoh and co-workers fullerene adduct"),
    ("prior_maeda_sub", "Maeda allenoate", "Maeda and co-workers allenoate"),
    ("prior_maeda_prod", "Maeda C60 adduct", "Maeda and co-workers fullerene adduct"),
    ("prior_wang_sub", "Wang substrate", "Wang and co-workers chlorinated sulfonamide substrate"),
    ("prior_wang_prod", "Wang C60 adduct", "Wang and co-workers fullerene adduct"),
]:
    add(make_molecule(f"{DOC}:mol:{key}", label, name, [E[2], I[2]]))

R = []
GEN = f"{DOC}:generic_scheme"

for suffix, sub, prod, source_name in [
    ("literature_itoh", "prior_itoh_sub", "prior_itoh_prod", "Itoh and co-workers"),
    ("literature_maeda", "prior_maeda_sub", "prior_maeda_prod", "Maeda and co-workers"),
    ("literature_wang", "prior_wang_sub", "prior_wang_prod", "Wang and co-workers"),
]:
    r = reaction(suffix, [E[2], I[2]], "literature_example", location="Scheme 1, cited prior work",
                 name=f"C60 cascade reaction cited from {source_name}", transformation="C60 cycloaddition")
    r["participants"] = [participant("substrate", f"{DOC}:mol:{sub}", "reactant", [E[2]]),
                         participant("c60", f"{DOC}:mol:C60", "reactant", [E[2]]),
                         participant("product", f"{DOC}:mol:{prod}", "product", [E[2]], product_role="reported cited product")]
    r["notes"] = [f("Cited literature example; not a current-paper experiment.",
                     "Cited literature example; not a current-paper experiment.", [E[2]])]
    classify(r, "prior_cycloaddition")
    R.append(r)

r = reaction("generic_scheme", [E[2], I[2]], "generic_scheme", location="Scheme 1, This work",
             name="Rh-catalyzed bisallene cycloisomerization / C60 Diels-Alder cascade",
             transformation="cascade cycloisomerization and [4+2] cycloaddition")
r["participants"] = [participant("generic_1", f"{DOC}:mol:generic_1", "reactant", [E[2]]),
                     participant("c60", f"{DOC}:mol:C60", "reactant", [E[2]]),
                     participant("generic_3", f"{DOC}:mol:generic_3", "product", [E[2]], product_role="generic product family")]
r["notes"] = [f("General R-group scheme; no experiment-specific yield.", "General R-group scheme; no experiment-specific yield.", [E[2]])]
classify(r, "cascade")
R.append(r)


ligands = {1: "DTBM", 2: "DTBM", 3: "DTBM", 4: "DTBM", 5: "DTBM", 6: "BINAP", 7: "H8BINAP", 8: "TolBINAP"}
temps = {1: 180, 2: 180, 3: 120, 4: 90, 5: 60, 6: 90, 7: 90, 8: 90}
times = {1: 30, 2: 15, 3: 15, 4: 20, 5: 30, 6: 240, 7: 240, 8: 240}
yields = {1: q(0, "%"), 2: q(5, "%", qualifier="lt", raw_text="<5%"),
          3: q(5, "%", qualifier="lt", raw_text="<5%"), 4: q(56, "%"),
          5: q(0, "%"), 6: q(42, "%"), 7: q(42, "%"), 8: q(20, "%")}

for i in range(1, 9):
    ev = [E[2], I[2], E[3], I[3]]
    if i == 4: ev += [E[4], I[4], E[6], I[6]]
    r = reaction(f"table1_entry{i}", ev, categories=("optimization", "substrate_scope") if i == 4 else ("optimization",),
                 location=f"Table 1, entry {i}" + ("; Figure 2, 3a; General Procedure" if i == 4 else ""),
                 name="Rhodium(I)-catalyzed cycloaddition of bisallene 1a with fullerene C60",
                 transformation="Rh-catalyzed cascade cycloisomerization / Diels-Alder cycloaddition", schemes=[GEN])
    one_ev = [E[2]]
    ratio = 0.5 if i == 1 else 1.1
    solvent_key = "toluene" if i == 5 else "oDCB"
    p1_amount = f("0.06 mmol", q(0.06, "mmol"), one_ev)
    c60_eq = f("C60:1a = 1:2", q(0.5, "equiv"), one_ev) if i == 1 else f("C60:1a = 1.1:1", q(1.1, "equiv"), one_ev)
    if i == 4:
        p1_amount = f(["0.06 mmol (Table 1 footnote)", "0.061 mmol (General Procedure)"], None, ev,
            status="conflict", alternatives=[
                f("0.06 mmol", q(0.06, "mmol"), [E[2]]), f("0.061 mmol", q(0.061, "mmol"), [E[6]])])
        c60_eq = f(["C60:1a = 1.1:1", "0.064 mmol C60 / 0.061 mmol 1a"], None, ev,
            status="conflict", alternatives=[
                f("C60:1a = 1.1:1", q(1.1, "equiv"), [E[2]]),
                derived(q(round(0.064/0.061, 4), "equiv"), [E[6]], "Calculated from the two reported mmol charges.",
                        raw="0.064 mmol / 0.061 mmol")])
    r["participants"] = [
        participant("bisallene_1a", f"{DOC}:mol:1a", "reactant", ev, amount=p1_amount,
                    mass=f("17.0 mg", q(17.0, "mg"), [E[6]]) if i == 4 else None),
        participant("c60", f"{DOC}:mol:C60", "reactant", ev, equivalents=c60_eq,
                    amount=f("0.064 mmol", q(0.064, "mmol"), [E[6]]) if i == 4 else None,
                    mass=f("46.3 mg", q(46.3, "mg"), [E[6]]) if i == 4 else None,
                    reference=f"{DOC}:mol:1a"),
        participant("rh", f"{DOC}:mol:Rh", "catalyst", ev,
                    mol_percent=f("10 mol%", q(10, "mol%"), one_ev), reference=f"{DOC}:mol:1a",
                    amount=f("0.006 mmol", q(0.006, "mmol"), [E[6]]) if i == 4 else None,
                    mass=f("2.2 mg", q(2.2, "mg"), [E[6]]) if i == 4 else None),
        participant("ligand", f"{DOC}:mol:{ligands[i]}", "ligand", ev,
                    mol_percent=f("10 mol%", q(10, "mol%"), one_ev), reference=f"{DOC}:mol:1a",
                    amount=f("0.006 mmol", q(0.006, "mmol"), [E[6]]) if i == 4 else None,
                    mass=f("7.5 mg", q(7.5, "mg"), [E[6]]) if i == 4 else None),
        participant("hydrogen", f"{DOC}:mol:H2", "reagent", ev),
        participant("ch2cl2", f"{DOC}:mol:CH2Cl2", "solvent", ev,
                    volume=f("4 mL", q(4, "mL"), [E[6]]) if i == 4 else None),
        participant("solvent", f"{DOC}:mol:{solvent_key}", "solvent", ev),
        participant("product_3a", f"{DOC}:mol:3a", "product", ev, product_role="desired cycloadduct"),
    ]
    if i in (1, 5):
        r["participants"][-1]["product_role"] = f(
            "intended desired product; not observed", "intended desired product; not observed", one_ev)
    r["reaction_scale"] = f(["0.06 mmol 1a", "0.061 mmol 1a"], None, ev, status="conflict",
        alternatives=[f("0.06 mmol 1a", {"reference_molecule_id": f"{DOC}:mol:1a", "amount": q(0.06, "mmol")}, [E[2]]),
                      f("0.061 mmol 1a", {"reference_molecule_id": f"{DOC}:mol:1a", "amount": q(0.061, "mmol")}, [E[6]])]) if i == 4 else \
        f("0.06 mmol of 1a", {"reference_molecule_id": f"{DOC}:mol:1a", "amount": q(0.06, "mmol")}, one_ev)
    r["solvent_conditions"]["composition"] = f("toluene" if i == 5 else "o-DCB", ["toluene"] if i == 5 else ["o-DCB"], one_ev)
    if i != 5:
        r["solvent_conditions"]["concentration"] = f("1.4 mM", q(1.4, "mM"), one_ev)
        r["solvent_conditions"]["concentration_basis"] = f("reaction in o-DCB (1.4 mM)", "basis not explicitly stated", one_ev)
    r["steps"] = [
        step("catalyst_activation", one_ev if i != 4 else [E[2], E[6]], ["rh", "ligand", "hydrogen", "ch2cl2"],
             operation="treat Rh/ligand solution with hydrogen before substrate addition",
             temperature="r.t.", time=q(30, "min") if i == 4 else None, atmosphere="hydrogen"),
        step("reaction", one_ev, ["bisallene_1a", "c60", "rh", "ligand", "solvent", "product_3a"],
             operation="heat/stir reaction mixture", temperature=q(temps[i], "°C"), time=q(times[i], "min")),
    ]
    r["procedure"]["pretreatment"] = f("10 mol% Rh/ligand mixture treated with hydrogen in CH2Cl2 before substrate addition",
                                        "Rh/ligand hydrogen pretreatment before substrate addition", one_ev)
    r["outcomes"] = [outcome("yield_3a", "yield", [f"{DOC}:mol:3a"], one_ev, yields[i],
                              role="desired product", yield_type="isolated yield" if i == 4 else None,
                              basis="isolated after column chromatography" if i == 4 else None)]
    if i == 1:
        r["outcomes"].append(outcome("starting_material_recovery", "recovery",
            [f"{DOC}:mol:1a", f"{DOC}:mol:C60"], [E[2]], "not recovered", scope="both starting materials"))
        r["notes"] = [f("No starting materials were recovered; complete retention at the TLC baseline suggested polyaddition.",
                         "No starting materials were recovered; complete retention at the TLC baseline suggested polyaddition.", [E[2]])]
    elif i == 5:
        r["outcomes"].append(outcome("no_reactivity", "no_reaction", [], [E[3]], True,
                                      scope="no reactivity observed"))
    else:
        r["notes"] = []
    if i == 4:
        r["participants"] += [
            participant("silica", f"{DOC}:mol:silica", "purification_agent", [E[6]], stage="purification"),
            participant("cs2", f"{DOC}:mol:CS2", "purification_agent", [E[6]], stage="purification"),
            participant("toluene_eluent", f"{DOC}:mol:toluene", "purification_agent", [E[6]], stage="purification"),
        ]
        r["steps"] = [
            step("prepare_catalyst_solution", [E[6]], ["rh", "ligand", "ch2cl2"], operation="prepare anhydrous CH2Cl2 solution under nitrogen", atmosphere="nitrogen"),
            step("hydrogen_activation", [E[6]], ["rh", "ligand", "hydrogen", "ch2cl2"], operation="bubble hydrogen through catalyst solution", temperature="r.t.", time=q(30, "min"), atmosphere="hydrogen"),
            step("concentrate_catalyst", [E[6]], ["rh", "ligand"], operation="concentrate catalyst solution to dryness"),
            step("add_to_preheated_substrates", [E[6]], ["bisallene_1a", "c60", "rh", "ligand", "solvent"], operation="redissolve catalyst in anhydrous o-DCB and introduce into preheated C60/1a solution", temperature=q(90, "°C"), atmosphere="nitrogen"),
            step("reaction", [E[2], E[3], E[6]], ["bisallene_1a", "c60", "rh", "ligand", "solvent", "product_3a"], operation="stir with TLC monitoring", temperature=q(90, "°C"), time=q(20, "min"), atmosphere="nitrogen"),
            step("cool", [E[6]], ["product_3a"], operation="allow reaction mixture to cool", temperature="room temperature"),
            step("column_cs2", [E[6]], ["product_3a", "c60", "silica", "cs2"], stage="purification", operation="silica-gel chromatography with CS2 to elute unreacted C60"),
            step("column_toluene", [E[6]], ["product_3a", "silica", "toluene_eluent"], stage="purification", operation="further elution with toluene to furnish 3a"),
        ]
        r["special_conditions"]["reactor"] = f("10 mL capped vial", "10 mL capped vial", [E[6]])
        r["procedure"]["original_text"] = f("General Procedure for the Synthesis of Fused Dihydroazepine Derivatives of C60 (3a)",
                                             "General procedure for 3a", [E[6]])
        r["procedure"]["addition_order"] = f("activate Rh/ligand with H2; dry; redissolve in o-DCB; add to preheated C60/1a solution",
                                              ["activate catalyst", "dry", "redissolve", "add to preheated substrates"], [E[6]])
        r["outcomes"] += [
            outcome("amount_3a", "product_amount", [f"{DOC}:mol:3a"], [E[6]], q(34.1, "mg"), role="isolated product"),
            outcome("recovered_c60", "recovery", [f"{DOC}:mol:C60"], [E[6]], "unreacted fullerene obtained; amount not reported", scope="amount not reported"),
        ]
        r["notes"] = [
            f("Table 1 entry 4, Figure 2 product 3a, and the page-6 General Procedure describe the same 56% experiment and are merged.",
              "Same experiment repeated across Table 1, Figure 2, and the General Procedure.", ev),
            f("Table 1 reports 0.06 mmol 1a and C60:1a = 1.1:1; the General Procedure reports 0.061 mmol 1a and 0.064 mmol C60 (calculated ratio 1.049:1). Both are retained without selecting one.",
              "Source quantity discrepancy retained as conflict.", ev),
        ]
        r["review"]["unresolved"] += ["Resolve the rounded Table 1 scale/ratio versus the General Procedure charges; neither value was selected."]
    classify(r, "not_applicable" if i in (1, 5) else "cascade")
    R.append(r)

# Prose-only reduced catalyst loading experiment.
r = reaction("optimization_5molpct", [E[3]], categories=("optimization",), location="Page 3 prose after Table 1",
             name="Attempted reduction of Rh/diphosphine catalyst loading to 5 mol%",
             transformation="Rh-catalyzed cascade cycloaddition", schemes=[GEN])
r["participants"] = [
    participant("bisallene_1a", f"{DOC}:mol:1a", "reactant", [E[3]]), participant("c60", f"{DOC}:mol:C60", "reactant", [E[3]]),
    participant("rh", f"{DOC}:mol:Rh", "catalyst", [E[3]], mol_percent=f("5 mol%", q(5, "mol%"), [E[3]]), reference=f"{DOC}:mol:1a"),
    participant("ligand", f"{DOC}:mol:DTBM", "ligand", [E[3]]), participant("product_3a", f"{DOC}:mol:3a", "product", [E[3]], product_role="desired cycloadduct")]
r["steps"] = [step("reaction", [E[3]], [p["participant_id"] for p in r["participants"]], operation="attempt reaction at reduced catalyst loading")]
r["outcomes"] = [outcome("qualitative_yield", "yield", [f"{DOC}:mol:3a"], [E[3]], "detrimental in terms of reaction yield", role="desired product")]
r["notes"] = [f("Temperature, time, ligand loading, and numerical yield were not stated for this catalyst-loading attempt.",
                 "Temperature, time, ligand loading, and numerical yield were not stated for this catalyst-loading attempt.", [E[3]])]
r["review"]["unresolved"] += ["The prose gives no numerical yield or complete conditions for the 5 mol% experiment."]
classify(r, "cascade")
R.append(r)

# C60 substrate scope; 3a is already merged into Table 1 entry 4.
scope_yields = {"b":56,"c":35,"d":37,"e":35,"f":25,"g":49,"h":36,"i":43,"j":53,"k":33,"l":36,"m":27,"n":32,"o":39,"p":22}
scope_masses = {"b":33.7,"c":21.6,"d":24.7,"e":19.8,"f":16.6,"g":28.3,"h":20.3,"i":27.8,"j":37.2,"k":23.3,"l":20.8,"m":19.7,"n":17.7,"o":22.0,"p":13.1}
scope_characterization_page = {"b":6,"c":6,"d":7,"e":7,"f":7,"g":7,"h":7,"i":7,"j":7,"k":7,"l":8,"m":8,"n":8,"o":8,"p":8}
for x, y in scope_yields.items():
    lab1, lab3 = f"1{x}", f"3{x}"
    char_page = scope_characterization_page[x]
    ev = [E[4], I[4], E[2], I[2], E[char_page], I[char_page]]
    r = reaction(f"scope_{lab3}", ev, categories=("substrate_scope",), location=f"Figure {'2' if x <= 'l' else '3'}, {lab3}",
                 name=f"Cycloaddition of bisallene {lab1} with C60", transformation="Rh-catalyzed cascade cycloaddition", schemes=[GEN])
    time = 120 if x in "nop" else 20
    r["participants"] = [
        participant("bisallene", f"{DOC}:mol:{lab1}", "reactant", ev),
        participant("c60", f"{DOC}:mol:C60", "reactant", ev, equivalents=f("optimized C60:bisallene ratio 1.1:1", q(1.1, "equiv"), [E[2]], method="inherited_general_condition"), reference=f"{DOC}:mol:{lab1}"),
        participant("rh", f"{DOC}:mol:Rh", "catalyst", ev, mol_percent=f("optimized 10 mol%", q(10, "mol%"), [E[2]], method="inherited_general_condition"), reference=f"{DOC}:mol:{lab1}"),
        participant("ligand", f"{DOC}:mol:DTBM", "ligand", ev, mol_percent=f("optimized 10 mol%", q(10, "mol%"), [E[2]], method="inherited_general_condition"), reference=f"{DOC}:mol:{lab1}"),
        participant("solvent", f"{DOC}:mol:oDCB", "solvent", ev),
        participant("product", f"{DOC}:mol:{lab3}", "product", ev, product_role="reported cycloadduct"),
    ]
    r["solvent_conditions"]["composition"] = f("optimized o-DCB conditions", ["o-DCB"], [E[2]], method="inherited_general_condition")
    r["solvent_conditions"]["concentration"] = f("1.4 mM optimized condition", q(1.4, "mM"), [E[2]], method="inherited_general_condition")
    r["solvent_conditions"]["concentration_basis"] = f("basis not explicitly stated", "basis not explicitly stated", [E[2]], method="inherited_general_condition")
    r["steps"] = [step("reaction", ev, [p["participant_id"] for p in r["participants"]], operation="apply optimized cascade conditions",
                        temperature=q(90, "°C"), time=q(time, "min"), inherited=True)]
    r["outcomes"] = [outcome("yield", "yield", [f"{DOC}:mol:{lab3}"], [E[4]], q(y, "%"), role="reported product")]
    r["outcomes"].append(outcome("product_amount", "product_amount", [f"{DOC}:mol:{lab3}"],
                                 [E[char_page], I[char_page]], q(scope_masses[x], "mg"), role="reported product mass"))
    r["procedure"]["general_procedure_reference"] = f("optimized conditions described for 3a", f"{DOC}:table1_entry4", [E[2]], method="inherited_general_condition")
    r["review"]["unresolved"] += ["Yield measurement basis is not stated in Figure 2/3; no isolated/NMR basis was assumed.",
                                      "Absolute substrate charge and experiment-specific workup are not reported in the supplied main-text pages for this scope member."]
    classify(r, "cascade")
    R.append(r)

# Oxygen-tethered failure.
r = reaction("scope_oxygen_tethered_no_adduct", [E[4], I[4]], categories=("substrate_scope",), location="Page 4 prose below Figure 2",
             name="Attempted reaction of oxygen-tethered 1,5-bisallene with C60", transformation="attempted cascade cycloaddition", schemes=[GEN])
r["participants"] = [participant("bisallene", f"{DOC}:mol:1O", "reactant", [E[4]]), participant("c60", f"{DOC}:mol:C60", "reactant", [E[4]])]
r["steps"] = [step("reaction", [E[4]], ["bisallene", "c60"], operation="set up the reaction; detailed conditions not restated")]
r["outcomes"] = [outcome("no_functionalized_fullerene", "no_reaction", [], [E[4]], True, scope="no hint of functionalized fullerene detected")]
r["review"]["unresolved"] += ["Exact substrate label, complete conditions, detection method, conversion, recovery, and byproducts are unreported."]
classify(r, "not_applicable")
R.append(r)

# Light/air exposure experiment.
r = reaction("light_air_exposure_3a", [E[3], I[3]], categories=("mechanistic_study",), location="Page 3 prose, light/oxygen exposure",
             name="Light and air exposure of cycloadduct 3a", transformation="nonspecific photooxidation/degradation study")
r["participants"] = [participant("substrate_3a", f"{DOC}:mol:3a", "reactant", [E[3]]),
                     participant("chloroform", f"{DOC}:mol:CHCl3", "solvent", [E[3]])]
r["solvent_conditions"]["composition"] = f("CHCl3 solution", ["CHCl3"], [E[3]])
r["steps"] = [step("exposure", [E[3]], ["substrate_3a", "chloroform"], operation="stir under air and light exposure", time=q(16, "h"), atmosphere="air")]
r["special_conditions"]["light_source"] = f("light exposure", "ambient light exposure; source/power not reported", [E[3]])
r["outcomes"] = [outcome("consumption_3a", "conversion", [f"{DOC}:mol:3a"], [E[3]], "complete consumption", basis="TLC monitoring", scope="3a consumption")]
r["notes"] = [f("A new retained brown TLC spot and major [M(3a)+O2+Na]+ HRMS peak were observed, but purification/NMR gave a complex mixture; no defined oxidation product was identified.",
                 "Undefined oxygenated complex mixture; author hypothesized nonspecific singlet-oxygen chemistry.", [E[3]])]
r["review"]["unresolved"] += ["No defined product identity, yield, temperature, light source, oxygen pressure, or composition was established."]
classify(r, "unresolved")
R.append(r)

# C70 cycloaddition.  The page-8 experimental procedure is the same experiment as Scheme 2.
c70_ev = [E[4], I[4], E[5], I[5], E[8], I[8], E[9], I[9]]
r = reaction("c70_cycloaddition_4a", c70_ev, categories=("substrate_scope",),
             location="Scheme 2; Experimental Procedure for 4a (page 8)",
             name="Cycloaddition of bisallene 1a with C70", transformation="Rh-catalyzed cascade cycloaddition", schemes=[GEN])
r["participants"] = [
    participant("bisallene_1a", f"{DOC}:mol:1a", "reactant", [E[8], I[8]], mass=f("31.1 mg", q(31.1, "mg"), [E[8], I[8]]), amount=f("0.11 mmol", q(0.11, "mmol"), [E[8], I[8]])),
    participant("c70", f"{DOC}:mol:C70", "reactant", [E[8], I[8]], mass=f("100 mg", q(100, "mg"), [E[8], I[8]]), amount=f("0.12 mmol", q(0.12, "mmol"), [E[8], I[8]])),
    participant("rh", f"{DOC}:mol:Rh", "catalyst", [E[8], I[8]], mass=f("5.3 mg", q(5.3, "mg"), [E[8], I[8]]), amount=f("0.01 mmol", q(0.01, "mmol"), [E[8], I[8]]), mol_percent=f("0.01 mmol relative to 0.11 mmol 1a", q(9.09, "mol%", raw_text="0.01/0.11 mmol"), [E[8], I[8]], method="inference", derivation="Calculated from the reported mmol charges; the source does not print mol%."), reference=f"{DOC}:mol:1a"),
    participant("ligand", f"{DOC}:mol:DTBM", "ligand", [E[8], I[8]], mass=f("13.5 mg", q(13.5, "mg"), [E[8], I[8]]), amount=f("0.01 mmol", q(0.01, "mmol"), [E[8], I[8]]), mol_percent=f("0.01 mmol relative to 0.11 mmol 1a", q(9.09, "mol%", raw_text="0.01/0.11 mmol"), [E[8], I[8]], method="inference", derivation="Calculated from the reported mmol charges; the source does not print mol%."), reference=f"{DOC}:mol:1a"),
    participant("hydrogen", f"{DOC}:mol:H2", "reagent", [E[8], I[8]]),
    participant("ch2cl2", f"{DOC}:mol:CH2Cl2", "solvent", [E[8], I[8]], volume=f("4 mL", q(4, "mL"), [E[8], I[8]])),
    participant("odcb", f"{DOC}:mol:oDCB", "solvent", [E[8], I[8]]),
    participant("alpha_product", f"{DOC}:mol:alpha4a", "product", c70_ev, product_role="major regioisomer"),
    participant("beta_product", f"{DOC}:mol:beta4a", "product", c70_ev, product_role="minor regioisomer"),
    participant("silica", f"{DOC}:mol:silica", "purification_agent", [E[8], I[8]], stage="purification"),
    participant("cs2", f"{DOC}:mol:CS2", "purification_agent", [E[8], I[8]], stage="purification"),
    participant("toluene_eluent", f"{DOC}:mol:toluene", "purification_agent", [E[8], I[8]], stage="purification"),
]
r["reaction_scale"] = f("0.11 mmol of bisallene 1a", {"reference_molecule_id": f"{DOC}:mol:1a", "amount": q(0.11, "mmol")}, [E[8], I[8]])
r["solvent_conditions"]["composition"] = f("anhydrous o-DCB; catalyst initially in 4 mL anhydrous CH2Cl2", ["anhydrous o-DCB", "4 mL anhydrous CH2Cl2 catalyst carrier"], [E[8], I[8]])
r["solvent_conditions"]["concentration"] = f("o-DCB (1.4 mM)", q(1.4, "mM"), [E[8], I[8]])
r["solvent_conditions"]["concentration_basis"] = f("solution of C70 and 1a in anhydrous o-DCB (1.4 mM)", "exact concentration basis not explicitly identified", [E[8], I[8]])
r["steps"] = [
    step("prepare_catalyst_solution", [E[8], I[8]], ["rh", "ligand", "ch2cl2"], operation="prepare anhydrous CH2Cl2 solution in a capped vial", atmosphere="nitrogen"),
    step("catalyst_activation", [E[8], I[8]], ["rh", "ligand", "hydrogen", "ch2cl2"], operation="bubble hydrogen gas into the catalyst solution", time=q(30, "min"), atmosphere="hydrogen"),
    step("concentrate_catalyst", [E[8], I[8]], ["rh", "ligand"], operation="concentrate catalyst solution to dryness"),
    step("add_to_preheated_substrates", [E[8], I[8]], ["bisallene_1a", "c70", "rh", "ligand", "odcb"], operation="redissolve catalyst in anhydrous o-DCB and introduce it into the C70/1a solution preheated to 90 °C", temperature=q(90, "°C"), atmosphere="nitrogen"),
    step("reaction", [E[8], I[8]], ["bisallene_1a", "c70", "rh", "ligand", "odcb", "alpha_product", "beta_product"], operation="stir with TLC monitoring", temperature=q(90, "°C"), time=q(40, "min"), atmosphere="nitrogen"),
    step("cool", [E[8], I[8]], ["alpha_product", "beta_product"], operation="allow mixture to cool", temperature="room temperature"),
    step("column_cs2", [E[8], I[8]], ["alpha_product", "beta_product", "c70", "silica", "cs2"], stage="purification", operation="silica-gel chromatography with CS2 to elute unreacted C70"),
    step("column_toluene", [E[8], I[8]], ["alpha_product", "beta_product", "silica", "toluene_eluent"], stage="purification", operation="further elution with toluene to furnish the 4a regioisomer mixture"),
]
r["special_conditions"]["reactor"] = f("10 mL capped vial", "10 mL capped vial", [E[8], I[8]])
r["procedure"]["original_text"] = f("Experimental Procedure for the Synthesis of Fused Dihydroazepine Derivatives of C70 (4a)", "Experimental procedure for 4a", [E[8], I[8]])
r["procedure"]["addition_order"] = f("prepare Rh/ligand solution; bubble H2; dry; redissolve in o-DCB; add to preheated C70/1a solution", ["prepare catalyst solution", "H2 activation", "dry", "redissolve", "add to preheated substrates"], [E[8], I[8]])
r["procedure"]["pretreatment"] = f("Hydrogen gas was bubbled into the catalyst solution for 30 min", "H2 activation of Rh/ligand solution for 30 min", [E[8], I[8]])
r["outcomes"] = [
    outcome("combined_yield", "yield", [f"{DOC}:mol:alpha4a", f"{DOC}:mol:beta4a"], [E[8], I[8]], q(46, "%"), role="inseparable regioisomer mixture", yield_type="isolated yield", basis="silica-gel chromatography followed by toluene elution", scope="alpha-4a + beta-4a; traces of other regioisomers"),
    outcome("product_amount", "product_amount", [f"{DOC}:mol:alpha4a", f"{DOC}:mol:beta4a"], [E[8], I[8]], q(57.4, "mg"), role="isolated regioisomer mixture"),
    outcome("regioisomer_ratio", "selectivity", [f"{DOC}:mol:alpha4a", f"{DOC}:mol:beta4a"], [E[4], E[5]], "3:1", basis="2D NMR spectra and HPLC chromatograms", scope={"metric":"regioisomer ratio","order":["alpha-4a","beta-4a"],"major":"alpha-4a"}),
]
r["notes"] = [
    f("The two major regioisomers could not be separated; UV-vis supported an alpha adduct as the major product, and traces of other regioisomers were present.",
      "Inseparable alpha/beta mixture with traces of other regioisomers.", [E[4], E[5]]),
    f("Scheme 2 and the page-8 Experimental Procedure report the same 46%, 3:1 experiment and are merged rather than counted twice.",
      "Repeated source locations merged into one experiment.", [E[5], I[5], E[8], I[8]], method="inference", derivation="Same substrate/product labels, conditions, yield, and ratio identify one repeated experiment."),
]
r["review"]["unresolved"] += ["The total o-DCB volume/concentration denominator and the identities/amounts of trace regioisomers are not reported."]
classify(r, "cascade")
R.append(r)

# Selective reductions.
for suffix, start_ids, prod_ids, y, ratio in [
    ("reduction_3a_to_5a", ["3a"], ["5a"], 82, None),
    ("reduction_4a_to_6a", ["alpha4a", "beta4a"], ["alpha6a", "beta6a"], 69, "3:1"),
]:
    reduction_ev = [E[5], I[5], E[9], I[9]]
    r = reaction(suffix, reduction_ev, categories=("derivatization",),
                 location="Scheme 3; page-9 General Procedure/compound entry",
                 name="Selective reduction of the double bond vicinal to the sulfonamide group",
                 transformation="selective enamine/double-bond reduction")
    r["participants"] = []
    for j, mid in enumerate(start_ids):
        r["participants"].append(participant(f"start_{j+1}", f"{DOC}:mol:{mid}", "reactant", reduction_ev))
    if suffix == "reduction_3a_to_5a":
        r["participants"][0]["charges"]["mass"] = f("100.0 mg", q(100.0, "mg"), [E[9], I[9]])
        r["participants"][0]["charges"]["amount"] = f("0.10 mmol", q(0.10, "mmol"), [E[9], I[9]])
        r["participants"] += [
            participant("et3sih", f"{DOC}:mol:Et3SiH", "reagent", [E[5], E[9], I[9]], volume=f("0.15 mL", q(0.15, "mL"), [E[9], I[9]]), amount=f("0.9 mmol", q(0.9, "mmol"), [E[9], I[9]]), equivalents=f("excess Et3SiH", "excess", [E[5]]), reference=f"{DOC}:mol:3a"),
            participant("bf3", f"{DOC}:mol:BF3Et2O", "reagent", [E[5], E[9], I[9]], volume=f("0.1 mL", q(0.1, "mL"), [E[9], I[9]]), amount=f("0.81 mmol", q(0.81, "mmol"), [E[9], I[9]]), reference=f"{DOC}:mol:3a"),
            participant("ch2cl2", f"{DOC}:mol:CH2Cl2", "solvent", [E[5], E[9], I[9]], volume=f("10 mL", q(10, "mL"), [E[9], I[9]])),
            participant("silica", f"{DOC}:mol:silica", "purification_agent", [E[9], I[9]], stage="purification"),
            participant("toluene_eluent", f"{DOC}:mol:toluene", "purification_agent", [E[9], I[9]], stage="purification"),
        ]
    else:
        r["participants"] += [
            participant("et3sih", f"{DOC}:mol:Et3SiH", "reagent", [E[5]], equivalents=f("excess Et3SiH", "excess", [E[5]], method="inherited_general_condition"), reference=f"{DOC}:mol:alpha4a"),
            participant("bf3", f"{DOC}:mol:BF3Et2O", "reagent", [E[5]]),
            participant("ch2cl2", f"{DOC}:mol:CH2Cl2", "solvent", [E[5]]),
        ]
    for j, mid in enumerate(prod_ids):
        r["participants"].append(participant(f"product_{j+1}", f"{DOC}:mol:{mid}", "product", reduction_ev, product_role="reduced product"))
    r["solvent_conditions"]["composition"] = f("CH2Cl2", ["CH2Cl2"], [E[5], E[9]] if suffix == "reduction_3a_to_5a" else [E[5]])
    ids = [p["participant_id"] for p in r["participants"]]
    if suffix == "reduction_3a_to_5a":
        r["reaction_scale"] = f("0.10 mmol of 3a", {"reference_molecule_id": f"{DOC}:mol:3a", "amount": q(0.10, "mmol")}, [E[9], I[9]])
        r["solvent_conditions"]["total_volume"] = f("10 mL anhydrous CH2Cl2", q(10, "mL"), [E[9], I[9]])
        r["steps"] = [
            step("prepare_and_cool", [E[9], I[9]], ["start_1", "ch2cl2"], operation="stir 3a in anhydrous CH2Cl2 and cool", temperature=q(0, "°C"), atmosphere="nitrogen"),
            step("low_temperature_addition", [E[9], I[9]], ["start_1", "ch2cl2", "et3sih", "bf3"], operation="add Et3SiH and BF3·Et2O", temperature=q(0, "°C"), atmosphere="nitrogen"),
            step("room_temperature_hold", [E[9], I[9]], ids, operation="stir resulting mixture with TLC monitoring", temperature="room temperature", time=q(5, "h"), atmosphere="nitrogen"),
            step("column", [E[9], I[9]], ["product_1", "silica", "toluene_eluent"], stage="purification", operation="subject directly to silica-gel column chromatography using toluene as eluent"),
        ]
        r["special_conditions"]["reactor"] = f("25 mL round-bottom flask equipped with magnetic stirrer", "25 mL round-bottom flask with magnetic stirrer", [E[9], I[9]])
        r["procedure"]["original_text"] = f("General Procedure for the Reduction of Fused Dihydroazepine Derivatives (5a)", "General reduction procedure for 5a", [E[9], I[9]])
        r["procedure"]["addition_order"] = f("prepare 3a solution; cool to 0 °C; add Et3SiH and BF3·Et2O; stir 5 h at room temperature; chromatograph", ["prepare solution", "cool", "add reagents", "room-temperature hold", "chromatography"], [E[9], I[9]])
        r["outcomes"] = [
            outcome("yield", "yield", [f"{DOC}:mol:5a"], [E[9], I[9]], q(y, "%"), role="reduced product", yield_type="isolated yield", basis="silica-gel column chromatography with toluene eluent"),
            outcome("product_amount", "product_amount", [f"{DOC}:mol:5a"], [E[9], I[9]], q(82.1, "mg"), role="isolated reduced product"),
        ]
        r["review"]["unresolved"] += ["The durations of cooling, reagent addition, and warming from 0 °C to room temperature are not reported; no separate quench/workup is described before direct chromatography."]
    else:
        r["steps"] = [
            step("low_temperature_start", [E[5], I[5]], ids, operation="apply the same reduction protocol to the C70 4a regioisomer mixture with Et3SiH and BF3·Et2O in CH2Cl2", temperature=q(0, "°C"), inherited=True),
            step("warm_to_room_temperature", [E[5], I[5]], ids, operation="warm from 0 °C to room temperature; exact hold time is not stated for 4a", temperature="room temperature", inherited=True),
        ]
        r["procedure"]["general_procedure_reference"] = f("The same reaction protocol was then applied to the corresponding isomeric mixture of C70 adducts", f"{DOC}:reduction_3a_to_5a", [E[5], I[5]], method="inherited_general_condition")
        r["procedure"]["modifications"] = f("C70 adduct mixture alpha-4a/beta-4a used instead of 3a", "starting material changed to the C70 4a regioisomer mixture", [E[5], I[5]], method="inherited_general_condition")
        r["outcomes"] = [
            outcome("yield", "yield", [f"{DOC}:mol:{x}" for x in prod_ids], [E[9], I[9]], q(y, "%"), role="combined reduced regioisomer mixture"),
            outcome("product_amount", "product_amount", [f"{DOC}:mol:{x}" for x in prod_ids], [E[9], I[9]], q(42.0, "mg"), role="reported reduced regioisomer mixture"),
        ]
    if ratio:
        r["outcomes"].append(outcome("regioisomer_ratio", "selectivity", [f"{DOC}:mol:{x}" for x in prod_ids], [E[5]], ratio,
            basis="2D NMR spectra", scope={"metric":"regioisomer ratio","order":["alpha-6a","beta-6a"],"major":"alpha-6a"}))
        r["notes"] = [f(
            "Scheme 3 labels the starting C70 mixture as α-4a:β-4a (3:1); the reduced products are separately reported as α-6a:β-6a (3:1).",
            {"starting_mixture": {"metric": "regioisomer ratio", "order": ["alpha-4a", "beta-4a"], "value": "3:1"},
             "product_mixture": {"metric": "regioisomer ratio", "order": ["alpha-6a", "beta-6a"], "value": "3:1"}},
            [E[5], I[5]])]
        r["review"]["unresolved"] += ["For the C70 reduction, substrate/reagent quantities, exact reaction time, warming duration, workup/purification, and yield measurement basis are not separately reported; 'same reaction protocol' was retained without copying 5a charges."]
    classify(r, "hydrogenation")
    R.append(r)

# Source context: non-synthetic observations and explicit boundaries.
template["source_context"] = [
    {"context_id": f"{DOC}:context:metadata", "evidence_ids": [E[1], I[1]],
     "kind": "bibliographic_metadata", "title": "Synthesis of Fused Dihydroazepine Derivatives of Fullerenes by a Rh-Catalyzed Cascade Process",
     "doi": "10.1002/adsc.202100644", "note": "The main text advertises online SI, but no SI file was explicitly supplied or associated."},
    {"context_id": f"{DOC}:context:3a_characterization", "evidence_ids": [E[3], I[3], E[6], I[6]], "kind": "characterization",
     "note": "3a formula/NMR/UV-vis/HRMS and compound characterization are linked to its synthesis; they are not separate synthesis experiments."},
    {"context_id": f"{DOC}:context:photooxidation_hrms", "evidence_ids": [E[3], I[3]], "kind": "characterization",
     "note": "The [M(3a)+O2+Na]+ HRMS observation supports an undefined complex oxidation mixture and was not made a separate synthesis record."},
    {"context_id": f"{DOC}:context:figure1_computation", "evidence_ids": [E[3], I[3]], "kind": "computational_mechanism",
     "note": "Figure 1 is a computed Gibbs-energy profile, not an executed reaction experiment."},
    {"context_id": f"{DOC}:context:figure4_computation", "evidence_ids": [E[5], I[5]], "kind": "computational_regioselectivity",
     "note": "Figure 4 contains computed C70 pathways/activation-strain analysis, not separate synthesis experiments."},
    {"context_id": f"{DOC}:context:electrochemistry", "evidence_ids": [E[6], I[6]], "kind": "physical_property",
     "note": "Tables 2-3 report cyclic-voltammetry and optical/electronic properties; no synthesis records were created from these tables."},
    {"context_id": f"{DOC}:context:scope_characterization_3c_3k", "evidence_ids": [E[7], I[7]], "kind": "characterization",
     "note": "Page 7 contains continuation/compound characterization for 3c and full characterization blocks for 3d-3k. Product masses were linked to the existing scope experiments; IR/NMR/UV-vis/HRMS and molecular weights were not converted into new reactions."},
    {"context_id": f"{DOC}:context:scope_characterization_3k_3p", "evidence_ids": [E[8], I[8]], "kind": "characterization",
     "note": "Page 8 contains continuation/compound characterization for 3k and characterization blocks for 3l-3p. Product masses were linked to the existing scope experiments; spectral data and molecular weights remain characterization context."},
    {"context_id": f"{DOC}:context:4a_5a_6a_characterization", "evidence_ids": [E[8], I[8], E[9], I[9]], "kind": "characterization",
     "note": "Compound-characterization blocks for the 4a, 5a, and 6a products support their existing reactions. The reported masses/yields were linked to those reactions; NMR/IR/UV-vis/HRMS and molecular weights remain non-synthetic context."},
    {"context_id": f"{DOC}:context:references_pages_9_10", "evidence_ids": [E[9], I[9], E[10], I[10]], "kind": "bibliographic_references",
     "note": "The right column of page 9 and all of page 10 are acknowledgements/references, not current-paper experiments; cited works were not promoted to new records from bibliography entries alone."},
    {"context_id": f"{DOC}:context:negative_spectral_inventory", "evidence_ids": [E[i] for i in range(1, 11)] + [I[i] for i in range(1, 11)], "kind": "coverage_boundary",
     "note": "Standalone spectral values, molecular weights, m/z values, and bibliography entries on pages 1-10 were not treated as reactions; product masses/yields were associated only with the already identified synthesis experiments."},
]

template["molecules"] = list(M.values())
template["reactions"] = R
dump(RUN / "submitted.json", template)

page_records = {
    1: [],
    2: [r["reaction_id"] for r in R if r["reaction_id"].endswith(("literature_itoh","literature_maeda","literature_wang","generic_scheme"))] +
       [f"{DOC}:table1_entry{i}" for i in range(1,9)],
    3: [f"{DOC}:table1_entry{i}" for i in range(1,9)] + [f"{DOC}:optimization_5molpct", f"{DOC}:light_air_exposure_3a"],
    4: [f"{DOC}:table1_entry4"] + [f"{DOC}:scope_3{x}" for x in "bcdefghijklmnop"] +
       [f"{DOC}:scope_oxygen_tethered_no_adduct", f"{DOC}:c70_cycloaddition_4a"],
    5: [f"{DOC}:c70_cycloaddition_4a", f"{DOC}:reduction_3a_to_5a", f"{DOC}:reduction_4a_to_6a"],
    6: [f"{DOC}:table1_entry4", f"{DOC}:scope_3b", f"{DOC}:scope_3c"],
}
regions = {
    1: [
        {"region":"title/abstract/metadata", "finding":"Bibliographic metadata, abstract, and SI availability notice; no current synthesis arrow or experiment.", "evidence_ids":[E[1], I[1]]},
        {"region":"introduction columns", "finding":"Background discussion only; cited reactions are not depicted as executable records on this page.", "evidence_ids":[E[1], I[1]]},
    ],
    2: [
        {"region":"Scheme 1, upper-left", "finding":"Three cited prior-work fullerene cascades retained as literature_example records.", "evidence_ids":[E[2], I[2]]},
        {"region":"Scheme 1, This work", "finding":"Current R-group reaction family retained as one generic_scheme, not an experiment-specific yield.", "evidence_ids":[E[2], I[2]]},
        {"region":"Table 1 and footnote", "finding":"Eight optimization rows retained independently; entry 4 is cross-page duplicate evidence for Figure 2 product 3a and the page-6 procedure.", "evidence_ids":[E[2], I[2]]},
        {"region":"results prose below Table 1", "finding":"Entry-1 non-recovery/polyaddition observation and entry-2/3 discussion associated with their table rows.", "evidence_ids":[E[2], I[2]]},
    ],
    3: [
        {"region":"left-column optimization discussion", "finding":"Table 1 rows linked; a separate qualitative 5 mol% catalyst-loading attempt was retained because its numerical yield/complete conditions are absent.", "evidence_ids":[E[3], I[3]]},
        {"region":"left-column 3a characterization", "finding":"NMR/UV-vis/HRMS characterize 3a and remain source_context, not extra synthesis records.", "evidence_ids":[E[3], I[3]]},
        {"region":"left-column light/air exposure", "finding":"Actual 16 h light/air exposure retained as a mechanistic study; product is an undefined complex oxidation mixture.", "evidence_ids":[E[3], I[3]]},
        {"region":"Figure 1/right-column computation", "finding":"DFT Gibbs-energy profile and mechanistic discussion retained as computational context, not an executed experiment.", "evidence_ids":[E[3], I[3]]},
    ],
    4: [
        {"region":"Figure 2", "finding":"Products 3a-l inventoried; 3a merged with Table 1 entry 4, while 3b-l are independent scope experiments.", "evidence_ids":[E[4], I[4]]},
        {"region":"Figure 3", "finding":"Products 3m-p inventoried as independent scope experiments; 1n-p retain the stated 2 h time.", "evidence_ids":[E[4], I[4]]},
        {"region":"prose below Figure 2", "finding":"O-tethered bisallene attempt retained as a real no-adduct experiment without inventing a product.", "evidence_ids":[E[4], I[4]]},
        {"region":"right-column C70 discussion", "finding":"C70 reaction identified as an inseparable alpha/beta regioisomer mixture and linked forward to Scheme 2.", "evidence_ids":[E[4], I[4]]},
    ],
    5: [
        {"region":"Scheme 2", "finding":"One C70 cycloaddition experiment gives alpha-4a/beta-4a in 46% combined yield and 3:1 ratio, with trace other regioisomers.", "evidence_ids":[E[5], I[5]]},
        {"region":"Scheme 3", "finding":"Two distinct reductions retained: 3a to 5a (82%) and the C70 4a mixture to 6a mixture (69%, 3:1).", "evidence_ids":[E[5], I[5]]},
        {"region":"Figure 4", "finding":"DFT energy/activation-strain plots are computational context, not synthesis experiments.", "evidence_ids":[E[5], I[5]]},
    ],
    6: [
        {"region":"right-column General Procedure for 3a", "finding":"Detailed 3a procedure merged into Table 1 entry 4, including activation, charges, staged operations, chromatography, 34.1 mg and 56% yield.", "evidence_ids":[E[6], I[6]]},
        {"region":"Tables 2-3", "finding":"Electrochemical and optical/electronic property measurements retained as source_context; no synthesis records created.", "evidence_ids":[E[6], I[6]]},
        {"region":"compound 3a-3c characterization", "finding":"Characterization remains context; 3b/3c isolated masses were linked to their existing scope records without creating duplicate reactions.", "evidence_ids":[E[6], I[6]]},
    ],
}
unresolved = {
    1: ["Online SI is advertised but was not supplied; no SI association was inferred."],
    2: ["Table 1 does not state a yield measurement basis for rows other than entry 4 as clarified by prose/procedure."],
    3: ["The 5 mol% catalyst-loading experiment lacks numerical yield and complete conditions; the oxidation product remains undefined."],
    4: ["Most compact fullerene-product structures were not converted to complete SMILES; scope-yield measurement bases are not stated."],
    5: ["Reduction reagent quantities, times, ramp duration, workup/purification, and yield bases are absent on this page; C70 trace regioisomers are unidentified."],
    6: ["Table 1 ratio/scale and General Procedure charges differ slightly and are retained as conflicts; pages beyond 6 were not claimed as covered."],
}
coverage = {"document_id": DOC, "semantic_accuracy_verified": False, "pages": [
    {"page": i, "state": "attempted", "regions": regions[i], "record_ids": page_records[i], "unresolved": unresolved[i]}
    for i in range(1, 7)
]}
dump(RUN / "coverage/pages_0001_0006.json", coverage)

raw_inventory = {
    "document_id": DOC, "page_range": [1, 6], "source_review_status": "assistant visual inventory retained as raw model output; formal records remain needs_review",
    "pages": [{"page": i, "visual_inspection_attempted": True,
               "image_path": f"runs/current/assets/{DOC}/page{i}.jpg",
               "image_evidence_id": I[i], "regions": regions[i],
               "candidate_record_ids": page_records[i], "unresolved": unresolved[i]}
              for i in range(1, 7)],
    "deduplication": [{"source_mentions":["Table 1 entry 4","Figure 2 product 3a","General Procedure for 3a"],
                       "retained_record_id":f"{DOC}:table1_entry4"}],
    "negative_inventory": "No standalone MS peak table, molecular-weight table, or fragment-relative-abundance block was converted into a reaction on pages 1-6."
}
dump(RUN / "raw_model/page_extraction.json", raw_inventory)

# This session's bounded visual inventory and coverage claim: exactly PDF pages 7-10.
session_page_records = {
    7: [f"{DOC}:scope_3{x}" for x in "cdefghijk"],
    8: [f"{DOC}:scope_3{x}" for x in "klmnop"] + [f"{DOC}:c70_cycloaddition_4a"],
    9: [f"{DOC}:c70_cycloaddition_4a", f"{DOC}:reduction_3a_to_5a", f"{DOC}:reduction_4a_to_6a"],
    10: [],
}
session_regions = {
    7: [
        {"region":"left column top continuation", "finding":"Continuation of compound 3c characterization only (UV-vis/ESI-HRMS); no new reaction arrow or procedure.", "evidence_ids":[E[7], I[7]]},
        {"region":"compound blocks 3d-3g", "finding":"Compound masses and the previously reported scope yields were linked to scope_3d through scope_3g. MW, Rf, IR, NMR, UV-vis, and HRMS remain characterization context, not separate synthesis experiments.", "evidence_ids":[E[7], I[7]]},
        {"region":"right-column compound blocks 3h-3k", "finding":"Compound masses and previously reported scope yields were linked to scope_3h through scope_3k; the 3k characterization continues onto page 8. No new conditions are stated.", "evidence_ids":[E[7], I[7]]},
    ],
    8: [
        {"region":"left-column top and compound blocks 3l-3o", "finding":"Continuation of 3k and characterization of 3l-3o. Product masses were linked to existing scope records; spectral data/MW were retained as source_context.", "evidence_ids":[E[8], I[8]]},
        {"region":"right-column compound 3p", "finding":"3p characterization reports 13.1 mg and 22%; the mass was linked to scope_3p without creating a duplicate experiment.", "evidence_ids":[E[8], I[8]]},
        {"region":"right-column Experimental Procedure for 4a", "finding":"Detailed C70 cycloaddition procedure (charges, H2 activation, 90 °C/40 min, chromatography, 57.4 mg, 46%, 3:1) was merged into the existing Scheme-2 c70_cycloaddition_4a record.", "evidence_ids":[E[8], I[8]]},
        {"region":"right-column Compound 4a block start", "finding":"Start of 4a mixture characterization; spectral data continue on page 9 and are not a new synthesis reaction.", "evidence_ids":[E[8], I[8]]},
    ],
    9: [
        {"region":"left-column top Compound 4a continuation", "finding":"4a mixture characterization continues; only the major alpha-4a carbon spectrum is identified. These spectral values remain source_context.", "evidence_ids":[E[9], I[9]]},
        {"region":"left-column General Procedure for 5a", "finding":"Detailed reduction of 3a to 5a was merged into reduction_3a_to_5a with all reported charges, ordered 0 °C addition and 5 h room-temperature hold, direct chromatography, 82.1 mg and 82% isolated yield.", "evidence_ids":[E[9], I[9]]},
        {"region":"left-column Compound 5a and 6a blocks", "finding":"5a/6a characterization was associated with the existing two reduction records. The 6a mixture mass (42.0 mg), 69% yield and 3:1 ratio were retained; no separate 6a procedure is printed here.", "evidence_ids":[E[9], I[9]]},
        {"region":"right column", "finding":"Acknowledgements and references [1]-[20]; no current-paper experiment or reaction arrow.", "evidence_ids":[E[9], I[9]]},
    ],
    10: [
        {"region":"both columns", "finding":"References [21]-[70] only. Bibliographic citations were not converted into current-paper experiments, and no reaction arrows, tables, footnotes, or procedures are present.", "evidence_ids":[E[10], I[10]]},
    ],
}
session_unresolved = {
    7: ["No experiment-specific procedure or yield measurement basis is supplied for scope members 3d-3k on this page; compact fullerene product structures remain not encoded."],
    8: ["For 4a, the total o-DCB volume/concentration denominator and identities/amounts of trace regioisomers are unreported; the alpha/beta product graphs remain not encoded."],
    9: ["For the C70 reduction to 6a, separate substrate/reagent charges, exact duration, workup/purification and yield measurement basis are not reported; absolute/relative stereochemical graphs remain unresolved."],
    10: ["No synthetic content occurs on this page; cited references were not inspected as source papers and no SI association was inferred."],
}
session_coverage = {"document_id": DOC, "semantic_accuracy_verified": False, "pages": [
    {"page": i, "state": "attempted", "regions": session_regions[i],
     "record_ids": session_page_records[i], "unresolved": session_unresolved[i]}
    for i in range(7, 11)
]}
dump(RUN / "coverage/pages_0007_0010.json", session_coverage)

session_raw_inventory = {
    "document_id": DOC,
    "page_range": [7, 10],
    "source_review_status": "assistant visual inventory retained as raw model output; formal records remain needs_review",
    "pages": [{"page": i, "visual_inspection_attempted": True,
               "image_path": f"runs/current/assets/{DOC}/page{i}.jpg",
               "image_evidence_id": I[i], "regions": session_regions[i],
               "candidate_record_ids": session_page_records[i], "unresolved": session_unresolved[i]}
              for i in range(7, 11)],
    "deduplication": [
        {"source_mentions":["Scheme 2", "page-8 Experimental Procedure for 4a", "page-8/9 Compound 4a characterization"],
         "retained_record_id":f"{DOC}:c70_cycloaddition_4a"},
        {"source_mentions":["Scheme 3 reduction of 3a", "page-9 General Procedure for 5a", "page-9 Compound 5a characterization"],
         "retained_record_id":f"{DOC}:reduction_3a_to_5a"},
        {"source_mentions":["Scheme 3 C70 reduction", "page-9 Compound 6a characterization"],
         "retained_record_id":f"{DOC}:reduction_4a_to_6a"},
    ],
    "negative_inventory": "Pages 7-9 spectral blocks and molecular weights were retained as characterization context; page 10 is references only. No standalone spectrum or bibliography entry was converted into a reaction."
}
dump(RUN / "raw_model/page_extraction_pages_0007_0010.json", session_raw_inventory)

stereo_checks = {}
for lab, smi in bis_smiles.items(): stereo_checks[lab] = st.stereo_inspect(smi)
dump(RUN / "raw_model/stereo_checks.json", stereo_checks)

print(json.dumps({"molecules": len(template["molecules"]), "reactions": len(R),
                  "specific_experiments": sum(x["record_type"] == "specific_experiment" for x in R),
                  "coverage_pages_this_session": [p["page"] for p in session_coverage["pages"]]}, ensure_ascii=False))
