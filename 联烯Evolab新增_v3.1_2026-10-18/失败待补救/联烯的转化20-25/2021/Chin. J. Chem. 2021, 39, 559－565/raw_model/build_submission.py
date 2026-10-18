#!/usr/bin/env python3
"""Build the evidence-linked v3.1 submission for PDF pages 1--6.

This is a machine-authored extraction artifact.  It deliberately preserves
unresolved primary structures instead of manufacturing stereochemical SMILES.
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
DATASET = json.loads((RUN / "dataset.json").read_text())
DOC = "doc_7cbfadf024ffc8d7"
EV = {p: f"{DOC}:page{p}:text" for p in range(1, 7)}


def f(raw=None, value=None, ev=(), method="original_extraction", derivation=None,
      status=None):
    if isinstance(ev, str):
        ev = [ev]
    return st.fact(raw, value, ev, method=method, derivation=derivation,
                   status=status)


def q(number, unit, qualifier=None, raw_text=None):
    out = {"number": number, "unit": unit}
    if qualifier:
        out["qualifier"] = qualifier
    if raw_text:
        out["raw_text"] = raw_text
    return out


def pct(number, qualifier=None, raw_text=None):
    return q(number, "%", qualifier, raw_text)


bundle = st.empty(DATASET["documents"])
molecules = {}


def molecule(mid, label, ev, *, chirality="unresolved", geometry="unresolved",
             representation="unresolved", assignment=None, axial=None,
             double_bonds=None, mixture=None, identity_note=None,
             smiles=None):
    evidence = [ev] if isinstance(ev, str) else list(ev)
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f(label, label, evidence)
    m["names"]["original"] = f(label, label, evidence)
    m["identity_note"] = f(identity_note, identity_note, evidence) if identity_note else f(ev=evidence)
    m["stereochemistry"]["chirality_status"] = f(chirality, chirality, evidence)
    m["stereochemistry"]["geometry_status"] = f(geometry, geometry, evidence)
    m["stereochemistry"]["representation_status"] = f(representation, representation, evidence)
    if assignment is not None:
        m["stereochemistry"]["assignment"] = f(assignment, assignment, evidence)
    if axial is not None:
        m["stereochemistry"]["axial"] = f(axial, axial, evidence)
    if double_bonds is not None:
        m["stereochemistry"]["double_bonds"] = f(double_bonds, double_bonds, evidence)
    if mixture is not None:
        m["stereochemistry"]["mixture"] = f(mixture, mixture, evidence)
    if smiles:
        diag = st.stereo_inspect(smiles)
        deriv = "Connectivity candidate derived from the explicit source chemical name; graph parsing is not source-identity verification."
        m["representation_type"] = "smiles"
        m["smiles"] = f(label, smiles, evidence, "structure_derivation", deriv)
        m["stereochemistry"]["chirality_status"] = f("graph has no stereogenic element", "achiral", evidence, "structure_derivation", deriv)
        m["stereochemistry"]["geometry_status"] = f("not applicable", "not_applicable", evidence, "structure_derivation", deriv)
        m["stereochemistry"]["representation_status"] = f("encoded", "encoded", evidence, "structure_derivation", deriv)
        m["stereochemistry"]["isomeric_smiles"] = f(smiles, diag["canonical_isomeric_smiles"], evidence, "structure_derivation", deriv)
        m["stereochemistry"]["identity_key"] = f(None, diag["stereo_graph_key"], evidence, "structure_derivation", deriv)
    molecules[mid] = m
    return m


def template_molecule(mid, label, ev, variables, fixed, *, chirality="unresolved",
                      geometry="partially_specified", assignment=None):
    m = molecule(mid, label, ev, chirality=chirality, geometry=geometry,
                 representation="not_encoded", assignment=assignment,
                 identity_note="Generic source structure; variable substituents are retained rather than expanded into a concrete molecule.")
    m["representation_type"] = "generic_template"
    m["generic_template"] = {
        "fixed_structure": f(fixed, {"description": fixed}, [ev]),
        "variables": [{
            "label": v,
            "attachment_points": f(f"source drawing attachment for {v}", {"source_label": v}, [ev]),
            "definition": f(f"variable {v} as drawn", f"variable {v} as drawn", [ev]),
            "allowed_values": f(ev=[ev]),
            "constraints": f(ev=[ev]),
        } for v in variables],
        "instance_mapping": [],
        "reason": f("R-group source drawing", "Concrete substituent identity is not fixed by this generic scheme.", [ev]),
    }
    return m


# Primary allene-dienes and products.  Source labels and stereochemical text are
# authoritative; no concrete primary-structure SMILES is guessed.
substituents = {
    "a": ("Me", "n-Bu", "H"), "b": ("Me", "Me", "H"),
    "c": ("Me", "Et", "H"), "d": ("Me", "n-Pr", "H"),
    "e": ("Me", "n-C5H11", "H"), "f": ("Me", "n-C6H13", "H"),
    "g": ("Me", "n-C7H15", "H"), "h": ("Me", "n-C8H17", "H"),
    "i": ("Me", "n-C9H19", "H"), "j": ("Me", "n-C10H21", "H"),
    "k": ("n-Pr", "n-Bu", "H"), "l": ("Me", "Cy", "H"),
    "m": ("Me", "(CH2)2OH", "H"), "n": ("H", "n-Bu", "Me"),
}
for suffix, groups in substituents.items():
    molecule(f"m_sub1{suffix}_rac", f"racemic (2E,4E)-1{suffix}", EV[3],
             chirality="racemic", geometry="specified", representation="not_encoded",
             axial=[{"element": "allene axis", "configuration": "racemic; individual Ra/Sa not separated"}],
             double_bonds=[{"source_assignment": "2E"}, {"source_assignment": "4E"}],
             mixture={"type": "racemate", "ratio": "not numerically reported"},
             identity_note=f"Table 3 instance: R1={groups[0]}, R2={groups[1]}, R3={groups[2]}; compact drawing not converted to a concrete SMILES.")
    for hand, cfg in (("S", "(3aS,7aS,E)"), ("R", "(3aR,7aR,E)")):
        molecule(f"m_prod3{suffix}_{hand}", f"{cfg}-3{suffix}", EV[3],
                 chirality="specified", geometry="specified", representation="not_encoded",
                 assignment={"ring_junction": "cis", "absolute_configuration": f"3a{hand},7a{hand}"},
                 double_bonds=[{"configuration": "E", "source": "compound label and drawing"}],
                 identity_note=f"Table 3 product for R1={groups[0]}, R2={groups[1]}, R3={groups[2]}; ring/alkene stereochemistry retained textually.")

# Specific stereoisomers and side products used in mechanistic experiments.
molecule("m_sub1a_S", "(Sa,2E,4E)-1a", EV[4], chirality="specified", geometry="specified",
         representation="not_encoded", axial=[{"element": "allene axis", "configuration": "Sa"}],
         double_bonds=[{"configuration": "2E"}, {"configuration": "4E"}])
molecule("m_sub1a_R", "(Ra,2E,4E)-1a", EV[4], chirality="specified", geometry="specified",
         representation="not_encoded", axial=[{"element": "allene axis", "configuration": "Ra"}],
         double_bonds=[{"configuration": "2E"}, {"configuration": "4E"}])
molecule("m_prod2a_SRS", "(3aS,5R,7aS,Z)-2a", EV[3], chirality="specified", geometry="specified",
         representation="not_encoded", assignment={"absolute_configuration": "3aS,5R,7aS", "ring_junction": "cis"},
         double_bonds=[{"configuration": "Z"}])
molecule("m_prod2a_rac", "racemic (3aR*,5S*,7aR*,Z)-2a", EV[4], chirality="racemic", geometry="specified",
         representation="not_encoded", assignment={"configuration_notation": "3aR*,5S*,7aR*"},
         double_bonds=[{"configuration": "Z"}], mixture={"type": "racemate"})
molecule("m_prod2a_opt", "normal [4+2]-cycloaddition product 2a (optimization sample)", EV[2],
         chirality="partially_specified", geometry="specified", representation="not_encoded",
         assignment={"source_table": "Table 1 drawing labelled (3aR,5S,7aR,Z)-2a; ee sign/order varies by ligand"},
         double_bonds=[{"configuration": "Z"}])
molecule("m_side1a_prime", "(2E,4E,9Z,11E)-1a'", EV[2], chirality="unresolved", geometry="specified",
         representation="not_encoded", double_bonds=[{"configuration": "2E"}, {"configuration": "4E"}, {"configuration": "9Z"}, {"configuration": "11E"}])
molecule("m_side1a_prime_D", "deuterated (2E,4E,9Z,11E)-1a''", EV[4], chirality="unresolved", geometry="specified",
         representation="not_encoded", double_bonds=[{"configuration": "2E"}, {"configuration": "4E"}, {"configuration": "9Z"}, {"configuration": "11E"}],
         identity_note="Scheme 2f deuterium-incorporation target; exact isotopic position is shown in the source drawing but no concrete SMILES is asserted.")
molecule("m_sub1o", "malonate-tethered (2E,4E)-1o", EV[4], chirality="unresolved", geometry="specified", representation="not_encoded",
         double_bonds=[{"configuration": "2E"}, {"configuration": "4E"}])
molecule("m_sub1p", "oxygen-tethered (2E,4E)-1p", EV[4], chirality="unresolved", geometry="specified", representation="not_encoded",
         double_bonds=[{"configuration": "2E"}, {"configuration": "4E"}])

# Derivatization products.
for label, note in {
    "4a": "double epoxide, n-Bu", "4c": "double epoxide, Et",
    "4f": "double epoxide, n-C6H13", "4h": "double epoxide, n-C8H17",
    "5a": "regioisomeric bis-dihydroxylation product",
    "6a": "regioisomeric bis-dihydroxylation product",
    "7a": "alpha-hydroxyl ketone", "8a": "ketal", "9a": "alpha,beta-unsaturated ketone",
}.items():
    molecule(f"m_prod{label}", label, EV[5], chirality="specified", geometry="partially_specified",
             representation="not_encoded", assignment={"source_drawing": "wedge/dash relative stereochemistry retained", "note": note},
             identity_note=f"Scheme 4 {note}; compact stereochemical drawing not converted to a concrete SMILES.")

# Generic structures for overview/mechanism records.
template_molecule("m_generic_sub1", "generic racemic (2E,4E)-1", EV[1], ["R1", "R2", "R3"],
                  "N-tosyl-tethered allene-1,3-diene shown in the comprehensive graphic",
                  chirality="racemic", assignment={"allene": "racemic"})
template_molecule("m_generic_prod3", "generic cis-(3aS,7aS,E)-3", EV[1], ["R1", "R2", "R3"],
                  "cis-aza-[4.3.0]bicyclic E-alkene product", chirality="specified", geometry="specified",
                  assignment={"ring_junction": "cis", "alkene": "E"})
template_molecule("m_generic_side1prime", "generic (2E,4E,9Z,11E)-1'", EV[1], ["R1", "R2", "R3", "R4"],
                  "conjugated tetra-ene side product", geometry="specified",
                  assignment={"double_bonds": "2E,4E,9Z,11E"})
template_molecule("m_generic_prod2", "generic (3aS,5R,7aS,Z)-2", EV[5], ["R1", "R2"],
                  "normal [4+2]-cycloaddition intermediate product", chirality="specified", geometry="specified",
                  assignment={"absolute_configuration": "3aS,5R,7aS", "alkene": "Z"})

# Agents.  Simple solvents receive graph candidates and are later checked by the
# standard stereo command; complex catalysts/ligands remain unresolved.
simple_smiles = {
    "m_EtOH": ("EtOH", "CCO"), "m_MeOH": ("MeOH", "CO"),
    "m_nPrOH": ("n-PrOH", "CCCO"), "m_toluene": ("toluene", "Cc1ccccc1"),
    "m_THF": ("THF", "C1CCOC1"), "m_DCM": ("DCM", "ClCCl"),
    "m_MeCN": ("MeCN", "CC#N"), "m_dioxane": ("1,4-dioxane", "O1CCOCC1"),
    "m_water": ("H2O", "O"), "m_acetone": ("acetone", "CC(=O)C"),
    "m_EtOAc": ("ethyl acetate", "CCOC(=O)C"),
    "m_iPrOH": ("i-PrOH", "CC(O)C"), "m_hexane": ("hexane", "CCCCCC"),
}
for mid, (label, smi) in simple_smiles.items():
    molecule(mid, label, EV[2] if label not in {"H2O", "acetone", "ethyl acetate", "i-PrOH", "hexane"} else EV[5], smiles=smi)
for mid, label, ev in [
    ("m_Rh_coe", "[Rh(coe)2Cl]2", EV[2]), ("m_Rh_eth", "[Rh(C2H4)2Cl]2", EV[3]),
    ("m_Rh_nbd", "[Rh(nbd)Cl]2", EV[3]), ("m_AgSbF6", "AgSbF6", EV[2]),
    ("m_AgOTf", "AgOTf", EV[3]), ("m_AgBF4", "AgBF4", EV[3]),
    ("m_AgNTf2", "AgNTf2", EV[3]), ("m_L1", "(R,R)-Norphos (L1)", EV[2]),
    ("m_L2", "(R,R)-Me-DUPHOS (L2)", EV[2]), ("m_L3", "ferrocene-based phosphine ligand L3", EV[2]),
    ("m_L4", "(S)-p-Tol-BINAP (L4)", EV[2]), ("m_L5", "menthol-derived ligand L5", EV[2]),
    ("m_L6", "Feringa ligand L6", EV[2]), ("m_L7_RR", "(R,R)-Ph-BPE (L7)", EV[2]),
    ("m_L7_SS", "(S,S)-Ph-BPE", EV[3]), ("m_EtOH_d6", "EtOH-d6", EV[4]),
    ("m_mCPBA", "mCPBA", EV[5]), ("m_NaHCO3", "NaHCO3", EV[5]),
    ("m_OsO4", "OsO4", EV[5]), ("m_NMO", "NMO", EV[5]),
    ("m_TsNH2", "TsNH2", EV[5]), ("m_DMP", "Dess-Martin periodinane (DMP)", EV[5]),
    ("m_TsOH", "TsOH", EV[5]), ("m_argon", "argon", EV[5]),
    ("m_silica", "silica gel", EV[5]), ("m_petether", "petroleum ether", EV[5]),
]:
    ch = "specified" if mid in {"m_L1", "m_L2", "m_L4", "m_L5", "m_L6", "m_L7_RR", "m_L7_SS"} else "unresolved"
    assign = {"source_name": label} if ch == "specified" else None
    molecule(mid, label, ev, chirality=ch, geometry="unresolved", representation="unresolved", assignment=assign)


def reaction(rid, ev, record_type="specific_experiment", categories=None, location=None,
             name=None, transformation=None):
    evidence = [ev] if isinstance(ev, str) else list(ev)
    r = st.new_reaction(f"{DOC}:{rid}", DOC, evidence)
    r["record_type"] = record_type
    r["experiment_categories"] = f(categories or [], categories or [], evidence)
    r["location"] = f(location, location, evidence)
    r["reaction_name"] = f(name, name, evidence) if name else f(ev=evidence)
    if transformation:
        r["transformation_type"] = f(transformation, transformation, evidence)
    r["review"]["unresolved"] = []
    r["completeness"] = "partial"
    bundle["reactions"].append(r)
    return r


def add_part(r, pid, mid, role, ev, *, stage="reaction", charges=None, ref=None, product_role=None):
    p = st.new_participant(pid, mid, role, [ev], stage)
    if product_role:
        p["product_role"] = f(product_role, product_role, [ev])
    for field, value in (charges or {}).items():
        raw = value.get("raw_text") if isinstance(value, dict) else str(value)
        if raw is None:
            raw = f"{value['number']} {value['unit']}"
        p["charges"][field] = f(raw, value, [ev])
    if ref:
        p["equivalence_reference"] = f(ref, ref, [ev])
    r["participants"].append(p)
    return p


def set_solvent(r, ev, composition, volume=None, ratio=None):
    r["solvent_conditions"]["composition"] = f(composition, composition, [ev])
    if volume is not None:
        r["solvent_conditions"]["total_volume"] = f(f"{volume['number']} {volume['unit']}", volume, [ev])
    if ratio is not None:
        r["solvent_conditions"]["mixture_ratio"] = f(ratio, ratio, [ev])


def add_step(r, sid, ev, *, op, temp=None, time=None, stage="reaction", pids=None,
             atmosphere=None):
    s = st.new_step(sid, [ev], stage)
    s["participant_ids"] = list(pids or [])
    s["operation"] = f(op, op, [ev])
    if temp is not None:
        raw = temp.get("raw_text") if isinstance(temp, dict) else str(temp)
        raw = raw or f"{temp['number']} {temp['unit']}"
        s["temperature"] = f(raw, temp, [ev])
    if time is not None:
        raw = time.get("raw_text") if isinstance(time, dict) else str(time)
        raw = raw or f"{time['number']} {time['unit']}"
        s["time"] = f(raw, time, [ev])
    if atmosphere is not None:
        s["atmosphere"] = f(atmosphere, atmosphere, [ev])
    r["steps"].append(s)
    return s


def add_outcome(r, oid, kind, targets, ev, value, *, basis=None, yield_type=None, scope=None,
                raw=None):
    o = st.new_outcome(oid, kind, targets, [ev])
    o["value"] = f(raw if raw is not None else str(value), value, [ev])
    if basis is not None:
        o["measurement_basis"] = f(basis, basis, [ev])
    if yield_type is not None:
        o["yield_type"] = f(yield_type, yield_type, [ev])
    if scope is not None:
        o["scope"] = f(scope, scope, [ev])
    r["outcomes"].append(o)
    return o


def set_class(r, ev, mode):
    """Record only source-legible net change classes; never claim full mapping."""
    changes = []
    if mode in {"cycloaddition", "cycloaddition_isomerization", "mixed_cyclo_side"}:
        for _ in range(2):
            changes.append(("formed", "C", "C", 0, 1))
    if mode in {"isomerization", "cycloaddition_isomerization", "mixed_cyclo_side"}:
        changes += [("hydrogen_gain", "C", "H", 0, 1), ("hydrogen_loss", "C", "H", 1, 0)]
    if mode == "epoxidation":
        changes += [("formed", "C", "O", 0, 1)] * 4
    if mode == "dihydroxylation":
        changes += [("formed", "C", "O", 0, 1)] * 4
    if mode == "oxidation":
        changes += [("order_increased", "C", "O", 1, 2), ("hydrogen_loss", "C", "H", 1, 0), ("hydrogen_loss", "O", "H", 1, 0)]
    if mode == "ketal_dehydration":
        changes += [("formed", "C", "O", 0, 1)] * 2 + [("order_increased", "C", "C", 1, 2)]
    if mode in {"none", "unresolved"}:
        r["classification"]["coverage"] = f(mode, "not_applicable" if mode == "none" else "unresolved", [ev], "inference",
                                                 "No observed transformation for this control." if mode == "none" else "Source-level atom mapping was not attempted.")
        r["classification"]["labels"] = f(mode, [] if mode == "none" else None, [ev], "inference",
                                               "No formed-bond labels apply." if mode == "none" else "Classification unresolved.")
        return
    facts = []
    labels = {}
    for i, (kind, a, b, before, after) in enumerate(changes, 1):
        change = {"change_id": f"b{i}", "atoms": [{"atom_map": None, "element": a}, {"atom_map": None, "element": b}],
                  "before_order": before, "after_order": after, "change_type": kind}
        if "hydrogen" in kind:
            change["hydrogen_count"] = 1
        facts.append(f(None, change, [ev], "structure_derivation",
                       "Partial net change inferred from the source reaction drawing; atom positions are not fully mapped."))
        if kind in {"formed", "hydrogen_gain", "order_increased"}:
            lab = st.bond_label(a, b, after)
            labels.setdefault(lab["code"], {**lab, "supporting_change_ids": []})["supporting_change_ids"].append(f"b{i}")
    r["classification"]["bond_changes"] = facts
    r["classification"]["labels"] = f(None, list(labels.values()), [ev], "structure_derivation",
                                          "Labels include all and only the explicitly listed partial formed/promoted bonds.")
    r["classification"]["coverage"] = f("partial source-structure derivation", "partial", [ev], "inference",
                                            "No complete atom mapping; omitted partners, exact atom identities and additional bond-order changes remain unresolved.")
    r["classification"]["limitations"] = [
        "Partial source-structure classification only; no complete atom mapping is claimed.",
        "Hydrogen relocation is represented only when the net allene/alkene isomerization is part of the observed transformation.",
        "Unlisted coproducts and exact bond-order decreases are unresolved.",
    ]


def add_standard_conditions(r, ev, substrate_mid, *, catalyst="m_Rh_eth", ligand="m_L7_RR",
                            silver="m_AgSbF6", solvent="m_EtOH", scale=None,
                            catalyst_pct=2.5, ligand_pct=5, silver_pct=5, volume=None,
                            temp=45, time=3, include_ligand=True):
    add_part(r, "substrate", substrate_mid, "reactant", ev,
             charges={"amount": q(scale, "mmol")} if scale is not None else {})
    if scale is not None:
        scale_value = {"reference_molecule_id": substrate_mid, "amount": q(scale, "mmol")}
        r["reaction_scale"] = f(f"{scale} mmol of {substrate_mid}", scale_value, [ev])
    add_part(r, "catalyst", catalyst, "catalyst", ev, charges={"mol_percent": q(catalyst_pct, "mol%")}, ref=substrate_mid)
    if include_ligand:
        add_part(r, "ligand", ligand, "ligand", ev, charges={"mol_percent": q(ligand_pct, "mol%")}, ref=substrate_mid)
    add_part(r, "silver", silver, "additive", ev, charges={"mol_percent": q(silver_pct, "mol%")}, ref=substrate_mid)
    add_part(r, "solvent", solvent, "solvent", ev, charges={"volume": volume} if volume else {})
    set_solvent(r, ev, molecules[solvent]["names"]["original"]["normalized_value"], volume)
    add_step(r, "reaction", ev, op="standard conditions", temp=q(temp, "°C"), time=q(time, "h"), pids=[p["participant_id"] for p in r["participants"]])


# Overview generic reaction (page 1, repeated conceptually in Scheme 1b-ii).
r = reaction("overview_generic", [EV[1], EV[2]], "generic_scheme", [], "page 1 comprehensive graphic; page 2 Scheme 1b-ii",
             "Rh-catalyzed kinetic resolution-based enantioselective [4+2]-cycloaddition-isomerization",
             "intramolecular [4+2]-cycloaddition followed by C=C isomerization")
for pid, mid, role, pr in [("generic_sub", "m_generic_sub1", "reactant", None),
                           ("generic_bicycle", "m_generic_prod3", "product", "main product"),
                           ("generic_tetraene", "m_generic_side1prime", "product", "parallel isomerization product")]:
    add_part(r, pid, mid, role, EV[1], product_role=pr)
add_part(r, "catalyst", "m_Rh_eth", "catalyst", EV[1], charges={"mol_percent": q(2.5, "mol%")}, ref="m_generic_sub1")
add_part(r, "ligand", "m_L7_RR", "ligand", EV[1], charges={"mol_percent": q(5, "mol%")}, ref="m_generic_sub1")
add_part(r, "silver", "m_AgSbF6", "additive", EV[1], charges={"mol_percent": q(5, "mol%")}, ref="m_generic_sub1")
add_part(r, "solvent", "m_EtOH", "solvent", EV[1])
set_solvent(r, EV[1], "EtOH")
add_step(r, "generic_reaction", EV[1], op="generic reaction arrow", temp=q(45, "°C"), pids=[p["participant_id"] for p in r["participants"]])
set_class(r, EV[1], "cycloaddition_isomerization")
r["review"]["unresolved"] = ["Generic scheme has variable R groups and no experiment-specific scale, time or yields."]

# Scheme 1b-i is explicitly prior work, not a current-paper experiment.
r = reaction("scheme1_prior_chirality_transfer", EV[2], "literature_example", [], "page 2 Scheme 1b-i",
             "Chirality transfer approach", "prior-literature Rh-catalyzed intramolecular [4+2] cycloaddition")
r["notes"].append(f("cat. Rh; ref. 7", "Prior work cited as ref. 7; not counted as a current-paper experiment.", [EV[2]]))
set_class(r, EV[2], "unresolved")
r["review"]["unresolved"] = ["Prior-work generic structures and conditions are not expanded from the cited paper."]

# Table 1: ligand screen.
table1 = [
    (1, "m_L1", 6, 0, 59, 42, 0, None), (2, "m_L2", 6, 0, 64, 56, 0, None),
    (3, "m_L3", 19, 90, "trace", None, 0, None), (4, "m_L4", 9, 58, "trace", None, 0, None),
    (5, "m_L5", 6, 0, 56, 1, 0, None), (6, "m_L6", 4, 0, 79, 3, 0, None),
    (7, "m_L7_RR", 3, 0, 0, None, 29, -95),
]
for entry, lig, hours, recovery, y2, ee2, y3, ee3 in table1:
    r = reaction(f"table1_entry{entry}", EV[2], categories=["optimization"], location=f"page 2 Table 1 entry {entry}",
                 name="Ligand screening", transformation="intramolecular [4+2]-cycloaddition with possible C=C isomerization")
    scale = 0.2 if entry == 7 else 0.1
    vol = q(2.0 if entry == 7 else 1.0, "mL")
    add_part(r, "substrate", "m_sub1a_rac", "reactant", EV[2], charges={"amount": q(scale, "mmol")})
    add_part(r, "catalyst", "m_Rh_coe", "catalyst", EV[2], charges={"mol_percent": q(3, "mol%")}, ref="m_sub1a_rac")
    add_part(r, "ligand", lig, "ligand", EV[2], charges={"mol_percent": q(12 if entry == 5 else 6, "mol%")}, ref="m_sub1a_rac")
    add_part(r, "silver", "m_AgSbF6", "additive", EV[2], charges={"mol_percent": q(6, "mol%")}, ref="m_sub1a_rac")
    add_part(r, "solvent", "m_EtOH", "solvent", EV[2], charges={"volume": vol})
    r["reaction_scale"] = f(f"{scale} mmol of (2E,4E)-1a", {"reference_molecule_id": "m_sub1a_rac", "amount": q(scale, "mmol")}, [EV[2]])
    add_part(r, "prod2", "m_prod2a_opt", "product", EV[2], product_role="normal [4+2] product")
    add_part(r, "prod3", "m_prod3a_S", "product", EV[2], product_role="cycloaddition-isomerization product")
    if entry == 7:
        add_part(r, "side", "m_side1a_prime", "product", EV[2], product_role="tetra-ene side product")
    set_solvent(r, EV[2], "EtOH", vol)
    add_step(r, "reaction", EV[2], op=f"Table 1 entry {entry}", temp=q(45 if entry == 7 else 60, "°C"), time=q(hours, "h"), pids=[p["participant_id"] for p in r["participants"]])
    add_outcome(r, "recovery_1a", "recovery", ["m_sub1a_rac"], EV[2], pct(recovery), basis="1H NMR using mesitylene as internal standard", raw=f"{recovery}%")
    if y2 == "trace":
        add_outcome(r, "trace_2a", "trace", ["m_prod2a_opt"], EV[2], "trace", basis="1H NMR using mesitylene as internal standard")
    else:
        add_outcome(r, "yield_2a", "yield", ["m_prod2a_opt"], EV[2], pct(y2), basis="1H NMR using mesitylene as internal standard", raw=f"{y2}%")
    if ee2 is not None:
        add_outcome(r, "ee_2a", "selectivity", ["m_prod2a_opt"], EV[2], pct(ee2), basis="HPLC analysis", scope={"metric": "ee", "reference_order": "as reported in Table 1"}, raw=f"{ee2}% ee")
    add_outcome(r, "yield_3a", "yield", ["m_prod3a_S"], EV[2], pct(y3), basis="1H NMR using mesitylene as internal standard", raw=f"{y3}%")
    if ee3 is not None:
        add_outcome(r, "ee_3a", "selectivity", ["m_prod3a_S"], EV[2], pct(abs(ee3)), basis="HPLC analysis", scope={"metric": "signed ee", "reported_signed_value": f"{ee3}%", "raw_sign_retained": True}, raw=f"{ee3}% ee")
    if entry == 7:
        add_outcome(r, "yield_1a_prime", "yield", ["m_side1a_prime"], EV[2], pct(47), basis="1H NMR using mesitylene as internal standard", raw="47% NMR yield")
    set_class(r, EV[2], "cycloaddition_isomerization" if entry == 7 else "cycloaddition")
    r["review"]["unresolved"] = ["Primary structures remain unencoded; Table 1 ee sign/reference is retained without reinterpretation."]

# Table 2: catalyst/salt/solvent optimization.
table2 = [
    (1, "m_Rh_coe", "m_AgSbF6", "m_EtOH", 0, 31, ">99", 45),
    (2, "m_Rh_coe", "m_AgOTf", "m_EtOH", 85, 0, None, 0),
    (3, "m_Rh_coe", "m_AgBF4", "m_EtOH", 0, 24, ">99", 44),
    (4, "m_Rh_coe", "m_AgNTf2", "m_EtOH", 0, 31, 98, 46),
    (5, "m_Rh_coe", "m_AgSbF6", "m_nPrOH", 0, 18, 98, 41),
    (6, "m_Rh_coe", "m_AgSbF6", "m_MeOH", 74, 0, None, 8),
    (7, "m_Rh_coe", "m_AgSbF6", "m_toluene", 70, 0, None, 0),
    (8, "m_Rh_nbd", "m_AgSbF6", "m_EtOH", 57, 0, None, 7),
    (9, "m_Rh_eth", "m_AgSbF6", "m_EtOH", 0, 36, ">99", 52),
]
for entry, cat, silver, solvent, recovery, y3, ee, yside in table2:
    r = reaction(f"table2_entry{entry}", EV[3], categories=["optimization"], location=f"page 3 Table 2 entry {entry}",
                 name="Optimization of reaction conditions", transformation="intramolecular [4+2]-cycloaddition-isomerization / parallel allene isomerization")
    scale = 0.1 if entry in {2, 3, 4} else 0.2
    volume = q(1.0 if entry in {2, 3, 4} else 2.0, "mL")
    add_standard_conditions(r, EV[3], "m_sub1a_rac", catalyst=cat, silver=silver, solvent=solvent, scale=scale, volume=volume)
    add_part(r, "prod3", "m_prod3a_S", "product", EV[3], product_role="desired product")
    add_part(r, "side", "m_side1a_prime", "product", EV[3], product_role="tetra-ene side product")
    add_outcome(r, "recovery_1a", "recovery", ["m_sub1a_rac"], EV[3], pct(recovery), basis="1H NMR using mesitylene as internal standard", raw=f"{recovery}%")
    add_outcome(r, "yield_3a", "yield", ["m_prod3a_S"], EV[3], pct(y3), basis="1H NMR using mesitylene as internal standard", raw=f"{y3}%")
    if ee is not None:
        val = pct(99, "gt", ">99%") if ee == ">99" else pct(ee)
        add_outcome(r, "ee_3a", "selectivity", ["m_prod3a_S"], EV[3], val, basis="HPLC analysis", scope={"metric": "ee"}, raw=f"{ee}% ee")
    add_outcome(r, "yield_side", "yield", ["m_side1a_prime"], EV[3], pct(yside), basis="1H NMR using mesitylene as internal standard", raw=f"{yside}%")
    set_class(r, EV[3], "cycloaddition_isomerization" if y3 else ("isomerization" if yside else "none"))
    r["review"]["unresolved"] = ["Primary reactant/product structures remain unencoded."]

# Table 3: each ligand column is a distinct experiment.  The 1a left-column
# record is also the detailed 1.0 mmol experimental procedure on pages 5--6.
table3 = [
    (1, "b", 34, "96", 41, "98"), (2, "c", 36, ">99", 36, ">99"),
    (3, "d", 30, "99", 26, "98"), (4, "a", 37, ">99", 34, "96"),
    (5, "e", 34, ">99", 33, ">99"), (6, "f", 35, ">99", 35, "96"),
    (7, "g", 34, ">99", 35, ">99"), (8, "h", 31, "99", 33, "98"),
    (9, "i", 36, ">99", 32, "98"), (10, "j", 33, "98", 31, "97"),
    (11, "k", 34, "98", 33, ">99"), (12, "l", 26, "97", 24, "93"),
    (13, "m", 34, ">99", 33, ">99"), (14, "n", 40, "92", 43, "92"),
]
for entry, suffix, yS, eeS, yR, eeR in table3:
    for hand, lig, scale, volume, yld, ee in (("S", "m_L7_RR", 1.0, 5.0, yS, eeS), ("R", "m_L7_SS", 0.5, 2.5, yR, eeR)):
        rid = f"table3_entry{entry}_{hand}"
        evs = [EV[3], EV[5], EV[6]] if entry == 4 and hand == "S" else [EV[3]]
        r = reaction(rid, evs, categories=["substrate_scope"] + (["scale_up"] if entry == 4 and hand == "S" else []),
                     location=f"page 3 Table 3 entry {entry}, {hand}-configured product column" + ("; pages 5-6 experimental procedure" if entry == 4 and hand == "S" else ""),
                     name="Rh-catalyzed intramolecular [4+2]-cycloaddition-C=C bond isomerization",
                     transformation="kinetic-resolution-based cycloaddition-isomerization")
        hours = 24 if entry == 11 or (entry == 4 and hand == "R") else 12 if entry in {12, 13} else 3
        temperature = 60 if entry == 6 else 45
        actual_scale = 0.2 if entry == 13 else scale
        if entry == 4 and hand == "S":
            # Detailed procedure supersedes the single generic reaction step.
            add_part(r, "substrate", "m_sub1a_rac", "reactant", EV[5], charges={"mass": q(359.6, "mg"), "amount": q(1.0, "mmol")})
            add_part(r, "catalyst", "m_Rh_eth", "catalyst", EV[5], charges={"mass": q(9.8, "mg"), "amount": q(0.025, "mmol"), "mol_percent": q(2.5, "mol%")}, ref="m_sub1a_rac")
            add_part(r, "ligand", lig, "ligand", EV[5], charges={"mass": q(25.4, "mg"), "amount": q(0.05, "mmol"), "mol_percent": q(5, "mol%")}, ref="m_sub1a_rac")
            add_part(r, "silver", "m_AgSbF6", "additive", EV[5], charges={"mass": q(17.1, "mg"), "amount": q(0.05, "mmol"), "mol_percent": q(5, "mol%")}, ref="m_sub1a_rac")
            add_part(r, "solvent", "m_EtOH", "solvent", EV[5], charges={"volume": q(5.0, "mL")})
            add_part(r, "workup_eluent", "m_EtOAc", "workup_agent", EV[5], stage="workup", charges={"volume": q(60, "mL")})
            add_part(r, "silica", "m_silica", "purification_agent", EV[5], stage="purification")
            add_part(r, "petether", "m_petether", "purification_agent", EV[5], stage="purification")
            r["reaction_scale"] = f("1.0 mmol of (2E,4E)-1a", {"reference_molecule_id": "m_sub1a_rac", "amount": q(1.0, "mmol")}, [EV[5]])
            set_solvent(r, EV[5], "EtOH", q(5.0, "mL"))
            add_step(r, "catalyst_preformation", EV[5], op="add catalyst, ligand and AgSbF6 in glove box; replace air with argon three times; add 2.5 mL EtOH and stir", temp="r.t.", time=q(10, "min"), atmosphere="argon")
            add_step(r, "substrate_addition_and_reaction", EV[5], op="add (2E,4E)-1a then 2.5 mL EtOH; place in preheated oil bath and stir, monitored by TLC", temp=q(45, "°C"), time=q(3, "h"), atmosphere="argon")
            add_step(r, "workup", EV[5], op="cool to rt; filter through a short silica pad eluted with ethyl acetate (60 mL); concentrate", stage="workup")
            add_step(r, "purification", EV[5], op="silica-gel chromatography; petroleum ether/ethyl acetate 40/1 (500 mL) to 20/1 (300 mL)", stage="purification")
            r["procedure"]["original_text"] = f("Asymmetric synthesis of (3aS,7aS,E)-3a ...", "Full printed procedure spans page 5 and characterization continues on page 6.", [EV[5], EV[6]])
            r["procedure"]["addition_order"] = f("catalyst/ligand/Ag salt -> EtOH -> substrate -> EtOH", ["catalyst", "ligand", "silver", "EtOH", "substrate", "EtOH"], [EV[5]])
            r["procedure"]["pretreatment"] = f("oven-dried Schlenk tube; glove box; three vacuum/argon cycles", "oven-dried Schlenk tube; glove box; three vacuum/argon cycles", [EV[5]])
            r["notes"].append(f("[Rh(C2H4)2Cl]2 (9.8 mg, 0.025 mmol, 99%)", "Catalyst purity is reported as 99%.", [EV[5]]))
        else:
            add_standard_conditions(r, EV[3], f"m_sub1{suffix}_rac", ligand=lig, scale=actual_scale, volume=q(volume if entry != 13 else 5.0, "mL"), temp=temperature, time=hours)
        add_part(r, "product", f"m_prod3{suffix}_{hand}", "product", EV[3], product_role="desired enantiomer")
        add_outcome(r, "isolated_yield", "yield", [f"m_prod3{suffix}_{hand}"], EV[3], pct(yld), basis="isolated yield", yield_type="isolated", raw=f"{yld}%")
        eeval = pct(99, "gt", ">99%") if ee == ">99" else pct(int(ee))
        add_outcome(r, "ee", "selectivity", [f"m_prod3{suffix}_{hand}"], EV[3], eeval, basis="HPLC analysis", scope={"metric": "ee", "product_configuration": hand}, raw=f"{ee}% ee")
        if entry == 4 and hand == "S":
            add_outcome(r, "product_mass", "product_amount", ["m_prod3a_S"], EV[5], q(133.1, "mg"), basis="weighed isolated product", raw="133.1 mg")
        set_class(r, EV[3], "cycloaddition_isomerization")
        r["review"]["unresolved"] = ["Primary substrate/product structures are source-identified but not encoded as concrete stereochemical SMILES."]


def mechanistic_record(rid, panel, substrate, hours=3, ligand="m_L7_RR", include_ligand=True,
                       solvent="m_EtOH", scale=None, temp=45, volume=None):
    r = reaction(rid, EV[4], categories=["mechanistic_study"], location=f"page 4 Scheme 2{panel}",
                 name="Mechanistic study", transformation="control/isomerization/cycloaddition experiment")
    add_standard_conditions(r, EV[4], substrate, ligand=ligand, include_ligand=include_ligand,
                            solvent=solvent, scale=scale, temp=temp, time=hours, volume=volume)
    return r


# Scheme 2a: ligand dependence of 2a -> 3a.
for tag, lig, include, y, ee, rec in [("no_ligand", "m_L7_RR", False, None, None, 100),
                                      ("RR", "m_L7_RR", True, 93, 95, 5),
                                      ("SS", "m_L7_SS", True, 96, 94, 4)]:
    r = mechanistic_record(f"scheme2a_{tag}", "a", "m_prod2a_SRS", hours=24, ligand=lig, include_ligand=include)
    if y is not None:
        add_part(r, "product", "m_prod3a_S", "product", EV[4], product_role="C=C-isomerized product")
        add_outcome(r, "yield_3a", "yield", ["m_prod3a_S"], EV[4], pct(y), basis="1H NMR (labelled NMR yield)", raw=f"{y}% NMR yield")
        add_outcome(r, "ee_3a", "selectivity", ["m_prod3a_S"], EV[4], pct(ee), scope={"metric": "ee"}, raw=f"{ee}% ee")
        set_class(r, EV[4], "isomerization")
    else:
        set_class(r, EV[4], "none")
        add_outcome(r, "no_reaction", "no_reaction", [], EV[4], True, raw="100% recovery; no reaction")
    add_outcome(r, "recovery_2a", "recovery", ["m_prod2a_SRS"], EV[4], pct(rec), raw=f"{rec}% recovery")
    r["review"]["unresolved"] = ["Recovery measurement basis is not stated in panel 2a."] + ([] if y is None else ["ee measurement method is not stated in panel 2a."])

# Scheme 2b.
r = mechanistic_record("scheme2b_racemic_2a", "b", "m_prod2a_rac", hours=24)
add_part(r, "product", "m_prod3a_S", "product", EV[4], product_role="C=C-isomerized product")
add_outcome(r, "yield", "yield", ["m_prod3a_S"], EV[4], pct(88), raw="88%")
add_outcome(r, "ee", "selectivity", ["m_prod3a_S"], EV[4], pct(8), scope={"metric": "ee"}, raw="8% ee")
set_class(r, EV[4], "isomerization")
r["review"]["unresolved"] = ["Yield and ee measurement bases are not stated for Scheme 2b."]

# Scheme 2c solvent controls.
for entry, solvent, rec, y2, ee2, y3, yside in [(1, "m_THF", 43, 13, 38, 0, 20), (2, "m_DCM", 32, 18, 66, 0, 7),
                                                 (3, "m_MeCN", 77, 0, None, 0, 0), (4, "m_dioxane", 84, 0, None, 0, 0)]:
    r = mechanistic_record(f"scheme2c_entry{entry}", "c", "m_sub1a_rac", solvent=solvent, scale=0.1,
                           volume=q(1.0, "mL"))
    add_part(r, "prod2", "m_prod2a_SRS", "product", EV[4], product_role="normal cycloaddition product")
    add_part(r, "prod3", "m_prod3a_S", "product", EV[4], product_role="cycloaddition-isomerization product")
    add_part(r, "side", "m_side1a_prime", "product", EV[4], product_role="tetra-ene product")
    add_outcome(r, "recovery", "recovery", ["m_sub1a_rac"], EV[4], pct(rec), basis="1H NMR using mesitylene as internal standard", raw=f"{rec}%")
    add_outcome(r, "yield_2a", "yield", ["m_prod2a_SRS"], EV[4], pct(y2), basis="1H NMR using mesitylene as internal standard", raw=f"{y2}%")
    if ee2 is not None:
        add_outcome(r, "ee_2a", "selectivity", ["m_prod2a_SRS"], EV[4], pct(ee2), basis="HPLC analysis", scope={"metric": "ee"}, raw=f"{ee2}% ee")
    add_outcome(r, "yield_3a", "yield", ["m_prod3a_S"], EV[4], pct(y3), basis="1H NMR using mesitylene as internal standard", raw=f"{y3}%")
    add_outcome(r, "yield_side", "yield", ["m_side1a_prime"], EV[4], pct(yside), basis="1H NMR using mesitylene as internal standard", raw=f"{yside}%")
    set_class(r, EV[4], "mixed_cyclo_side" if y2 or yside else "none")
    r["review"]["unresolved"] = ["Primary structures remain unencoded."]

# Scheme 2d time course.
for entry, hours, rec, rec_ee, y2, ee2, y3, ee3, yside in [
    (1, 0.5, 47, 31, 11, ">99", 12, ">99", 29), (2, 1.0, 17, 63, 10, ">99", 23, ">99", 44),
    (3, 1.5, 9, 71, 9, ">99", 29, ">99", 46), (4, 3.0, 0, None, 0, None, 36, ">99", 52)]:
    r = mechanistic_record(f"scheme2d_time{entry}", "d", "m_sub1a_rac", hours=hours, scale=0.2,
                           volume=q(1.0, "mL"))
    for pid, mid, pr in [("prod2", "m_prod2a_SRS", "normal cycloaddition product"), ("prod3", "m_prod3a_S", "final product"), ("side", "m_side1a_prime", "tetra-ene product")]:
        add_part(r, pid, mid, "product", EV[4], product_role=pr)
    add_outcome(r, "recovery", "recovery", ["m_sub1a_rac"], EV[4], pct(rec), basis="1H NMR using mesitylene as internal standard", raw=f"{rec}%")
    if rec_ee is not None:
        add_outcome(r, "recovered_ee", "selectivity", ["m_sub1a_rac"], EV[4], pct(rec_ee), basis="HPLC analysis", scope={"metric": "ee", "major": "Sa recovered starting material"}, raw=f"{rec_ee}% ee")
    add_outcome(r, "yield_2a", "yield", ["m_prod2a_SRS"], EV[4], pct(y2), basis="1H NMR using mesitylene as internal standard", raw=f"{y2}%")
    if ee2:
        add_outcome(r, "ee_2a", "selectivity", ["m_prod2a_SRS"], EV[4], pct(99, "gt", ">99%"), basis="HPLC analysis", scope={"metric": "ee"}, raw=ee2 + "% ee")
    add_outcome(r, "yield_3a", "yield", ["m_prod3a_S"], EV[4], pct(y3), basis="1H NMR using mesitylene as internal standard", raw=f"{y3}%")
    add_outcome(r, "ee_3a", "selectivity", ["m_prod3a_S"], EV[4], pct(99, "gt", ">99%"), basis="HPLC analysis", scope={"metric": "ee"}, raw=ee3 + "% ee")
    add_outcome(r, "yield_side", "yield", ["m_side1a_prime"], EV[4], pct(yside), basis="1H NMR using mesitylene as internal standard", raw=f"{yside}%")
    set_class(r, EV[4], "cycloaddition_isomerization")
    r["review"]["unresolved"] = ["Primary structures remain unencoded."]

# Scheme 2e enantiopure starting materials.
for hand, sub, y3, ee3, yside, basis3, basisside in [
    ("S", "m_sub1a_S", 61, ">99", 30, None, None),
    ("R", "m_sub1a_R", 0, None, 84, None, "isolated yield")]:
    r = mechanistic_record(f"scheme2e_{hand}", "e", sub, hours=3)
    add_part(r, "prod3", "m_prod3a_S", "product", EV[4], product_role="cycloaddition-isomerization product")
    add_part(r, "side", "m_side1a_prime", "product", EV[4], product_role="tetra-ene product")
    add_outcome(r, "yield_3a", "yield", ["m_prod3a_S"], EV[4], pct(y3), basis=basis3, raw=f"{y3}%")
    if ee3:
        add_outcome(r, "ee_3a", "selectivity", ["m_prod3a_S"], EV[4], pct(99, "gt", ">99%"), scope={"metric": "ee"}, raw=ee3 + "% ee")
    add_outcome(r, "yield_side", "yield", ["m_side1a_prime"], EV[4], pct(yside), basis=basisside, yield_type="isolated" if basisside else None, raw=f"{yside}%")
    set_class(r, EV[4], "cycloaddition_isomerization" if hand == "S" else "isomerization")
    r["review"]["unresolved"] = ["Measurement basis is unreported for Scheme 2e outcomes except the explicitly isolated 84% tetra-ene yield."]

# Scheme 2f EtOH-d6 experiments.
for hand, sub, y3, ee3, yside, yD in [("R", "m_sub1a_R", 0, None, 79, 0), ("S", "m_sub1a_S", 54, ">99", 26, None)]:
    r = mechanistic_record(f"scheme2f_{hand}_EtOHd6", "f", sub, hours=3, solvent="m_EtOH_d6")
    add_part(r, "prod3", "m_prod3a_S", "product", EV[4], product_role="cycloaddition-isomerization product")
    add_part(r, "side", "m_side1a_prime", "product", EV[4], product_role="non-deuterated tetra-ene")
    if yD is not None:
        add_part(r, "sideD", "m_side1a_prime_D", "product", EV[4], product_role="deuterated tetra-ene target")
    add_outcome(r, "yield_3a", "yield", ["m_prod3a_S"], EV[4], pct(y3), raw=f"{y3}%")
    if ee3:
        add_outcome(r, "ee_3a", "selectivity", ["m_prod3a_S"], EV[4], pct(99, "gt", ">99%"), scope={"metric": "ee"}, raw=ee3 + "% ee")
    add_outcome(r, "yield_side", "yield", ["m_side1a_prime"], EV[4], pct(yside), raw=f"{yside}%")
    if yD is not None:
        add_outcome(r, "yield_sideD", "yield", ["m_side1a_prime_D"], EV[4], pct(yD), raw=f"{yD}%")
    set_class(r, EV[4], "cycloaddition_isomerization" if hand == "S" else "isomerization")
    r["review"]["unresolved"] = ["Outcome measurement bases are not stated in Scheme 2f; isotopic structure remains textually represented."]

# Scheme 2g negative tether controls.
for suffix, sub, temp, hours in [("1p", "m_sub1p", 45, 12), ("1o", "m_sub1o", 60, 24)]:
    r = mechanistic_record(f"scheme2g_{suffix}", "g", sub, hours=hours, temp=temp)
    add_outcome(r, "no_expected_product", "no_reaction", [], EV[4], True, raw="expected product not obtained")
    add_outcome(r, "recovery", "recovery", [sub], EV[4], pct(0), raw="no recovery")
    set_class(r, EV[4], "none")
    r["review"]["unresolved"] = ["No product identity is fabricated for this failed control; recovery measurement basis is unreported."]

# Scheme 3 proposed mechanism (not an experimental record).
r = reaction("scheme3_proposed_mechanism", EV[5], "proposed_mechanism", [], "page 5 Scheme 3",
             "Proposed mechanism", "author-proposed Rh-catalyzed pathways a and b")
for pid, mid, role in [("sub", "m_generic_sub1", "reactant"), ("normal", "m_generic_prod2", "product"),
                       ("final", "m_generic_prod3", "product"), ("side", "m_generic_side1prime", "product")]:
    add_part(r, pid, mid, role, EV[5], product_role="proposed pathway species" if role == "product" else None)
set_class(r, EV[5], "unresolved")
r["notes"].append(f("possible mechanism is proposed", "Author proposal; intermediates A-E and Rh-H pathway are not treated as experimentally isolated species.", [EV[5]]))
r["review"]["unresolved"] = ["Elementary-step atom mapping and transient intermediate identities are not asserted."]

# Scheme 4a double epoxidations.
for suffix, start, product, yld in [("a", "m_prod3a_S", "m_prod4a", 61), ("c", "m_prod3c_S", "m_prod4c", 65),
                                     ("f", "m_prod3f_S", "m_prod4f", 56), ("h", "m_prod3h_S", "m_prod4h", 58)]:
    r = reaction(f"scheme4a_{suffix}", [EV[1], EV[5]] if suffix == "a" else [EV[5]], categories=["derivatization"],
                 location=f"page 5 Scheme 4a product 4{suffix}", name="Double epoxidation", transformation="mCPBA double epoxidation")
    add_part(r, "substrate", start, "reactant", EV[5])
    add_part(r, "mcpba", "m_mCPBA", "reagent", EV[5], charges={"equivalents": q(2.6, "equiv")}, ref=start)
    add_part(r, "bicarbonate", "m_NaHCO3", "additive", EV[5], charges={"equivalents": q(2.6, "equiv")}, ref=start)
    add_part(r, "dcm", "m_DCM", "solvent", EV[5])
    add_part(r, "water", "m_water", "solvent", EV[5])
    add_part(r, "product", product, "product", EV[5], product_role="double epoxide")
    set_solvent(r, EV[5], ["DCM", "H2O"], ratio="1:1")
    add_step(r, "reaction", EV[5], op="double epoxidation", temp=q(0, "°C"), time=q(24, "h"), pids=[p["participant_id"] for p in r["participants"]])
    add_outcome(r, "yield", "yield", [product], EV[5], pct(yld), raw=f"{yld}%")
    add_outcome(r, "dr", "selectivity", [product], EV[5], ">20:1", scope={"metric": "dr"}, raw="dr >20:1")
    set_class(r, EV[5], "epoxidation")
    r["review"]["unresolved"] = ["Yield and dr measurement bases are not stated; primary structures remain unencoded."]

# Scheme 4b first bis-dihydroxylation.
r = reaction("scheme4b_dihydroxylation", [EV[1], EV[5]], categories=["derivatization"], location="page 5 Scheme 4b, central green products",
             name="Bishydroxylation", transformation="OsO4/NMO bishydroxylation")
add_part(r, "substrate", "m_prod3a_S", "reactant", EV[5])
add_part(r, "osmium", "m_OsO4", "catalyst", EV[5], charges={"mol_percent": q(5, "mol%")}, ref="m_prod3a_S")
add_part(r, "nmo", "m_NMO", "reagent", EV[5], charges={"equivalents": q(1.2, "equiv")}, ref="m_prod3a_S")
add_part(r, "tsnh2", "m_TsNH2", "reagent", EV[5], charges={"equivalents": q(1.2, "equiv")}, ref="m_prod3a_S")
add_part(r, "acetone", "m_acetone", "solvent", EV[5])
add_part(r, "p5", "m_prod5a", "product", EV[5], product_role="regioisomeric dihydroxylation product")
add_part(r, "p6", "m_prod6a", "product", EV[5], product_role="regioisomeric dihydroxylation product")
set_solvent(r, EV[5], "acetone")
add_step(r, "reaction", EV[5], op="bishydroxylation", temp=q(0, "°C"), time=q(12, "h"), pids=[p["participant_id"] for p in r["participants"]])
add_outcome(r, "yield_5a", "yield", ["m_prod5a"], EV[5], pct(50), raw="5a, 50%")
add_outcome(r, "yield_6a", "yield", ["m_prod6a"], EV[5], pct(25), raw="6a, 25%")
set_class(r, EV[5], "dihydroxylation")
r["review"]["unresolved"] = ["Product-yield measurement bases are not stated; exact regioisomer structures remain unencoded."]

# Scheme 4b DMP oxidation of the regioisomer mixture.
r = reaction("scheme4b_DMP_oxidation", [EV[1], EV[5]], categories=["derivatization"], location="page 5 Scheme 4b upper branch",
             name="Dess-Martin oxidation of dihydroxylation mixture", transformation="selective alcohol oxidation")
add_part(r, "p5_input", "m_prod5a", "reactant", EV[5])
add_part(r, "p6_input", "m_prod6a", "reactant", EV[5])
add_part(r, "dmp", "m_DMP", "reagent", EV[5], charges={"equivalents": q(2.5, "equiv")}, ref="m_prod6a")
add_part(r, "bicarbonate", "m_NaHCO3", "additive", EV[5], charges={"equivalents": q(10, "equiv")}, ref="m_prod6a")
add_part(r, "dcm", "m_DCM", "solvent", EV[5])
add_part(r, "p7", "m_prod7a", "product", EV[5], product_role="alpha-hydroxyl ketone")
set_solvent(r, EV[5], "DCM")
add_step(r, "reaction", EV[5], op="DMP oxidation of the 5a/6a mixture", temp="r.t.", time=q(18, "h"), pids=[p["participant_id"] for p in r["participants"]])
add_outcome(r, "recovery_5a", "recovery", ["m_prod5a"], EV[5], pct(61), raw="5a, 61% recovered")
add_outcome(r, "yield_7a", "yield", ["m_prod7a"], EV[5], "quantitative", raw="7a, quant. yield")
add_outcome(r, "dr_7a", "selectivity", ["m_prod7a"], EV[5], ">20:1", scope={"metric": "dr"}, raw="dr >20:1")
set_class(r, EV[5], "oxidation")
r["review"]["unresolved"] = ["Yield/recovery/dr measurement bases and the exact amount of each regioisomer charged are unreported."]

# Scheme 4b acid treatment of the mixture.
r = reaction("scheme4b_TsOH_treatment", [EV[1], EV[5]], categories=["derivatization"], location="page 5 Scheme 4b lower branch",
             name="TsOH treatment of dihydroxylation mixture", transformation="ketalization and dehydration")
add_part(r, "p5_input", "m_prod5a", "reactant", EV[5])
add_part(r, "p6_input", "m_prod6a", "reactant", EV[5])
add_part(r, "tsoh", "m_TsOH", "catalyst", EV[5], charges={"mol_percent": q(10, "mol%")}, ref="m_prod5a")
add_part(r, "acetone", "m_acetone", "solvent", EV[5])
add_part(r, "p8", "m_prod8a", "product", EV[5], product_role="ketalization product from 5a")
add_part(r, "p9", "m_prod9a", "product", EV[5], product_role="dehydration product from 6a")
set_solvent(r, EV[5], "acetone")
add_step(r, "reaction", EV[5], op="TsOH treatment of the 5a/6a mixture", temp="r.t.", time=q(2, "h"), pids=[p["participant_id"] for p in r["participants"]])
add_outcome(r, "yield_8a", "yield", ["m_prod8a"], EV[5], "quantitative", raw="8a, quant. yield")
add_outcome(r, "dr_8a", "selectivity", ["m_prod8a"], EV[5], ">20:1", scope={"metric": "dr"}, raw="dr >20:1")
add_outcome(r, "yield_9a", "yield", ["m_prod9a"], EV[5], pct(60), raw="9a, 60%")
add_outcome(r, "dr_9a", "selectivity", ["m_prod9a"], EV[5], ">20:1", scope={"metric": "dr"}, raw="dr >20:1")
set_class(r, EV[5], "ketal_dehydration")
r["review"]["unresolved"] = ["Outcome measurement bases and exact amounts of each regioisomer charged are unreported."]

bundle["molecules"] = list(molecules.values())
bundle["source_context"] = [
    {"context_id": f"{DOC}:page6_3a_characterization", "document_id": DOC, "evidence_ids": [EV[5], EV[6]],
     "kind": "characterization", "description": "3a NMR/IR/EI-MS/HRMS, HPLC retention, optical rotation and physical form; not a separate synthesis."},
    {"context_id": f"{DOC}:page3_xray_3k", "document_id": DOC, "evidence_ids": [EV[3]],
     "kind": "structure_assignment", "description": "ORTEP of (3aS,7aS,E)-3k, CCDC 1953833; supports the paper's absolute-configuration assignment."},
    {"context_id": f"{DOC}:page5_xray_products", "document_id": DOC, "evidence_ids": [EV[5]],
     "kind": "structure_assignment", "description": "X-ray depictions for 4a (CCDC 1953831) and a dihydroxylation derivative (CCDC 1953832); not separate syntheses."},
    {"context_id": f"{DOC}:page6_si_availability", "document_id": DOC, "evidence_ids": [EV[6]],
     "kind": "supplementary_information_statement", "description": "Article states SI is available online; no SI file was supplied or explicitly associated in this workspace."},
]

(RUN / "submitted.json").write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")

# Raw visual inventory and exact requested page coverage.
all_ids = [r["reaction_id"] for r in bundle["reactions"]]
def ids(prefix):
    return [x for x in all_ids if x.split(":", 1)[1].startswith(prefix)]

coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {"page": 1, "state": "attempted", "regions": [
            {"region": "title/abstract/keywords", "finding": "article metadata and summary; no specific experiment", "evidence_ids": [EV[1]]},
            {"region": "Comprehensive Graphic Content", "finding": "generic current-work reaction and overview of products later detailed in Scheme 4", "evidence_ids": [EV[1]]},
        ], "record_ids": ids("overview_generic") + ids("scheme4"), "unresolved": ["Overview contains generic R groups and no experiment-specific yields/scale/time."]},
        {"page": 2, "state": "attempted", "regions": [
            {"region": "Scheme 1a", "finding": "biological structures only; not synthetic experiments in this paper", "evidence_ids": [EV[2]]},
            {"region": "Scheme 1b-i and 1b-ii", "finding": "one prior-literature generic reaction and one current-work generic concept", "evidence_ids": [EV[2]]},
            {"region": "Table 1 rows 1-7 and footnotes a-f", "finding": "seven ligand-screen experiments, including trace/zero products, recovery, ee and entry-7 tetra-ene side product", "evidence_ids": [EV[2]]},
        ], "record_ids": ids("scheme1_") + ids("overview_generic") + ids("table1_"), "unresolved": ["Primary structures are not safely encoded; Table 1 signed ee reference is retained verbatim."]},
        {"page": 3, "state": "attempted", "regions": [
            {"region": "Table 2 rows 1-9 and footnotes", "finding": "nine independent catalyst/silver-salt/solvent optimization experiments", "evidence_ids": [EV[3]]},
            {"region": "Table 3 rows 1-14, both ligand columns and footnotes", "finding": "28 independent scope experiments (14 with (R,R)-Ph-BPE and 14 with (S,S)-Ph-BPE)", "evidence_ids": [EV[3]]},
            {"region": "Figure 1", "finding": "3k ORTEP/CCDC structure assignment, not a synthesis", "evidence_ids": [EV[3]]},
            {"region": "mechanistic discussion", "finding": "textual context for Scheme 2 experiments on the next page", "evidence_ids": [EV[3]]},
        ], "record_ids": ids("table2_") + ids("table3_"), "unresolved": ["Compact Table 3 structures are retained by label/R-group assignment rather than guessed SMILES."]},
        {"page": 4, "state": "attempted", "regions": [
            {"region": "Scheme 2a-b", "finding": "three ligand-dependence controls and one racemic-2a isomerization experiment", "evidence_ids": [EV[4]]},
            {"region": "Scheme 2c table rows 1-4", "finding": "four solvent controls with NMR/HPLC outcomes", "evidence_ids": [EV[4]]},
            {"region": "Scheme 2d table rows 1-4", "finding": "four time-course experiments", "evidence_ids": [EV[4]]},
            {"region": "Scheme 2e-f", "finding": "two enantiopure-substrate experiments and two EtOH-d6 isotope experiments", "evidence_ids": [EV[4]]},
            {"region": "Scheme 2g", "finding": "two real failed tether controls; no expected products fabricated", "evidence_ids": [EV[4]]},
        ], "record_ids": ids("scheme2"), "unresolved": ["Several Scheme 2 panels omit measurement bases; these remain unreported rather than borrowed from other panels."]},
        {"page": 5, "state": "attempted", "regions": [
            {"region": "Scheme 3", "finding": "author-proposed mechanism, recorded separately from experiments", "evidence_ids": [EV[5]]},
            {"region": "Scheme 4a", "finding": "four specific double-epoxidation experiments", "evidence_ids": [EV[5]]},
            {"region": "Scheme 4b", "finding": "one bishydroxylation and two subsequent mixture-transformations", "evidence_ids": [EV[5]]},
            {"region": "Experimental: asymmetric synthesis of 3a", "finding": "detailed 1.0 mmol procedure merged with Table 3 entry 4 left-column record", "evidence_ids": [EV[5]]},
        ], "record_ids": ids("scheme3_") + ids("scheme4") + [f"{DOC}:table3_entry4_S"], "unresolved": ["Scheme 4 outcome measurement bases are not printed; no SI procedure was associated."]},
        {"page": 6, "state": "attempted", "regions": [
            {"region": "top-left continuation of 3a characterization", "finding": "NMR/IR/EI-MS/HRMS continuation; characterization context, not a new synthesis", "evidence_ids": [EV[6]]},
            {"region": "Supporting Information statement", "finding": "SI availability statement only; no SI file supplied or associated", "evidence_ids": [EV[6]]},
            {"region": "Acknowledgement and references", "finding": "publication back matter and cited prior work; no current-paper experiments", "evidence_ids": [EV[6]]},
        ], "record_ids": [f"{DOC}:table3_entry4_S"], "unresolved": ["The detailed procedure continues only as characterization on this page; no additional synthetic arrow/table is present."]},
    ],
}
(RUN / "coverage").mkdir(parents=True, exist_ok=True)
(RUN / "coverage/pages_0001_0006.json").write_text(json.dumps(coverage, ensure_ascii=False, indent=2) + "\n")

inventory = {
    "document_id": DOC, "pages_viewed_at_original_resolution": [1, 2, 3, 4, 5, 6],
    "crops": ["page2_table1.png", "page3_tables2_3.png", "page4_scheme2.png", "page5_schemes3_4_and_experimental.png"],
    "notes": "Raw visual inventory only; chemical decisions are in submitted.json, coverage, and runs/workflow.md.",
}
(RUN / "raw_model/page1_6_visual_inventory.json").write_text(json.dumps(inventory, ensure_ascii=False, indent=2) + "\n")

# Diagnostics for every concrete SMILES candidate.
diagnostics = {mid: st.stereo_inspect(st.value(m["smiles"])) for mid, m in molecules.items() if m["representation_type"] == "smiles"}
(RUN / "raw_model/stereo_diagnostics.json").write_text(json.dumps(diagnostics, ensure_ascii=False, indent=2) + "\n")

print(json.dumps({"molecules": len(bundle["molecules"]), "reactions": len(bundle["reactions"]),
                  "specific_experiments": sum(r["record_type"] == "specific_experiment" for r in bundle["reactions"]),
                  "generic_or_contextual_reactions": sum(r["record_type"] != "specific_experiment" for r in bundle["reactions"])}, ensure_ascii=False))
