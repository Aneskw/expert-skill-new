#!/usr/bin/env python3
"""Build the source-derived v3.1 submission for doc_923862b0df8c4625.

This is a machine extraction.  It intentionally leaves dense structures,
unreported procedures, measurement bases, and atom mapping unresolved.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

DOC = "doc_923862b0df8c4625"
P = {n: f"{DOC}:page{n}:text" for n in range(1, 5)}
I = {n: f"{DOC}:page{n}:image" for n in range(1, 5)}


def f(raw=None, normalized=None, ev=(), method="original_extraction", derivation=None,
      status=None):
    return st.fact(raw, normalized, ev, method=method, derivation=derivation,
                   status=status)


def inferred(value, ev, derivation):
    return f(None, value, ev, method="inference", derivation=derivation)


def derived(value, ev, derivation):
    return f(None, value, ev, method="structure_derivation", derivation=derivation)


def molecule(mid, label, pages, name=None, smiles=None, stereo_note=None,
             chirality="unresolved", geometry="unresolved", axial_note=None):
    ev = [I[p] for p in pages]
    m = st.new_molecule(mid, DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"] = [{
        "document_id": DOC,
        "compound_label": f(label, label, ev),
        "evidence_ids": ev,
    }]
    display = name or label
    m["names"]["original"] = f(display, display, ev)
    if name:
        m["names"]["english"] = f(name, name, ev)
    m["cas_numbers"] = f()
    m["identity_note"] = f(
        f"Source identity: {display}; no CAS number is reported in pages 1–4.",
        {"source_label": label, "name": name, "cas_reported": False}, ev)
    if smiles:
        diag = st.stereo_inspect(smiles)
        canonical = diag["canonical_isomeric_smiles"]
        m["representation_type"] = "smiles"
        m["smiles"] = derived(
            canonical, ev,
            "Candidate SMILES for an explicitly named simple reagent/solvent; RDKit parsing is not source-structure verification.")
        m["stereochemistry"]["atom_centers"] = derived(
            diag["atom_centers"], ev, "RDKit graph diagnostic on the candidate SMILES.")
        m["stereochemistry"]["double_bonds"] = derived(
            diag["double_bonds"], ev, "RDKit graph diagnostic on the candidate SMILES.")
        potential = diag["potential_stereo"]
        has_atom = any(x["type"].startswith("Atom_") for x in potential)
        has_bond = any(x["type"].startswith("Bond_") for x in potential)
        m["stereochemistry"]["chirality_status"] = derived(
            "unspecified" if has_atom else "achiral", ev,
            "Status is limited to the candidate molecular graph; absence of a toolkit flag is not used for dense source structures.")
        m["stereochemistry"]["geometry_status"] = derived(
            "unspecified" if has_bond else "not_applicable", ev,
            "Status is limited to configurable double bonds in the candidate molecular graph.")
        m["stereochemistry"]["representation_status"] = derived(
            "partial" if (has_atom or has_bond) else "encoded", ev,
            "Candidate SMILES was inspected with standard stereo; any unassigned potential element prevents an encoded claim.")
        m["stereochemistry"]["isomeric_smiles"] = derived(
            canonical, ev, "Canonical isomeric SMILES returned by standard stereo diagnostics.")
        m["stereochemistry"]["identity_key"] = derived(
            diag["stereo_graph_key"], ev, "Stereo graph key returned by standard stereo diagnostics.")
    else:
        m["representation_type"] = "unresolved"
        if stereo_note:
            m["stereochemistry"]["atom_centers"] = f(stereo_note, stereo_note, ev)
            m["stereochemistry"]["assignment"] = f(stereo_note, stereo_note, ev)
        if axial_note:
            m["stereochemistry"]["axial"] = f(axial_note, axial_note, ev)
        m["stereochemistry"]["chirality_status"] = f(
            stereo_note or chirality, chirality, ev,
            method="original_extraction" if stereo_note else "inference",
            derivation=None if stereo_note else "The source label/drawing does not establish a fully encoded stereochemical identity.")
        m["stereochemistry"]["geometry_status"] = f(
            geometry if geometry != "unresolved" else "not fully encoded from source drawing",
            geometry, ev, method="inference",
            derivation="E/Z and ring-relative geometry were reviewed visually; the dense drawing is retained without a concrete graph.")
        m["stereochemistry"]["representation_status"] = inferred(
            "not_encoded", ev,
            "No concrete SMILES is asserted because the dense drawing or exact stereochemical graph was not safely encoded.")
    return m


def set_charge(participant, key, raw, value, ev, reference=None,
               method="original_extraction", derivation=None):
    participant["charges"][key] = f(raw, value, ev, method=method,
                                     derivation=derivation)
    if key in ("equivalents", "mol_percent"):
        participant["equivalence_reference"] = inferred(
            reference, ev,
            "The loading/equivalent denominator is the stated limiting or reference substrate in this experiment.")


def participant(pid, mid, role, ev, product_role=None, **charges):
    p = st.new_participant(pid, mid, role, ev)
    if product_role:
        p["product_role"] = f(product_role, product_role, ev)
    for key, spec in charges.items():
        raw, value, reference = spec
        set_charge(p, key, raw, value, ev, reference)
    return p


def classification(ev, changes, note):
    bond_changes = []
    labels = {}
    for idx, (a, b, before, after, kind) in enumerate(changes, 1):
        cid = f"b{idx}"
        change = {
            "change_id": cid,
            "atoms": [{"atom_map": None, "element": a},
                      {"atom_map": None, "element": b}],
            "before_order": before,
            "after_order": after,
            "change_type": kind,
        }
        bond_changes.append(derived(change, ev, note))
        if kind in ("formed", "hydrogen_gain", "order_increased"):
            lab = st.bond_label(a, b, after)
            labels.setdefault(lab["code"], dict(lab, supporting_change_ids=[]))[
                "supporting_change_ids"].append(cid)
    return {
        "standard_id": st.STANDARD,
        "mapped_reaction_smiles": f(),
        "bond_changes": bond_changes,
        "labels": derived(list(labels.values()), ev, note),
        "coverage": inferred("partial", ev,
                             "Category-level net changes were derived from source drawings; no atom-complete mapping was asserted."),
        "limitations": [
            "Category-level source-drawing interpretation only; atom identities are not mapped.",
            "Unreported coproducts and proton sources are outside the partial classification.",
        ],
        "model_votes": [],
    }


HYDRO = [("C", "O", 0, 1, "formed"), ("C", "H", 0, 1, "hydrogen_gain"),
         ("C", "C", 2, 1, "order_decreased"), ("O", "H", 1, 0, "broken")]
RCM = [("C", "C", 0, 2, "formed")]
DIOL = [("C", "O", 0, 1, "formed"), ("C", "O", 0, 1, "formed"),
        ("C", "C", 2, 1, "order_decreased"),
        ("O", "H", 0, 1, "hydrogen_gain"), ("O", "H", 0, 1, "hydrogen_gain")]


def reaction(suffix, ev, record_type, categories, location, name, transformation,
             parts, steps, outcomes, changes, procedure, unresolved,
             general=()):
    rid = f"{DOC}:{suffix}"
    r = st.new_reaction(rid, DOC, ev)
    r["record_type"] = record_type
    r["experiment_categories"] = f(location, categories, ev)
    r["location"] = f(location, {"source_location": location}, ev)
    r["general_scheme_ids"] = list(general)
    r["reaction_name"] = f(name, name, ev)
    r["transformation_type"] = derived(
        transformation, ev, "Transformation type derived from the source-described reagents and depicted connectivity change.")
    r["participants"] = parts
    r["steps"] = steps
    r["procedure"]["original_text"] = f(procedure, procedure, ev)
    r["outcomes"] = outcomes
    r["classification"] = classification(
        ev, changes,
        "Partial net-bond assignment from the depicted starting material and product; no atom mapping was supplied.")
    r["evidence_ids"] = list(dict.fromkeys(ev))
    r["review"] = {
        "status": "needs_review", "method": "model", "reviewer": None,
        "checked_fields": [],
        "unresolved": list(unresolved) + [
            "Supporting Information was mentioned but not supplied or explicitly associated.",
            "Atom-complete structure mapping and bond-change verification remain unresolved.",
        ],
    }
    r["completeness"] = "partial"
    return r


def step(sid, ev, pids, operation, temp=None, time=None):
    s = st.new_step(sid, ev)
    s["participant_ids"] = pids
    s["operation"] = f(operation, operation, ev)
    if temp is not None:
        s["temperature"] = f(temp[0], temp[1], ev)
    if time is not None:
        s["time"] = f(time[0], time[1], ev)
    return s


def outcome(oid, kind, targets, ev, raw, value, role=None, yield_type=None,
            basis=None, scope=None):
    o = st.new_outcome(oid, kind, targets, ev)
    o["value"] = f(raw, value, ev)
    if role:
        o["product_role"] = f(role, role, ev)
    if yield_type:
        o["yield_type"] = f(yield_type, yield_type, ev)
    if basis:
        o["measurement_basis"] = f(basis, basis, ev)
    if scope:
        o["scope"] = f(scope, scope, ev)
    return o


dataset = json.loads((ROOT / "runs/current/dataset.json").read_text())
bundle = st.empty(dataset["documents"])

# Source-labelled dense structures.  Wedge/dash information is preserved as text;
# exact graphs are intentionally not guessed.
dense = [
    ("m_scheme_A", "A", [1], "generic homoallylic alcohol A", None, "unspecified", "not_applicable"),
    ("m_scheme_B", "B", [1], "generic ene-alkoxyallene B", "allene C=C=C is drawn; no axial descriptor is given", "unspecified", "not_applicable"),
    ("m_scheme_C", "C", [1], "generic dihydropyran intermediate C", None, "specified", "specified"),
    ("m_scheme_D", "D", [1], "generic saccharosamine–rhamnose product D", None, "specified", "specified"),
    ("m_scheme_E", "E", [1], "generic unprotected β-D-fucose E", None, "specified", "specified"),
    ("m_scheme_F", "F", [1], "generic ene-alkoxyallene F", "allene C=C=C is drawn; no axial descriptor is given", "partially_specified", "not_applicable"),
    ("m_scheme_G", "G", [1], "generic rhamnose–fucose product G", None, "specified", "specified"),
    ("m1", "1", [2], "alcohol 1", None, "specified", "not_applicable"),
    ("m2", "2", [2], "alkoxyallene 2", "allene C=C=C is drawn; no axial descriptor is given", "unspecified", "not_applicable"),
    ("m4a", "4a", [2], "acyclic O,O-acetal 4a", None, "specified", "not_applicable"),
    ("m4b", "4b", [2], "acyclic O,O-acetal 4b", None, "specified", "not_applicable"),
    ("m5", "5", [2, 3], "ene-alkoxyallene 5", "allene C=C=C is drawn; no axial descriptor is given", "partially_specified", "not_applicable"),
    ("m7a", "7a", [2], "dihydropyran 7a", None, "specified", "specified"),
    ("m7b", "7b", [2], "dihydropyran 7b", None, "specified", "specified"),
    ("m8", "8", [2], "alcohol 8", None, "specified", "specified"),
    ("m9", "9", [2], "carbamate 9", None, "specified", "specified"),
    ("m10", "10", [2], "β-D-saccharosamine precursor 10", None, "specified", "specified"),
    ("m11", "11", [2, 3], "allylic alcohol 11", None, "specified", "specified"),
    ("m12a", "12a", [2, 3], "acyclic O,O-acetal 12a", None, "specified", "specified"),
    ("m12b", "12b", [2, 3], "acyclic O,O-acetal 12b", None, "specified", "specified"),
    ("m13a", "13a", [2, 3], "biscyclic intermediate 13a", None, "specified", "specified"),
    ("m13b", "13b", [2, 3], "biscyclic intermediate 13b", None, "specified", "specified"),
    ("m14", "14", [2, 3], "diol 14", None, "specified", "specified"),
    ("m15", "15", [2, 3], "MOM-protected compound 15", None, "specified", "specified"),
    ("m15oh", "unnumbered alcohol from 15", [2, 3], "debenzylated alcohol from 15", None, "specified", "specified"),
    ("m16", "16", [2, 3], "carbamate 16", None, "specified", "specified"),
    ("m17", "17", [2, 3], "saccharosamine–rhamnose disaccharide 17", None, "specified", "specified"),
    ("m18", "18", [3], "protected β-D-fucoside 18", None, "specified", "specified"),
    ("m18_adduct", "trace adduct from 18 and 19", [3], "unassigned trace adduct", None, "unresolved", "unresolved"),
    ("m19", "19", [3], "ene-alkoxyallene 19", "allene C=C=C is drawn; no axial descriptor is given", "partially_specified", "not_applicable"),
    ("m20", "20", [3], "unprotected β-D-fucose triol 20", None, "specified", "specified"),
    ("m21a", "21a", [3], "fucose–rhamnose disaccharide 21a", None, "specified", "specified"),
    ("m21b", "21b", [3], "fucose–rhamnose disaccharide 21b", None, "specified", "specified"),
    ("m21regio", "regioisomers of 21", [3], "unassigned regioisomer mixture of 21", None, "mixture", "mixture"),
    ("m22", "22", [3], "bisbenzoate 22", None, "specified", "specified"),
    ("m23", "23", [3], "tetrol 23", None, "specified", "specified"),
    ("lig_ss", "(S,S)-3", [2, 3], "(S,S)-Trost-type ligand 3", None, "specified", "specified"),
    ("lig_rr", "(R,R)-3", [2, 3], "(R,R)-Trost-type ligand 3", None, "specified", "specified"),
]
for mid, label, pages, name, axial, chiral, geom in dense:
    stereo = ("Source drawing specifies stereogenic centers/relative configuration"
              if chiral in ("specified", "partially_specified", "mixture") else None)
    bundle["molecules"].append(molecule(mid, label, pages, name=name,
                                        stereo_note=stereo, chirality=chiral,
                                        geometry=geom, axial_note=axial))

simple = [
    ("m6", "6", [2], "cyclohexanol", "OC1CCCCC1"),
    ("pd2dba3", "Pd2(dba)3", [2, 3], "tris(dibenzylideneacetone)dipalladium(0)", None),
    ("et3n", "Et3N", [2, 3], "triethylamine", "CCN(CC)CC"),
    ("toluene", "toluene", [2], "toluene", "Cc1ccccc1"),
    ("ru_cat", "Ru cat.", [2, 3], "second-generation Hoveyda–Grubbs catalyst", None),
    ("na", "Na", [2, 3], "sodium", "[Na]"),
    ("nh3", "NH3", [2, 3], "ammonia", "N"),
    ("pmp_nco", "PMP–NCO", [2, 3], "p-methoxyphenyl isocyanate", "COc1ccc(N=C=O)cc1"),
    ("dbu", "DBU", [2, 3], "1,8-diazabicyclo[5.4.0]undec-7-ene", "C1CCC2=NCCCN2CC1"),
    ("ch2cl2", "CH2Cl2", [2, 3], "dichloromethane", "ClCCl"),
    ("ibx", "IBX", [2, 3], "2-iodoxybenzoic acid", None),
    ("nahco3", "NaHCO3", [2, 3], "sodium bicarbonate", "[Na+].[O-]C(=O)O"),
    ("dmso", "DMSO", [2, 3], "dimethyl sulfoxide", "CS(C)=O"),
    ("thf", "THF", [2, 3], "tetrahydrofuran", "C1CCOC1"),
    ("oso4", "OsO4", [2, 3], "osmium tetroxide", None),
    ("nmo", "NMO", [2, 3], "N-methylmorpholine N-oxide", "C[N+]1(CCOCC1)[O-]"),
    ("momcl", "MOMCl", [2, 3], "methoxymethyl chloride", "COCCl"),
    ("dipea", "N,N-diisopropylethylamine", [3], "N,N-diisopropylethylamine", "CCN(C(C)C)C(C)C"),
    ("benzoyl_chloride", "benzoyl chloride", [3], "benzoyl chloride", "O=C(Cl)c1ccccc1"),
    ("pyridine", "pyridine", [3], "pyridine", "n1ccccc1"),
    ("acetone", "acetone", [3], "acetone", "CC(C)=O"),
]
for mid, label, pages, name, smiles in simple:
    bundle["molecules"].append(molecule(mid, label, pages, name=name, smiles=smiles))

M = {m["molecule_id"]: m for m in bundle["molecules"]}


def parts_for(specs, ev):
    parts = []
    for idx, spec in enumerate(specs, 1):
        mid, role, product_role, charges = spec
        parts.append(participant(f"p{idx}", mid, role, ev,
                                 product_role=product_role, **charges))
    return parts


def q(number, unit):
    return {"number": number, "unit": unit}


GS_A = f"{DOC}:scheme1_A_generic"
GS_B = f"{DOC}:scheme1_B_generic"

# Scheme 1 current-work concepts.  The graphical abstract is duplicate evidence,
# not another experiment.
ev = [I[1], P[1]]
pa = parts_for([
    ("m_scheme_A", "reactant", None, {}), ("m_scheme_B", "reactant", None, {}),
    ("pd2dba3", "catalyst", None, {}), ("lig_ss", "ligand", None, {}),
    ("ru_cat", "catalyst", None, {}), ("m_scheme_C", "product", "intermediate", {}),
    ("m_scheme_D", "product", "target", {}),
], ev)
bundle["reactions"].append(reaction(
    "scheme1_A_generic", ev, "generic_scheme", ["other"], "Page 1 Scheme 1A and graphical abstract",
    "activation-free synthesis of β-(D)-saccharoside", "Pd-catalyzed hydroalkoxylation followed by RCM and oxidative cyclization",
    pa, [], [], HYDRO + RCM + [("C", "N", 0, 1, "formed")],
    "A + B → C → D; i) Pd/chiral ligand; ii) Ru (RCM)",
    ["Generic labels A–D and catalyst identities do not define a reproducible experiment."], ()))
pb = parts_for([
    ("m_scheme_E", "reactant", None, {}), ("m_scheme_F", "reactant", None, {}),
    ("pd2dba3", "catalyst", None, {}), ("lig_ss", "ligand", None, {}),
    ("ru_cat", "catalyst", None, {}), ("oso4", "reagent", None, {}),
    ("m_scheme_G", "product", "target", {}),
], ev)
bundle["reactions"].append(reaction(
    "scheme1_B_generic", ev, "generic_scheme", ["other"], "Page 1 Scheme 1B and graphical abstract",
    "regioselective reaction of β-(D)-fucoside", "Pd-catalyzed hydroalkoxylation/RCM followed by dihydroxylation",
    pb, [], [], HYDRO + RCM + DIOL,
    "E + F → G; i) Pd/chiral ligand; ii) Ru; iii) OsO4",
    ["Generic labels E–G and abbreviated reagents do not define a reproducible experiment."], ()))


def hydro_parts(allene, alcohol, product, allene_eq, alcohol_eq, ligand, pd, lig,
                ev, et=0.1, ru=None, minor=None):
    specs = [
        (allene, "reactant", None, {"equivalents": (f"{allene_eq} equiv", q(allene_eq, "equiv"), alcohol if alcohol_eq == 1 else allene)}),
        (alcohol, "reactant", None, {"equivalents": (f"{alcohol_eq} equiv", q(alcohol_eq, "equiv"), alcohol if alcohol_eq == 1 else allene)}),
        ("pd2dba3", "catalyst", None, {"mol_percent": (f"{pd} mol %", q(pd, "mol%"), allene)}),
        (ligand, "ligand", None, {"mol_percent": (f"{lig} mol %", q(lig, "mol%"), allene)}),
        ("et3n", "additive", None, {"equivalents": (f"{et} equiv", q(et, "equiv"), allene)}),
        ("toluene", "solvent", None, {}),
        (product, "product", "major product", {}),
    ]
    if ru is not None:
        specs.insert(-1, ("ru_cat", "catalyst", None,
                          {"mol_percent": (f"{ru} mol %", q(ru, "mol%"), allene)}))
    if minor:
        specs.append((minor, "product", "minor stereoisomer", {}))
    return parts_for(specs, ev)


# Preliminary route 1 + 2 -> 4a.
ev = [I[2], P[2]]
parts = hydro_parts("m2", "m1", "m4a", 1, 1, "lig_ss", 2.5, 5, ev)
steps = [step("coupling", ev, [p["participant_id"] for p in parts],
              "Pd-catalyzed coupling of alcohol 1 and alkoxyallene 2",
              ("40 °C", q(40, "°C")), ("24 h", q(24, "h")))]
outs = [outcome("o1", "yield", ["m4a"], ev, "<20% yield",
                {"number": 20, "unit": "%", "qualifier": "lt", "raw_text": "<20%"},
                "adduct", scope="single coupling step")]
bundle["reactions"].append(reaction(
    "preliminary_1_2_to_4a", ev, "specific_experiment", ["optimization"],
    "Page 2 Table 1 scheme and opening optimization prose",
    "Pd-catalyzed asymmetric coupling", "allene hydroalkoxylation",
    parts, steps, outs, HYDRO,
    "Alcohol 1 and alkoxyallene 2 gave 4a in low (<20%) yield under the depicted standard conditions.",
    ["Yield measurement basis and experimental details referred to SI are unavailable."], (GS_A,)))

# Table 1 entries 1–5.
entries = [
    (1, "lig_ss", 1, 1, 2.5, 5, "m4a", 36, None, None, None),
    (2, "lig_ss", 2, 1, 2.5, 5, "m4a", 48, None, None, None),
    (3, "lig_ss", 1, 2, 2.5, 5, "m4a", 59, None, None, None),
    (4, "lig_ss", 1, 2, 5, 10, "m4a", 65, "m7a", 52, "8:1"),
    (5, "lig_rr", 1, 2, 3.5, 7, "m4b", 94, "m7b", 90, "<1:25"),
]
for ent, ligand, ae, oe, pd, lig, first, y1, final, y2, ratio in entries:
    ev = [I[2], P[2]]
    minor = ("m7b" if final == "m7a" else "m7a") if final else None
    parts = hydro_parts("m5", "m6", first, ae, oe, ligand, pd, lig, ev,
                        ru=5 if final else None, minor=minor)
    if final:
        parts.append(participant(f"p{len(parts)+1}", final, "product", ev,
                                 product_role="isolated major second-step product"))
    pids = [p["participant_id"] for p in parts]
    steps = [step("coupling", ev, pids,
                  "Pd-catalyzed coupling of ene-alkoxyallene 5 and alcohol 6",
                  ("40 °C", q(40, "°C")), ("24 h", q(24, "h")))]
    changes = list(HYDRO)
    if final:
        steps.append(step("rcm", ev, pids,
                          "Subsequent ring-closing metathesis with second-generation Hoveyda–Grubbs catalyst"))
        changes += RCM
    y1_basis = "isolated yield" if ent in (4, 5) else "integration of the crude 1H NMR"
    outs = [outcome("o1", "yield", [first], ev, f"{first[1:]} ({y1}%)",
                    q(y1, "%"), "first-step product",
                    yield_type="isolated" if ent in (4, 5) else "crude NMR",
                    basis=y1_basis, scope="first Pd-catalyzed coupling step")]
    if final:
        outs.append(outcome("o2", "yield", [final], ev, f"{final[1:]} ({y2}%)",
                            q(y2, "%"), "pure major stereoisomer",
                            yield_type="isolated", basis="isolated yield of the pure major stereoisomer",
                            scope="overall yield over Pd coupling and RCM"))
        outs.append(outcome("o3", "selectivity", ["m7a", "m7b"], ev,
                            f"ratio of 7a:7b = {ratio}", ratio,
                            scope="diastereomer ratio ordered 7a:7b; crude 1H NMR integration"))
        if ent == 4:
            # Table footnote e supplies the only explicit scale in the optimization table.
            pass
    r = reaction(
        f"table1_entry{ent}", ev, "specific_experiment", ["optimization"] + (["scale_up"] if ent == 4 else []),
        f"Page 2 Table 1 entry {ent}", "Pd-catalyzed asymmetric coupling" + (" and RCM" if final else ""),
        "allene hydroalkoxylation" + (" followed by ring-closing metathesis" if final else ""),
        parts, steps, outs, changes,
        f"Table 1 entry {ent}: ligand, 5/6 ratio, first-step product and any second-step product as tabulated.",
        ["Workup and purification are unreported in the main-text pages."], (GS_A,))
    if ent == 4:
        r["reaction_scale"] = f("1 mmol scale", {"reference_molecule_id": "m5", "amount": q(1, "mmol")}, ev)
    bundle["reactions"].append(r)

# Scheme 2, three separate arrows.
ev = [I[2], P[2]]
specs = [("m7a", "reactant", None, {}), ("na", "reagent", None, {}),
         ("nh3", "solvent", None, {}), ("m8", "product", "alcohol product", {})]
parts = parts_for(specs, ev)
bundle["reactions"].append(reaction(
    "scheme2_7a_to_8", ev, "specific_experiment", ["other"], "Page 2 Scheme 2, arrow i",
    "benzyl-group removal", "Na/NH3 debenzylation", parts,
    [step("debenzylation", ev, [p["participant_id"] for p in parts], "Treatment with Na/NH3")],
    [outcome("o1", "yield", ["m8"], ev, "93%", q(93, "%"), "alcohol 8", scope="7a to 8")],
    [("O", "H", 0, 1, "hydrogen_gain"), ("C", "O", 1, 0, "broken")],
    "The benzyl group in 7a was removed by Na/NH3 to give alcohol 8.",
    ["Temperature, time, workup, purification, and yield basis are unreported in the main text."], (GS_A,)))

specs = [("m8", "reactant", None, {}),
         ("pmp_nco", "reagent", None, {"equivalents": ("1.1 equiv", q(1.1, "equiv"), "m8")}),
         ("dbu", "reagent", None, {"equivalents": ("0.1 equiv", q(0.1, "equiv"), "m8")}),
         ("ch2cl2", "solvent", None, {}), ("m9", "product", "carbamate product", {})]
parts = parts_for(specs, ev)
bundle["reactions"].append(reaction(
    "scheme2_8_to_9", ev, "specific_experiment", ["other"], "Page 2 Scheme 2, arrow ii",
    "carbamate formation", "alcohol addition to aryl isocyanate", parts,
    [step("carbamoylation", ev, [p["participant_id"] for p in parts],
          "Treatment with p-methoxyphenyl isocyanate and DBU in CH2Cl2",
          ("rt", "rt"))],
    [outcome("o1", "yield", ["m9"], ev, "90%", q(90, "%"), "carbamate 9", scope="8 to 9")],
    [("C", "O", 0, 1, "formed"), ("N", "H", 0, 1, "hydrogen_gain"),
     ("C", "N", 2, 1, "order_decreased"), ("O", "H", 1, 0, "broken")],
    "Alcohol 8 was reacted with PMP–NCO and DBU in CH2Cl2 at rt.",
    ["Reaction time, workup, purification, and yield basis are unreported in the main text."], (GS_A,)))

specs = [("m9", "reactant", None, {}),
         ("ibx", "reagent", None, {"equivalents": ("6 equiv", q(6, "equiv"), "m9")}),
         ("nahco3", "reagent", None, {"equivalents": ("12 equiv", q(12, "equiv"), "m9")}),
         ("thf", "solvent", None, {}), ("dmso", "solvent", None, {}),
         ("m10", "product", "single-diastereomer target", {})]
parts = parts_for(specs, ev)
r = reaction(
    "scheme2_9_to_10", ev, "specific_experiment", ["other"], "Page 2 Scheme 2, arrow iii",
    "IBX-mediated oxidative cyclization", "oxidative intramolecular C–N bond formation", parts,
    [step("oxidative_cyclization", ev, [p["participant_id"] for p in parts],
          "Heat with IBX and NaHCO3 in THF:DMSO in a sealed tube",
          ("90 °C", q(90, "°C")))],
    [outcome("o1", "yield", ["m10"], ev, "61%", q(61, "%"), "target 10", scope="9 to 10"),
     outcome("o2", "selectivity", ["m10"], ev, "single diastereomer", "single diastereomer",
             scope="diastereoselectivity; measurement method not reported")],
    [("C", "N", 0, 1, "formed"), ("C", "C", 2, 1, "order_decreased")],
    "Carbamate 9 was heated with IBX and NaHCO3 in DMSO:THF (10:1) at 90 °C.",
    ["Reaction time, workup, purification, and yield basis are unreported in the main text."], (GS_A,))
r["solvent_conditions"]["composition"] = f("DMSO:THF", ["DMSO", "THF"], ev)
r["solvent_conditions"]["mixture_ratio"] = f("DMSO:THF (10:1)", {"DMSO": 10, "THF": 1}, ev)
r["special_conditions"]["reactor"] = f("sealed tube", "sealed tube", ev)
bundle["reactions"].append(r)

# Scheme 3 coupling/RCM pairs.
def scheme3_coupling(suffix, ligand, allene_eq, alcohol_eq, product, y):
    ev = [I[3], P[2]]
    specs = [
        ("m5", "reactant", None, {"equivalents": (f"{allene_eq} equiv", q(allene_eq, "equiv"), "m11")}),
        ("m11", "reactant", None, {"equivalents": (f"{alcohol_eq} equiv", q(alcohol_eq, "equiv"), "m11")}),
        ("pd2dba3", "catalyst", None, {"mol_percent": ("3.5 mol %", q(3.5, "mol%"), "m5")}),
        (ligand, "ligand", None, {"mol_percent": ("7 mol %", q(7, "mol%"), "m5")}),
        ("et3n", "additive", None, {"equivalents": ("0.1 equiv", q(0.1, "equiv"), "m5")}),
        ("m12a" if product == "m12a" else "m12b", "product", "single stereoisomer", {}),
    ]
    parts = parts_for(specs, ev)
    outcomes = [outcome(
        "o1", "yield", [product], ev, f"{y}%", q(y, "%"), "O,O-acetal intermediate",
        yield_type="isolated" if product == "m12a" else None,
        basis="isolated yield" if product == "m12a" else None,
        scope="Pd-catalyzed coupling step")]
    unresolved = ["Temperature, time, solvent, workup, and purification are unreported in the main text."]
    if product == "m12a":
        outcomes.append(outcome(
            "o2", "selectivity", [product], ev, "single stereoisomer", "single stereoisomer",
            basis="within the detection limit of 1H NMR", scope="stereoisomeric purity"))
    else:
        unresolved.append("The 99% yield basis for 12b and its intermediate stereoisomer ratio are not reported; the 13b result is recorded separately.")
    return reaction(
        suffix, ev, "specific_experiment", ["other"], f"Page 3 Scheme 3, {product} formation",
        "Pd-catalyzed asymmetric coupling", "allene hydroalkoxylation", parts,
        [step("coupling", ev, [p["participant_id"] for p in parts],
              "Pd-catalyzed coupling of allene 5 and allylic alcohol 11")],
        outcomes,
        HYDRO, (f"Allene 5 and alcohol 11 gave {product[1:]} as a single stereoisomer."
                if product == "m12a" else
                f"Allene 5 and alcohol 11 gave the depicted intermediate {product[1:]} in 99% yield."),
        unresolved, (GS_A,))

bundle["reactions"].append(scheme3_coupling("scheme3_5_11_to_12a", "lig_ss", 1, 2, "m12a", 99))
bundle["reactions"].append(scheme3_coupling("scheme3_5_11_to_12b", "lig_rr", 2, 1, "m12b", 99))

for suffix, start, product, y, selectivity in [
    ("scheme3_12a_to_13a", "m12a", "m13a", 79, ">25:1 selectivity"),
    ("scheme3_12b_to_13b", "m12b", "m13b", 87, "<1:25 selectivity"),
]:
    ev = [I[3], P[2]]
    specs = [(start, "reactant", None, {}),
             ("ru_cat", "catalyst", None, {"mol_percent": ("5 mol %", q(5, "mol%"), start)}),
             (product, "product", "biscyclic product", {})]
    parts = parts_for(specs, ev)
    outs = [outcome("o1", "yield", [product], [I[3]], f"{y}%", q(y, "%"),
                    "biscyclic product", scope=f"single RCM step from {start[1:]}")]
    if product == "m13b":
        outs.append(outcome("o2", "yield", [product], [P[2]], "86% yield over two steps",
                            q(86, "%"), "biscyclic product",
                            scope="overall yield over Pd coupling and RCM from 5 and 11"))
        soid = "o3"
    else:
        soid = "o2"
    outs.append(outcome(soid, "selectivity", [product], [I[3]], selectivity, selectivity,
                        scope="reported diastereoselectivity; isomer order is not explicitly printed beside the ratio"))
    unresolved = ["Temperature, time, solvent, workup, purification, and step-yield basis are unreported."]
    if product == "m13b":
        unresolved.append(
            "The 13b scheme step yield (87%) and prose overall two-step yield (86%) have different scopes and are not merged.")
    bundle["reactions"].append(reaction(
        suffix, ev, "specific_experiment", ["other"], f"Page 3 Scheme 3, {start[1:]} to {product[1:]}",
        "ring-closing metathesis", "ring-closing metathesis", parts,
        [step("rcm", ev, [p["participant_id"] for p in parts],
              "Treatment with second-generation Hoveyda–Grubbs catalyst")],
        outs, RCM, f"{start[1:]} was treated with Ru catalyst to give {product[1:]}.",
        unresolved,
        (GS_A,)))

# Scheme 3 elaboration 13a -> 17.
ev = [I[3], P[2]]
specs = [("m13a", "reactant", None, {}),
         ("oso4", "catalyst", None, {"mol_percent": ("3 mol %", q(3, "mol%"), "m13a")}),
         ("nmo", "reagent", None, {"equivalents": ("1.3 equiv", q(1.3, "equiv"), "m13a")}),
         ("m14", "product", "single diol product", {})]
parts = parts_for(specs, ev)
bundle["reactions"].append(reaction(
    "scheme3_13a_to_14", ev, "specific_experiment", ["other"], "Page 3 Scheme 3, arrow i",
    "Upjohn dihydroxylation", "alkene dihydroxylation", parts,
    [step("dihydroxylation", ev, [p["participant_id"] for p in parts], "OsO4/NMO dihydroxylation")],
    [outcome("o1", "yield", ["m14"], ev, "81%", q(81, "%"), "diol 14", scope="13a to 14"),
     outcome("o2", "selectivity", ["m14"], ev,
             "single compound with complete chemo- and stereoselectivity",
             "single compound; complete chemo- and stereoselectivity",
             scope="chemo- and diastereoselectivity; analytical basis not reported")],
    DIOL, "13a was subjected to OsO4 (3 mol %) and NMO (1.3 equiv).",
    ["Solvent, temperature, time, workup, purification, and yield basis are unreported."], (GS_A,)))

specs = [("m14", "reactant", None, {}),
         ("momcl", "reagent", None, {"equivalents": ("3 equiv", q(3, "equiv"), "m14")}),
         ("dipea", "reagent", None, {"equivalents": ("5 equiv", q(5, "equiv"), "m14")}),
         ("m15", "product", "MOM-protected product", {})]
parts = parts_for(specs, ev)
bundle["reactions"].append(reaction(
    "scheme3_14_to_15", ev, "specific_experiment", ["other"], "Page 3 Scheme 3, arrow ii",
    "MOM protection", "O-alkylation/protection", parts,
    [step("mom_protection", ev, [p["participant_id"] for p in parts],
          "Treatment with MOMCl and N,N-diisopropylethylamine")],
    [outcome("o1", "yield", ["m15"], ev, "89%", q(89, "%"), "protected product 15", scope="14 to 15")],
    [("C", "O", 0, 1, "formed"), ("C", "O", 0, 1, "formed")],
    "Diol 14 was treated with MOMCl (3 equiv) and N,N-diisopropylethylamine (5 equiv).",
    ["Solvent, temperature, time, workup, purification, and yield basis are unreported."], (GS_A,)))

specs = [("m15", "reactant", None, {}), ("na", "reagent", None, {}),
         ("nh3", "solvent", None, {}), ("m15oh", "product", "debenzylated alcohol", {})]
parts = parts_for(specs, ev)
bundle["reactions"].append(reaction(
    "scheme3_15_to_unnumbered_alcohol", ev, "specific_experiment", ["other"],
    "Page 3 Scheme 3, operation iii", "benzyl-group removal", "Na/NH3 debenzylation", parts,
    [step("debenzylation", ev, [p["participant_id"] for p in parts], "Treatment with Na/NH3",
          ("−78 °C", q(-78, "°C")))],
    [outcome("o1", "yield", ["m15oh"], ev, "93%", q(93, "%"), "unnumbered alcohol", scope="debenzylation of 15")],
    [("O", "H", 0, 1, "hydrogen_gain"), ("C", "O", 1, 0, "broken")],
    "Removal of the benzyl group using Na/NH3 at −78 °C gave the unnumbered alcohol.",
    ["Time, workup, purification, and yield basis are unreported; the product is unnumbered in the main text."], (GS_A,)))

specs = [("m15oh", "reactant", None, {}),
         ("pmp_nco", "reagent", None, {"equivalents": ("1.1 equiv", q(1.1, "equiv"), "m15oh")}),
         ("dbu", "reagent", None, {"equivalents": ("0.1 equiv", q(0.1, "equiv"), "m15oh")}),
         ("ch2cl2", "solvent", None, {}), ("m16", "product", "carbamate product", {})]
parts = parts_for(specs, ev)
r = reaction(
    "scheme3_unnumbered_alcohol_to_16", ev, "specific_experiment", ["other"],
    "Page 3 Scheme 3, operation iv", "carbamate formation", "alcohol addition to aryl isocyanate", parts,
    [step("carbamoylation", ev, [p["participant_id"] for p in parts],
          "Treatment with PMP–NCO and DBU in CH2Cl2")],
    [outcome("o1", "yield", ["m16"], ev, "83%", q(83, "%"), "carbamate 16", scope="unnumbered alcohol to 16")],
    [("C", "O", 0, 1, "formed"), ("N", "H", 0, 1, "hydrogen_gain"),
     ("C", "N", 2, 1, "order_decreased"), ("O", "H", 1, 0, "broken")],
    "The unnumbered alcohol was treated with PMP–NCO (1.1 equiv) and DBU (0.1 equiv) in CH2Cl2.",
    ["Time, workup, purification, and yield basis are unreported."], (GS_A,))
r["steps"][0]["temperature"] = f(
    "using the protocol developed in Scheme 2", "rt", ev,
    method="inherited_general_condition",
    derivation="The prose explicitly says the Scheme 2 protocol was used; Scheme 2 states rt for carbamate formation.")
bundle["reactions"].append(r)

specs = [("m16", "reactant", None, {}),
         ("ibx", "reagent", None, {"equivalents": ("Scheme 2 protocol: 6 equiv", q(6, "equiv"), "m16")}),
         ("nahco3", "reagent", None, {"equivalents": ("Scheme 2 protocol: 12 equiv", q(12, "equiv"), "m16")}),
         ("thf", "solvent", None, {}), ("dmso", "solvent", None, {}),
         ("m17", "product", "disaccharide target", {})]
parts = parts_for(specs, ev)
r = reaction(
    "scheme3_16_to_17", ev, "specific_experiment", ["other"], "Page 3 Scheme 3, operation v",
    "IBX-mediated oxidative cyclization", "oxidative intramolecular C–N bond formation", parts,
    [step("oxidative_cyclization", ev, [p["participant_id"] for p in parts],
          "Heat with IBX and NaHCO3 in THF:DMSO in a sealed tube",
          ("90 °C", q(90, "°C")))],
    [outcome("o1", "yield", ["m17"], ev, "75%", q(75, "%"), "disaccharide 17", scope="16 to 17")],
    [("C", "N", 0, 1, "formed"), ("C", "C", 2, 1, "order_decreased")],
    "Compound 16 was heated with IBX in DMSO:THF (10:1) with NaHCO3 in a sealed tube at 90 °C.",
    ["Reaction time, workup, purification, and yield basis are unreported."], (GS_A,))
r["solvent_conditions"]["composition"] = f("DMSO:THF", ["DMSO", "THF"], ev)
r["solvent_conditions"]["mixture_ratio"] = f("DMSO:THF (10:1)", {"DMSO": 10, "THF": 1}, ev)
r["special_conditions"]["reactor"] = f("sealed tube", "sealed tube", ev)
# The loadings are explicitly inherited from the stated Scheme 2 protocol.
for p in r["participants"]:
    if p["molecule_id"] in ("ibx", "nahco3"):
        p["charges"]["equivalents"]["acquisition_method"] = "inherited_general_condition"
        p["charges"]["equivalents"]["derivation"] = "The prose states that the Scheme 2 protocol was used; the loadings are printed in Scheme 2."
bundle["reactions"].append(r)

# Scheme 4 preliminary protected-substrate test and two ligand-controlled runs.
ev = [I[3], P[3]]
specs = [("m18", "reactant", None, {"equivalents": ("1.0–2.0 equiv", {"minimum": 1.0, "maximum": 2.0, "unit": "equiv", "raw_text": "1.0–2.0 equiv"}, "m19")}),
         ("m19", "reactant", None, {"equivalents": ("1.1 equiv", q(1.1, "equiv"), "m19")}),
         ("pd2dba3", "catalyst", None, {"mol_percent": ("2.5 mol %", q(2.5, "mol%"), "m18")}),
         ("lig_ss", "ligand", None, {"mol_percent": ("5 mol %", q(5, "mol%"), "m18")}),
         ("ru_cat", "catalyst", None, {"mol_percent": ("10 mol %", q(10, "mol%"), "m18")}),
         ("m18_adduct", "product", "unassigned trace adduct", {})]
parts = parts_for(specs, ev)
bundle["reactions"].append(reaction(
    "scheme4_18_19_trace", ev, "specific_experiment", ["optimization"],
    "Page 3 Scheme 4 and preliminary-test prose", "protected-substrate coupling test",
    "allene hydroalkoxylation/RCM test", parts,
    [step("coupling_rcm", ev, [p["participant_id"] for p in parts],
          "Pd-catalyzed coupling followed by Ru-catalyzed RCM")],
    [outcome("o1", "trace", ["m18_adduct"], ev, "only trace yield", "trace",
             "unassigned adduct", scope="protected substrate 18 plus allene 19")],
    HYDRO + RCM, "Protected substrate 18 and allene 19 formed an adduct only in trace yield.",
    ["The exact adduct identity, stereochemistry, measurement basis, temperature, time, solvent, workup, and purification are unreported."], (GS_B,)))

def scheme4_run(suffix, ligand, major, major_y, minor, minor_raw, regio_raw):
    ev = [I[3], P[3]]
    specs = [("m20", "reactant", None, {"equivalents": ("1.0 equiv", q(1, "equiv"), "m20")}),
             ("m19", "reactant", None, {"equivalents": ("1.1 equiv", q(1.1, "equiv"), "m20")}),
             ("pd2dba3", "catalyst", None, {"mol_percent": ("2.5 mol %", q(2.5, "mol%"), "m20")}),
             (ligand, "ligand", None, {"mol_percent": ("5 mol %", q(5, "mol%"), "m20")}),
             ("et3n", "additive", None, {"equivalents": ("0.1 equiv", q(0.1, "equiv"), "m20")}),
             ("ru_cat", "catalyst", None, {"mol_percent": ("10 mol %", q(10, "mol%"), "m20")}),
             (major, "product", "major stereoisomer", {}),
             (minor, "product", "minor stereoisomer", {}),
             ("m21regio", "product", "regioisomer mixture", {})]
    parts = parts_for(specs, ev)
    outs = [outcome("o1", "yield", [major], ev, f"{major_y}%", q(major_y, "%"),
                    "desired stereoisomer", scope="overall yield after Pd coupling and RCM"),
            outcome("o2", "trace", [minor], ev, minor_raw, "trace", "minor stereoisomer",
                    scope="stereoisomer amount described qualitatively")]
    if "%" in regio_raw:
        outs.append(outcome("o3", "yield", ["m21regio"], ev, regio_raw,
                            {"number": 13, "unit": "%", "qualifier": "approx", "raw_text": "~13% combined yield"},
                            "regioisomer mixture", scope="combined yield of regioisomers"))
    else:
        outs.append(outcome("o3", "trace", ["m21regio"], ev, regio_raw, "trace",
                            "regioisomer mixture", scope="regioisomer amount described qualitatively"))
    return reaction(
        suffix, ev, "specific_experiment", ["other"], f"Page 3 Scheme 4, {ligand} row",
        "protective-group-free ligand-controlled coupling", "allene hydroalkoxylation followed by RCM",
        parts, [step("coupling", ev, [p["participant_id"] for p in parts],
                     "Pd-catalyzed coupling of triol 20 and allene 19"),
                step("rcm", ev, [p["participant_id"] for p in parts],
                     "Subsequent Ru-catalyzed ring-closing metathesis")],
        outs, HYDRO + RCM,
        "Triol 20 and allene 19 were subjected to Pd/chiral ligand/Et3N and subsequent Ru-catalyzed RCM.",
        ["Temperature, time, solvent, workup, purification, and all yield measurement bases are unreported."], (GS_B,))

bundle["reactions"].append(scheme4_run(
    "scheme4_20_19_ss", "lig_ss", "m21a", 76, "m21b", "trace", "~13% combined yield"))
bundle["reactions"].append(scheme4_run(
    "scheme4_20_19_rr", "lig_rr", "m21b", 87, "m21a", "trace", "trace regioisomer"))

# 21a derivatization and final dihydroxylation.
specs = [("m21a", "reactant", None, {}),
         ("benzoyl_chloride", "reagent", None, {"equivalents": ("2 equiv", q(2, "equiv"), "m21a")}),
         ("pyridine", "solvent", None, {}), ("m22", "product", "bisbenzoate", {})]
parts = parts_for(specs, ev)
bundle["reactions"].append(reaction(
    "scheme4_21a_to_22", ev, "specific_experiment", ["derivatization"],
    "Page 3 Scheme 4, 21a to 22", "bisbenzoylation", "diol benzoyl protection", parts,
    [step("benzoylation", ev, [p["participant_id"] for p in parts],
          "Treatment with benzoyl chloride in pyridine",
          ("rt", "rt"), ("12 h", q(12, "h")))],
    [outcome("o1", "yield", ["m22"], ev, "86%", q(86, "%"), "bisbenzoate 22", scope="21a to 22")],
    [("C", "O", 0, 1, "formed"), ("C", "O", 0, 1, "formed")],
    "21a was treated with benzoyl chloride (2 equiv) in pyridine at rt for 12 h.",
    ["Workup, purification, and yield basis are unreported; detailed NMR analysis is in unavailable SI."], (GS_B,)))

specs = [("m21a", "reactant", None, {}),
         ("oso4", "catalyst", None, {"mol_percent": ("3 mol %", q(3, "mol%"), "m21a")}),
         ("nmo", "reagent", None, {"equivalents": ("1.3 equiv", q(1.3, "equiv"), "m21a")}),
         ("acetone", "solvent", None, {}), ("thf", "solvent", None, {}),
         ("m23", "product", "single-diastereomer tetrol", {})]
parts = parts_for(specs, ev)
r = reaction(
    "scheme4_21a_to_23", ev, "specific_experiment", ["other"],
    "Page 3 Scheme 4, 21a to 23", "dihydroxylation", "alkene dihydroxylation", parts,
    [step("dihydroxylation", ev, [p["participant_id"] for p in parts],
          "OsO4/NMO dihydroxylation in acetone/THF")],
    [outcome("o1", "yield", ["m23"], ev, "90%", q(90, "%"), "tetrol 23", scope="21a to 23"),
     outcome("o2", "selectivity", ["m23"], ev, "single diastereomer", "single diastereomer",
             scope="diastereoselectivity; measurement method not reported")],
    DIOL, "21a was treated with OsO4 (3 mol %) and NMO (1.3 equiv) in acetone/THF (1:1).",
    ["Temperature, time, workup, purification, and yield basis are unreported."], (GS_B,))
r["solvent_conditions"]["composition"] = f("acetone/THF", ["acetone", "THF"], ev)
r["solvent_conditions"]["mixture_ratio"] = f("acetone/THF (1:1)", {"acetone": 1, "THF": 1}, ev)
bundle["reactions"].append(r)

# Non-synthetic observations and limitations retained separately.
bundle["source_context"] = [
    {
        "context_id": f"{DOC}:context:saccharomicin_structures",
        "kind": "structure_context",
        "evidence_ids": [I[1]],
        "note": "Scheme 1 shows saccharomicins A/B and the Sac–Rha–Fuc motif; this is structural context, not a separate experiment.",
    },
    {
        "context_id": f"{DOC}:context:nmr_assignments",
        "kind": "characterization_reference",
        "evidence_ids": [P[2], P[3]],
        "note": "The main text mentions NMR confirmation for 7a, 13a/13b, and 22, but no independent peak table is present on pages 1–4 and the cited SI was not supplied.",
    },
    {
        "context_id": f"{DOC}:context:si_unassociated",
        "kind": "supplementary_material_notice",
        "evidence_ids": [I[3]],
        "note": "Page 3 states that Supporting Information exists. It was not provided or explicitly associated with this session and was not used.",
    },
    {
        "context_id": f"{DOC}:context:references",
        "kind": "non_experimental_page",
        "evidence_ids": [I[4]],
        "note": "Page 4 contains references and explanatory footnotes only; no current-work reaction arrow or experimental table occurs on that page.",
    },
]

inventory = {
    "document_id": DOC,
    "pages_reviewed": [1, 2, 3, 4],
    "machine_extraction": True,
    "semantic_accuracy_verified": False,
    "record_count": len(bundle["reactions"]),
    "molecule_count": len(bundle["molecules"]),
    "records": [{
        "reaction_id": r["reaction_id"],
        "record_type": r["record_type"],
        "location": r["location"]["raw_value"],
        "evidence_ids": r["evidence_ids"],
    } for r in bundle["reactions"]],
    "notes": [
        "Graphical abstract duplicates Scheme 1 concepts and was not counted as an additional experiment.",
        "All records remain needs_review; image inspection does not establish independent chemical correctness.",
        "No SI was associated; unavailable procedures and characterization were not inferred.",
    ],
}
st.write(ROOT / "runs/current/raw_model/inventory.json", inventory)
st.write(ROOT / "runs/current/submitted.json", bundle)
print(json.dumps({"molecules": len(bundle["molecules"]),
                  "reactions": len(bundle["reactions"])}, ensure_ascii=False))
