#!/usr/bin/env python3
"""Reproducible visual transcription for Org. Lett. 2020, 22, 395–399.

This generator records only facts visible on PDF pages 1–5 or in their extracted
text.  It deliberately leaves dense structures, missing SI procedures, and
unprinted measurement bases unresolved.
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
DOC = "doc_d394fafa02fef0e2"
E = {i: f"{DOC}:page{i}:image" for i in range(1, 6)}
T = {i: f"{DOC}:page{i}:text" for i in range(1, 6)}


def f(raw=None, normalized=None, ev=(), method="original_extraction", derivation=None,
      status=None):
    return st.fact(raw, normalized, ev, method=method, derivation=derivation, status=status)


def inferred(normalized, ev, derivation):
    return f(None, normalized, ev, method="inference", derivation=derivation)


def derived(normalized, ev, derivation):
    return f(None, normalized, ev, method="structure_derivation", derivation=derivation)


def missing(status="unreported"):
    return f(method="not_applicable" if status == "not_applicable" else "original_extraction",
             status=status)


def q(number, unit, qualifier=None, raw_text=None):
    out = {"number": number, "unit": unit}
    if qualifier:
        out["qualifier"] = qualifier
    if raw_text:
        out["raw_text"] = raw_text
    return out


def mid(label):
    safe = label.replace(" ", "_").replace("/", "_").replace("·", "_").replace("+", "plus")
    return f"{DOC}:mol_{safe}"


bundle = st.empty(DATASET["documents"])
molecules = {}


def set_statuses(m, ev, chirality="unresolved", geometry="unresolved", representation="not_encoded",
                 reason="Source drawing was not safely encoded as an isomeric molecular graph."):
    m["stereochemistry"]["chirality_status"] = inferred(chirality, ev, reason)
    m["stereochemistry"]["geometry_status"] = inferred(geometry, ev, reason)
    m["stereochemistry"]["representation_status"] = inferred(representation, ev, reason)


def add_unresolved(label, ev, description, *, chirality="unresolved", geometry="unresolved",
                   representation="not_encoded", atom_centers=None, double_bonds=None, axial=None,
                   mixture=None, assignment=None, aliases=None):
    identity = mid(label)
    if identity in molecules:
        occ = molecules[identity]["occurrences"][0]
        for eid in ev:
            if eid not in occ["evidence_ids"]:
                occ["evidence_ids"].append(eid)
            if eid not in occ["compound_label"]["evidence_ids"]:
                occ["compound_label"]["evidence_ids"].append(eid)
        return identity
    m = st.new_molecule(identity, DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f(label, label, ev)
    m["names"]["original"] = f(label, label, ev)
    if aliases:
        m["names"]["aliases"] = f(aliases, aliases, ev)
    m["identity_note"] = derived(description, ev, "Description derived from the source label and drawing; exact graph remains pending review.")
    set_statuses(m, ev, chirality, geometry, representation)
    if atom_centers is not None:
        m["stereochemistry"]["atom_centers"] = f(atom_centers, atom_centers, ev)
    if double_bonds is not None:
        m["stereochemistry"]["double_bonds"] = f(double_bonds, double_bonds, ev)
    if axial is not None:
        m["stereochemistry"]["axial"] = f(axial, axial, ev)
    if mixture is not None:
        m["stereochemistry"]["mixture"] = f(mixture, mixture, ev)
    if assignment is not None:
        m["stereochemistry"]["assignment"] = f(assignment, assignment, ev)
    molecules[identity] = m
    return identity


stereo_diagnostics = {}


def add_smiles(label, smiles, ev, description, *, chirality="achiral", geometry="not_applicable",
               representation="encoded", axial=None):
    identity = mid(label)
    diag = st.stereo_inspect(smiles)
    stereo_diagnostics[label] = diag
    m = st.new_molecule(identity, DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f(label, label, ev)
    m["names"]["original"] = f(label, label, ev)
    m["representation_type"] = "smiles"
    m["smiles"] = derived(smiles, ev, "Candidate connectivity transcribed from the source drawing/name; parsing is not source verification.")
    m["identity_note"] = derived(description, ev, "Identity description derived from the source drawing/name.")
    set_statuses(m, ev, chirality, geometry, representation,
                 "Stereo status inferred from the source depiction plus graph diagnostics; toolkit output is not source verification.")
    iso = diag["canonical_isomeric_smiles"]
    m["stereochemistry"]["isomeric_smiles"] = derived(iso, ev, "Canonical isomeric SMILES from bundled standard stereo diagnostics.")
    m["stereochemistry"]["identity_key"] = derived(diag["stereo_graph_key"], ev, "Stereo graph key from bundled standard stereo diagnostics.")
    m["stereochemistry"]["atom_centers"] = derived(diag["atom_centers"], ev, "RDKit graph diagnostics; not a source assignment.")
    m["stereochemistry"]["double_bonds"] = derived(diag["double_bonds"], ev, "RDKit graph diagnostics; source E/Z must be reviewed separately.")
    if axial is not None:
        m["stereochemistry"]["axial"] = f(axial, axial, ev)
    molecules[identity] = m
    return identity


def add_template(label, ev, fixed, variables, reason, *, chirality="unresolved", geometry="unresolved"):
    identity = mid(label)
    m = st.new_molecule(identity, DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f(label, label, ev)
    m["names"]["original"] = f(label, label, ev)
    m["representation_type"] = "generic_template"
    m["generic_template"] = {
        "fixed_structure": f(fixed, fixed, ev),
        "variables": [{
            "label": v[0],
            "attachment_points": f(v[1], v[1], ev),
            "definition": f(v[2], v[2], ev),
            "allowed_values": f(v[3], v[3], ev) if v[3] is not None else missing(),
            "constraints": missing(),
        } for v in variables],
        "instance_mapping": [],
        "reason": f(reason, reason, ev),
    }
    m["identity_note"] = f(reason, reason, ev)
    set_statuses(m, ev, chirality, geometry, "not_encoded", reason)
    if chirality in {"specified", "partially_specified"}:
        m["stereochemistry"]["assignment"] = f(
            "wedge/dash stereochemistry in the generic source drawing",
            {"source": "generic source drawing", "exact_R_S": "not printed"}, ev)
    molecules[identity] = m
    return identity


# Named reagents and solvents. Xylene is deliberately not collapsed to one isomer.
M_TOL = add_smiles("toluene", "Cc1ccccc1", [E[2], E[4]], "Toluene solvent named in Table 1/Figure 2.")
M_MEOH = add_smiles("methanol", "CO", [E[4]], "Methanol solvent named in Figure 2c.")
M_CS2CO3 = add_smiles("Cs2CO3", "[Cs+].[Cs+].[O-]C([O-])=O", [E[2]], "Cesium carbonate base.")
M_K2CO3 = add_smiles("K2CO3", "[K+].[K+].[O-]C([O-])=O", [E[2]], "Potassium carbonate base.")
M_K3PO4 = add_smiles("K3PO4", "[K+].[K+].[K+].[O-]P([O-])([O-])=O", [E[2], E[3], E[4]], "Tripotassium phosphate base.")
M_D2O = add_smiles("D2O", "[2H]O[2H]", [E[4]], "Deuterium oxide additive in Figure 2b.")
M_XYL = add_template("xylene", [E[2], E[3], E[4]], "dimethylbenzene solvent",
                     [("isomer/composition", "benzene substitution positions", "unspecified xylene isomer or mixture", ["o-xylene", "m-xylene", "p-xylene", "mixture"])],
                     "The source says only xylene; a unique constitutional isomer or mixture composition is not reported.",
                     chirality="not_applicable", geometry="not_applicable")

# Generic source drawings used by the three tabular schemes.
G1 = add_template("generic_1", [E[3]], "3-alkenyl-oxindole scaffold with Ar, EWG1 and Pg variables",
                  [("Ar", "oxindole arene", "substitution pattern shown in Table 2", None),
                   ("EWG1", "exocyclic alkene substituent", "electron-withdrawing group", None),
                   ("Pg", "oxindole nitrogen", "protecting group", None)],
                  "Table 2 uses variable labels rather than one concrete molecule.", geometry="unspecified")
G2 = add_template("generic_2", [E[3]], "electron-deficient allene scaffold with EWG2, EWG3 and R",
                  [("EWG2", "one allene terminus", "electron-withdrawing group", None),
                   ("EWG3", "opposite allene terminus", "electron-withdrawing group", None),
                   ("R", "allene terminus", "substituent", None)],
                  "Tables 2 and 3 use variable labels for the allenoate.", chirality="unspecified", geometry="not_applicable")
G34 = add_template("generic_3_and_4", [E[3]], "multifunctionalized 4H-pyran fused/spiro oxindole product scaffold",
                   [("Ar", "oxindole arene", "scope substituent", None),
                    ("EWG1", "pyran substituent", "electron-withdrawing group", None),
                    ("EWG2", "pyran substituent", "electron-withdrawing group", None),
                    ("EWG3", "pyran substituent", "electron-withdrawing group", None),
                    ("R", "stereogenic carbon substituent", "scope substituent", None)],
                   "Table 2 generic product contains variable groups and drawn relative stereochemistry.", chirality="specified", geometry="not_applicable")
G5 = add_template("generic_5", [E[3]], "ketone-activated 3-alkenyl-oxindole scaffold with Ar and Pg",
                  [("Ar", "ketone aryl group", "aryl substituent", None),
                   ("Pg", "oxindole nitrogen", "protecting group", None)],
                  "Table 3 uses variable labels rather than one concrete molecule.", geometry="unspecified")
G6 = add_template("generic_6", [E[3]], "ketone-derived multifunctionalized 4H-pyran product scaffold",
                  [("Ar", "ketone aryl group", "aryl substituent", None),
                   ("Pg", "oxindole nitrogen", "protecting group", None),
                   ("R", "stereogenic carbon substituent", "allenoate-derived group", None)],
                  "Table 3 generic product contains variables and drawn relative stereochemistry.", chirality="specified", geometry="not_applicable")

# Main reactants. 2a has a candidate connectivity; its unassigned allene axis is retained explicitly.
M_1A = add_unresolved("1a", [E[2], E[3]], "N-Boc 3-alkenyl-oxindole 1a bearing the source-drawn CO2Bn substituent.", geometry="unspecified",
                      double_bonds={"exocyclic_alkene": "drawn without an E/Z descriptor on pages 2–3"})
M_2A = add_smiles("2a", "CCOC(=O)C(C)=C=CC(=O)OC(C)(C)C", [E[2], E[3], E[4]],
                  "Electron-deficient allenoate 2a; candidate connectivity only.", chirality="unspecified",
                  geometry="not_applicable", representation="partial",
                  axial={"axis": "C=C=C", "configuration": "unspecified", "source": "no Ra/Sa or P/M descriptor printed"})

# Catalysts shown on pages 2 and 4. Exact isomeric graphs are retained as unresolved drawings.
cat_ids = {}
for label in [f"P{i}" for i in range(1, 11)] + ["P10-1", "P10-2"]:
    ev = [E[2]] if label in {f"P{i}" for i in range(1, 11)} else [E[4]]
    if label == "P10":
        ev = [E[2], E[3], E[4]]
    cat_ids[label] = add_unresolved(label, ev,
        f"Source-drawn bifunctional phosphonium salt catalyst {label}; counterion and substituents are retained only in the drawing.",
        chirality="specified", geometry="not_applicable", representation="not_encoded",
        assignment={"source": "described/drawn as a chiral amino-acid-derived phosphonium salt", "absolute_descriptor": "not printed on pages 1–5"})


def product_stereo(label, ev):
    assignment = {
        "source": "wedge/dash product drawing",
        "configuration": "relative/absolute stereochemical depiction present; R/S descriptors not printed in the scope panel",
    }
    if label == "3l":
        assignment = {
            "source": "single-crystal X-ray diffraction",
            "compound": "3l",
            "ccdc": "1947826",
            "configuration": "source drawing specified; exact R/S text not printed on pages 1–5",
        }
    return assignment


# Individual Table 2 molecules.
scope_3 = {
    "3a": (98, 93, "EWG1=CO2Bn"), "3b": (98, 99, "EWG1=CO2Me"),
    "3c": (96, 93, "EWG1=CO2Et"), "3d": (96, 93, "EWG1=CO2iPr"),
    "3e": (90, 93, "EWG1=CO2tBu"), "3f": (92, 86, "EWG1=Fmoc"),
    "3g": (88, 84, "N-Cbz product"), "3h": (96, 92, "R=Cl"),
    "3i": (91, 91, "R=Br"), "3j": (93, 90, "R=F"),
    "3k": (92, 87, "R=Cl"), "3l": (90, 92, "R=Br; X-ray reference"),
    "3m": (70, 85, "R=NO2"), "3n": (89, 91, "R=Me"),
    "3o": (96, 93, "R=OMe"), "3p": (91, 90, "R=Cl"),
    "3q": (93, 95, "R=Br"), "3r": (95, 90, "R=OMe"),
    "3s": (72, 93, "R=F"), "3t": (91, 96, "R=Cl"),
    "3u": (92, 94, "R=Br"), "3v": (99, 91, "R=Me"),
    "3w": (92, 85, "multi-substituted arene as drawn"),
    "3x": (92, 88, "multi-methyl arene as drawn"),
    "3y": (85, 90, "difluoro arene as drawn"),
}
substrate_1 = {}
product_3 = {}
for label, (_, _, desc) in scope_3.items():
    slabel = "1" + label[1:]
    substrate_1[slabel] = add_unresolved(slabel, [E[3]], f"Table 2 substrate {slabel} corresponding to product {label}; {desc}.", geometry="unspecified")
    product_3[label] = add_unresolved(label, [E[3]], f"Table 2 multifunctionalized 4H-pyran product {label}; {desc}.",
        chirality="specified", geometry="not_applicable", representation="not_encoded",
        atom_centers={"source": "bold/wedged bonds in product panel", "exact_R_S": "not printed"},
        assignment=product_stereo(label, [E[3]]))
# Same identities also occur in optimization/control regions.
add_unresolved("1a", [E[2]], "N-Boc 3-alkenyl-oxindole 1a bearing the source-drawn CO2Bn substituent.", geometry="unspecified")
add_unresolved("3a", [E[2]], "Optimization product 3a, a multifunctionalized 4H-pyran.", chirality="specified", geometry="not_applicable", assignment=product_stereo("3a", [E[2]]))
add_unresolved("1b", [E[4]], "3-Alkenyl-oxindole 1b used in mechanistic/control experiments.", geometry="unspecified")
add_unresolved("3b", [E[4]], "4H-Pyran product 3b in mechanistic/control experiments.", chirality="specified", geometry="not_applicable", assignment=product_stereo("3b", [E[4]]))

scope_4 = {
    "4a": (98, 82, "R=Me"), "4b": (97, 85, "R=Et"),
    "4c": (88, 80, "CO2tBu variant"), "4d": (96, 93, "R=nBu"),
    "4e": (97, 90, "R=Bn"),
}
allene_4 = {}
product_4 = {}
for label, (_, _, desc) in scope_4.items():
    allene_4[label] = add_unresolved(f"allenoate_for_{label}", [E[3]], f"Table 2 allenoate used for {label}; {desc}.",
        chirality="unspecified", geometry="not_applicable", axial={"axis": "C=C=C", "configuration": "unspecified"})
    product_4[label] = add_unresolved(label, [E[3]], f"Table 2 4H-pyran product {label}; {desc}.",
        chirality="specified", geometry="not_applicable", atom_centers={"source": "wedged product drawing", "exact_R_S": "not printed"},
        assignment=product_stereo(label, [E[3]]))

# Individual Table 3 molecules and values.
scope_6 = {
    "6a": (90, "19:1", 98, "Ar=Ph, Pg=Boc, R=Me"),
    "6b": (96, ">20:1", 97, "ortho-OMe aryl"),
    "6c": (92, ">20:1", 95, "ortho-F aryl"),
    "6d": (96, "16:1", 96, "meta-Cl aryl"),
    "6e": (93, "12:1", 95, "meta-Br aryl"),
    "6f": (94, ">20:1", 98, "meta-OMe aryl"),
    "6g": (95, ">20:1", 95, "para-Me aryl"),
    "6h": (97, ">20:1", 96, "para-OMe aryl"),
    "6i": (75, ">20:1", 93, "para-CN aryl"),
    "6j": (98, ">20:1", 97, "naphthyl ketone"),
    "6k": (93, "9:1", 88, "Pg=Cbz"),
    "6l": (92, ">20:1", ">99.9", "R=nBu allenoate variant"),
}
substrate_5 = {}
product_6 = {}
for label, (_, _, _, desc) in scope_6.items():
    if label != "6l":
        slabel = "5" + label[1:]
        substrate_5[slabel] = add_unresolved(slabel, [E[3]], f"Table 3 ketone-activated 3-alkenyl-oxindole {slabel}; {desc}.", geometry="unspecified")
    product_6[label] = add_unresolved(label, [E[3]], f"Table 3 4H-pyran product {label}; {desc}.",
        chirality="specified", geometry="not_applicable", atom_centers={"source": "wedged product drawing", "exact_R_S": "not printed"},
        assignment=product_stereo(label, [E[3]]))
allene_6l = add_unresolved("allenoate_for_6l", [E[3]], "n-Butyl-substituted allenoate used for 6l.",
    chirality="unspecified", geometry="not_applicable", axial={"axis": "C=C=C", "configuration": "unspecified"})

# Mechanistic-labeling product and incomplete prose-only products.
M_D3B = add_unresolved("d-3b", [E[4]], "Deuterium-labeled 3b with D at the alpha-prime position, 62% D incorporation.",
    chirality="specified", geometry="not_applicable", assignment=product_stereo("d-3b", [E[4]]))
M_7 = add_unresolved("7", [E[3]], "Optically enriched derivatization product 7; structure is only in unavailable Figure S2b.",
    chirality="partially_specified", geometry="unresolved", assignment={"source": "described as optically enriched", "configuration": "unreported on pages 1–5"})
M_8 = add_unresolved("8", [E[3]], "Optically enriched derivatization product 8; structure is only in unavailable Figure S2b.",
    chirality="partially_specified", geometry="unresolved", assignment={"source": "described as optically enriched", "configuration": "unreported on pages 1–5"})

# Literature-context and mechanism-only identities.
prior_mols = {}
for label in ["prior_Jorgensen_inputs", "prior_Jorgensen_products", "prior_Maruoka_inputs", "prior_Maruoka_products"]:
    prior_mols[label] = add_unresolved(label, [E[2]], f"Generic prior-literature structure set labeled in Figure 1a: {label}.")
mech_mols = {}
for label in ["intermediate_A", "intermediate_B", "intermediate_C", "intermediate_D", "intermediate_E", "phosphonium_ion_pair", "HBr"]:
    mech_mols[label] = add_unresolved(label, [E[4]], f"Proposed mechanism species {label} in Figure 2a; not directly observed on pages 1–5.")


def part(pid, molecule, role, ev, *, amount=None, equivalents=None, mol_percent=None, volume=None,
         reference=None, product_role=None, stage="reaction"):
    p = st.new_participant(pid, molecule, role, ev, stage=stage)
    if product_role is not None:
        p["product_role"] = f(product_role, product_role, ev)
    else:
        p["product_role"] = missing("not_applicable")
    if amount is not None:
        raw = f"{amount:g} mmol"
        p["charges"]["amount"] = f(raw, q(amount, "mmol", raw_text=raw), ev)
    if equivalents is not None:
        raw = f"{equivalents:g} equiv"
        p["charges"]["equivalents"] = f(raw, q(equivalents, "equiv", raw_text=raw), ev)
    if mol_percent is not None:
        raw = f"{mol_percent:g} mol%"
        p["charges"]["mol_percent"] = f(raw, q(mol_percent, "mol%", raw_text=raw), ev)
    if volume is not None:
        raw = f"{volume:g} mL"
        p["charges"]["volume"] = f(raw, q(volume, "mL", raw_text=raw), ev)
    if reference is not None:
        p["equivalence_reference"] = inferred(reference, ev, "Loading/equivalent denominator inferred as the limiting oxindole substrate from the displayed scheme and charges.")
    return p


def class_annulation(ev, deuterium=False):
    changes = [
        {"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}],
         "before_order": 0, "after_order": 1, "change_type": "formed"},
        {"change_id": "b2", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "O"}],
         "before_order": 0, "after_order": 1, "change_type": "formed"},
    ]
    labels = [
        {"category_number": 2, "code": "C-C", "supporting_change_ids": ["b1"]},
        {"category_number": 3, "code": "C-O", "supporting_change_ids": ["b2"]},
    ]
    limitations = [
        "Partial source-drawing classification: the two ring-closing connections are represented abstractly without atom maps.",
        "Allene/pyran bond-order changes and atom correspondence remain unresolved; no complete mapping is claimed.",
    ]
    if deuterium:
        changes += [
            {"change_id": "b3", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "H"}],
             "before_order": 1, "after_order": 0, "change_type": "hydrogen_loss", "hydrogen_count": 1},
            {"change_id": "b4", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "H"}],
             "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain", "hydrogen_count": 1},
        ]
        labels.append({"category_number": 1, "code": "C-H", "supporting_change_ids": ["b4"]})
        limitations.append("The C-H label represents source-observed H/D exchange; D is hydrogen isotope 2, and exact atom mapping is unresolved.")
    return {
        "standard_id": st.STANDARD,
        "mapped_reaction_smiles": missing(),
        "bond_changes": [derived(c, ev, "Bond change inferred from the source-drawn annulation/product connectivity.") for c in changes],
        "labels": derived(labels, ev, "Taxonomy labels computed from the partial formed-bond list."),
        "coverage": inferred("partial", ev, "Only confidently source-drawn new ring connections are classified."),
        "limitations": limitations,
        "model_votes": [],
    }


def class_empty(ev, coverage="unresolved", reason="No atom mapping or reliable bond-difference transcription is available."):
    labels = [] if coverage == "not_applicable" else None
    return {
        "standard_id": st.STANDARD,
        "mapped_reaction_smiles": missing(),
        "bond_changes": [],
        "labels": inferred(labels, ev, reason) if labels is not None else missing(),
        "coverage": inferred(coverage, ev, reason),
        "limitations": [reason],
        "model_votes": [],
    }


def base_reaction(suffix, ev, record_type, categories, location, *, general=(), classification=None,
                  reaction_name="formal [4 + 2] annulation"):
    rid = f"{DOC}:{suffix}"
    r = st.new_reaction(rid, DOC, ev)
    r["record_type"] = record_type
    r["experiment_categories"] = f(location, categories, ev)
    r["location"] = f(location, {"source": location}, ev)
    r["general_scheme_ids"] = list(general)
    r["reaction_name"] = f(reaction_name, reaction_name, ev) if reaction_name else missing()
    r["transformation_type"] = inferred(reaction_name, ev, "Normalized from the paper's reaction name and source drawing.") if reaction_name else missing()
    r["classification"] = classification or class_empty(ev)
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None, "checked_fields": [],
                   "unresolved": []}
    r["completeness"] = "partial"
    return r


def set_solvent(r, molecule, ev, *, volume=None, raw_name=None):
    name = raw_name or molecules[molecule]["names"]["original"]["normalized_value"]
    r["solvent_conditions"]["composition"] = f(name, [molecule], ev)
    if volume is not None:
        raw = f"{volume:g} mL"
        r["solvent_conditions"]["total_volume"] = f(raw, q(volume, "mL", raw_text=raw), ev)


def add_step(r, ev, pids, raw, *, temperature="room temperature", hours=8):
    s = st.new_step(f"{r['reaction_id']}:reaction", ev)
    s["participant_ids"] = list(pids)
    s["operation"] = f(raw, raw, ev)
    s["temperature"] = f(temperature, "r.t." if temperature == "room temperature" else temperature, ev)
    if hours is not None:
        raw_h = f"{hours:g} h"
        s["time"] = f(raw_h, q(hours, "h", raw_text=raw_h), ev)
    r["steps"] = [s]


def outcome(r, oid, kind, targets, ev, raw, normalized, *, role=None, yield_type=None,
            basis=None, scope=None):
    o = st.new_outcome(f"{r['reaction_id']}:{oid}", kind, targets, ev)
    o["value"] = f(raw, normalized, ev)
    if role is not None:
        o["product_role"] = f(role, role, ev)
    else:
        o["product_role"] = missing("not_applicable")
    if yield_type is not None:
        o["yield_type"] = f(yield_type, yield_type, ev)
    if basis is not None:
        o["measurement_basis"] = f(basis, basis, ev)
    if scope is not None:
        o["scope"] = f(scope, scope, ev)
    r["outcomes"].append(o)


reactions = []

# Generic schemes corresponding to the three tabular arrows.
SCHEME_T1 = f"{DOC}:table1_general_scheme"
r = base_reaction("table1_general_scheme", [E[2]], "generic_scheme", ["optimization"], "Page 2, Table 1 general arrow",
                  classification=class_annulation([E[2]]))
r["participants"] = [part("t1_g_1a", M_1A, "reactant", [E[2]]), part("t1_g_2a", M_2A, "reactant", [E[2]]),
                     part("t1_g_3a", product_3["3a"], "product", [E[2]], product_role="main product")]
r["review"]["unresolved"] = ["Generic optimization arrow is not an additional experiment; catalyst/base/solvent vary by entry."]
reactions.append(r)

SCHEME_T2 = f"{DOC}:table2_general_scheme"
r = base_reaction("table2_general_scheme", [E[3]], "generic_scheme", ["substrate_scope"], "Page 3, Table 2 general arrow",
                  classification=class_annulation([E[3]]))
r["participants"] = [part("t2_g_1", G1, "reactant", [E[3]]), part("t2_g_2", G2, "reactant", [E[3]]),
                     part("t2_g_p", G34, "product", [E[3]], product_role="generic products 3 and 4")]
r["review"]["unresolved"] = ["Generic variables are not a concrete experimental identity."]
reactions.append(r)

SCHEME_T3 = f"{DOC}:table3_general_scheme"
r = base_reaction("table3_general_scheme", [E[3]], "generic_scheme", ["substrate_scope"], "Page 3, Table 3 general arrow",
                  classification=class_annulation([E[3]]))
r["participants"] = [part("t3_g_5", G5, "reactant", [E[3]]), part("t3_g_2", G2, "reactant", [E[3]]),
                     part("t3_g_6", G6, "product", [E[3]], product_role="generic product 6")]
r["review"]["unresolved"] = ["Generic variables are not a concrete experimental identity."]
reactions.append(r)

# Table 1: every optimization row is a separate experiment.
table1 = [
    (1, "P1", M_TOL, M_CS2CO3, 76, 15, 20), (2, "P2", M_TOL, M_CS2CO3, 91, 5, 20),
    (3, "P3", M_TOL, M_CS2CO3, 90, 3, 20), (4, "P4", M_TOL, M_CS2CO3, 92, 7, 20),
    (5, "P5", M_TOL, M_CS2CO3, 97, 50, 20), (6, "P6", M_TOL, M_CS2CO3, 93, 25, 20),
    (7, "P7", M_TOL, M_CS2CO3, 96, 14, 20), (8, "P8", M_TOL, M_CS2CO3, 95, 45, 20),
    (9, "P9", M_TOL, M_CS2CO3, 98, 73, 20), (10, "P10", M_TOL, M_CS2CO3, 98, 81, 20),
    (11, "P10", M_TOL, M_K2CO3, 89, 83, 20), (12, "P10", M_TOL, M_K3PO4, 98, 91, 20),
    (13, "P10", M_XYL, M_K3PO4, 98, 93, 20), (14, "P10", M_XYL, M_K3PO4, 98, 93, 10),
    (15, "P10", M_XYL, None, None, None, 10),
]
table1_ids = []
for entry, catalyst, solvent, base, yld, ee, loading in table1:
    ev = [E[2]]
    r = base_reaction(f"table1_entry_{entry:02d}", ev, "specific_experiment", ["optimization"],
                      f"Page 2, Table 1, entry {entry}", general=[SCHEME_T1],
                      classification=class_annulation(ev) if yld is not None else class_empty(ev, "not_applicable", "No reaction was observed; intended bond formation is not classified."))
    ps = [part(f"e{entry}_1a", M_1A, "reactant", ev, amount=0.1),
          part(f"e{entry}_2a", M_2A, "reactant", ev, amount=0.12),
          part(f"e{entry}_cat", cat_ids[catalyst], "catalyst", ev, amount=0.02, mol_percent=loading, reference=M_1A),
          part(f"e{entry}_solv", solvent, "solvent", ev, volume=2.0)]
    if base is not None:
        ps.append(part(f"e{entry}_base", base, "reagent", ev, amount=0.2, equivalents=2.0, reference=M_1A))
    if yld is not None:
        ps.append(part(f"e{entry}_3a", product_3["3a"], "product", ev, product_role="main product"))
    r["participants"] = ps
    r["reaction_scale"] = f("1a (0.1 mmol)", {"reference_molecule_id": M_1A, "amount": q(0.1, "mmol")}, ev)
    set_solvent(r, solvent, ev, volume=2.0)
    raw_proc = "1a (0.1 mmol), 2a (0.12 mmol), catalyst (0.02 mmol), base (0.2 mmol) in solvent (2.0 mL) at room temperature for 8 h"
    add_step(r, ev, [p["participant_id"] for p in ps], raw_proc)
    r["procedure"]["original_text"] = f(raw_proc, raw_proc, ev)
    if yld is not None:
        outcome(r, "yield", "yield", [product_3["3a"]], ev, f"{yld}%", q(yld, "%", raw_text=f"{yld}%"),
                role="main product", yield_type="isolated yield", basis="isolated")
        outcome(r, "ee", "selectivity", [product_3["3a"]], ev, f"{ee}% ee", q(ee, "%", raw_text=f"{ee}% ee"),
                role="main product", basis="HPLC on a chiral stationary phase", scope={"metric": "ee", "major_configuration": "not named"})
        outcome(r, "dr", "selectivity", [product_3["3a"]], ev, ">20:1 dr", {"ratio": ">20:1", "metric": "dr"},
                role="main product", basis="1H NMR of crude product", scope={"metric": "dr", "major_diastereomer": "not named"})
    else:
        outcome(r, "no_reaction", "no_reaction", [], ev, "n.r.", "no reaction", scope={"detection_method": "unreported"})
        r["notes"].append(f("ee cell: –", "ee not applicable because no reaction was reported", ev))
    if entry in (14, 15):
        r["notes"].append(f("footnote a: catalyst (0.02 mmol); footnote d: catalyst loading was 10 mol %",
                            "The printed 0.02 mmol and 10 mol% are numerically inconsistent for 0.1 mmol 1a; both are preserved.", ev))
        r["review"]["unresolved"].append("Catalyst charge conflict: general footnote says 0.02 mmol, entry footnote says 10 mol%; no silent correction to 0.01 mmol was made.")
    r["review"]["unresolved"] += ["No SI was supplied; addition order, atmosphere, workup and purification are unreported in the provided pages."]
    reactions.append(r)
    table1_ids.append(r["reaction_id"])


def add_scope_reaction(suffix, substrate, allene, product, yld, ee, dr, scheme, table, ev, desc):
    r = base_reaction(suffix, ev, "specific_experiment", ["substrate_scope"], f"Page 3, {table}, product {suffix.split('_')[-1]}",
                      general=[scheme], classification=class_annulation(ev))
    ps = [part(f"{suffix}_sub", substrate, "reactant", ev, amount=0.2),
          part(f"{suffix}_allene", allene, "reactant", ev, amount=0.24),
          part(f"{suffix}_cat", cat_ids["P10"], "catalyst", ev, amount=0.02, mol_percent=10, reference=substrate),
          part(f"{suffix}_base", M_K3PO4, "reagent", ev, amount=0.4, equivalents=2.0, reference=substrate),
          part(f"{suffix}_solv", M_XYL, "solvent", ev, volume=4.0),
          part(f"{suffix}_prod", product, "product", ev, product_role="main product")]
    r["participants"] = ps
    r["reaction_scale"] = f("substrate (0.2 mmol)", {"reference_molecule_id": substrate, "amount": q(0.2, "mmol")}, ev)
    set_solvent(r, M_XYL, ev, volume=4.0)
    raw_proc = "substrate (0.2 mmol), allenoate (0.24 mmol), P10 (0.02 mmol), K3PO4 (0.4 mmol) in xylene (4.0 mL) at room temperature for 8 h"
    add_step(r, ev, [p["participant_id"] for p in ps], raw_proc)
    r["procedure"]["original_text"] = f(raw_proc, raw_proc, ev)
    outcome(r, "yield", "yield", [product], ev, f"{yld}%", q(yld, "%", raw_text=f"{yld}%"),
            role="main product", yield_type="isolated yield", basis="isolated")
    ee_norm = q(float(ee[1:]) if isinstance(ee, str) and ee.startswith(">") else ee, "%",
                qualifier="gt" if isinstance(ee, str) and ee.startswith(">") else None,
                raw_text=f"{ee}% ee")
    outcome(r, "ee", "selectivity", [product], ev, f"{ee}% ee", ee_norm,
            role="main product", basis="HPLC on a chiral stationary phase", scope={"metric": "ee", "major_configuration": "not named"})
    outcome(r, "dr", "selectivity", [product], ev, f"{dr} dr", {"ratio": dr, "metric": "dr"},
            role="main product", basis="HPLC analysis of crude product", scope={"metric": "dr", "major_diastereomer": "not named"})
    r["notes"].append(f(desc, desc, ev))
    r["review"]["unresolved"] = ["Dense source structures and drawn stereocenters are not encoded as exact isomeric SMILES; exact R/S descriptors remain pending review.",
                                      "No SI was supplied; addition order, atmosphere, workup and purification are unreported in the provided pages."]
    reactions.append(r)
    return r["reaction_id"]


table2_ids = []
for label, (yld, ee, desc) in scope_3.items():
    slabel = "1" + label[1:]
    table2_ids.append(add_scope_reaction(f"table2_{label}", substrate_1[slabel], M_2A, product_3[label], yld, ee, ">20:1",
                                         SCHEME_T2, "Table 2", [E[3]], desc))
for label, (yld, ee, desc) in scope_4.items():
    table2_ids.append(add_scope_reaction(f"table2_{label}", M_1A, allene_4[label], product_4[label], yld, ee, ">20:1",
                                         SCHEME_T2, "Table 2", [E[3]], desc))

table3_ids = []
for label, (yld, dr, ee, desc) in scope_6.items():
    if label == "6l":
        substrate, allene = substrate_5["5a"], allene_6l
    else:
        substrate, allene = substrate_5["5" + label[1:]], M_2A
    table3_ids.append(add_scope_reaction(f"table3_{label}", substrate, allene, product_6[label], yld, ee, dr,
                                         SCHEME_T3, "Table 3", [E[3]], desc))

# Prose-only scale-up and derivatization experiments referred to unavailable Figure S2.
prose_ids = []
for label in ["3b", "3l"]:
    substrate = substrate_1["1" + label[1:]]
    r = base_reaction(f"scale_up_{label}", [E[3], T[3]], "specific_experiment", ["scale_up"],
                      f"Page 3 prose referring to Figure S2a, gram-scale {label}", classification=class_annulation([E[3]]))
    r["participants"] = [part(f"scale_{label}_sub", substrate, "reactant", [E[3]]),
                         part(f"scale_{label}_allene", M_2A, "reactant", [E[3]]),
                         part(f"scale_{label}_prod", product_3[label], "product", [E[3]], product_role="main product")]
    outcome(r, "yield", "yield", [product_3[label]], [E[3], T[3]], "comparable yields", "comparable yield",
            role="main product", scope={"comparison": "standard-scale reaction; exact value unreported"})
    outcome(r, "stereoselectivity", "selectivity", [product_3[label]], [E[3], T[3]],
            "excellent diastereoselectivities and enantioselectivities", "excellent stereoselectivity",
            role="main product", scope={"metrics": ["dr", "ee"], "exact_values": "unreported"})
    r["review"]["unresolved"] = ["Figure S2a/SI was not supplied; exact gram-scale charge, conditions, yield, dr, ee and measurement methods are unavailable."]
    reactions.append(r); prose_ids.append(r["reaction_id"])

for target, product in [("7", M_7), ("8", M_8)]:
    r = base_reaction(f"derivatization_to_{target}", [E[3], T[3]], "specific_experiment", ["derivatization"],
                      f"Page 3 prose referring to Figure S2b, 3l to {target}", reaction_name="product 3l elaboration",
                      classification=class_empty([E[3]], "unresolved", "Product structure and reagents occur only in unavailable Figure S2b/SI."))
    r["participants"] = [part(f"deriv_{target}_3l", product_3["3l"], "reactant", [E[3]]),
                         part(f"deriv_{target}_prod", product, "product", [E[3]], product_role="derivatization product")]
    outcome(r, "yield", "yield", [product], [E[3], T[3]], "high yields", "high yield",
            role="derivatization product", scope={"exact_value": "unreported"})
    outcome(r, "optical_enrichment", "selectivity", [product], [E[3], T[3]], "optically enriched", "optically enriched",
            role="derivatization product", scope={"metric/value": "unreported"})
    r["review"]["unresolved"] = ["Figure S2b/SI was not supplied; product structure, reagents, conditions, exact yield and selectivity are unavailable."]
    reactions.append(r); prose_ids.append(r["reaction_id"])

# Figure 2b deuterium-labeling experiment.
r = base_reaction("figure2b_deuterium_labeling", [E[4], T[3]], "specific_experiment", ["mechanistic_study"],
                  "Page 4, Figure 2b", classification=class_annulation([E[4]], deuterium=True))
ps = [part("d_1b", substrate_1["1b"], "reactant", [E[4]]), part("d_2a", M_2A, "reactant", [E[4]]),
      part("d_cat", cat_ids["P10"], "catalyst", [E[4]], mol_percent=10, reference=substrate_1["1b"]),
      part("d_base", M_K3PO4, "reagent", [E[4]], equivalents=2.0, reference=substrate_1["1b"]),
      part("d_d2o", M_D2O, "additive", [E[4]], equivalents=3.0, reference=substrate_1["1b"]),
      part("d_xylene", M_XYL, "solvent", [E[4]]), part("d_prod", M_D3B, "product", [E[4]], product_role="deuterium-labeled product")]
r["participants"] = ps
set_solvent(r, M_XYL, [E[4]])
add_step(r, [E[4]], [p["participant_id"] for p in ps], "P10 (10 mol%), K3PO4 (2.0 equiv), D2O (3.0 equiv), xylene, RT, 12 h", hours=12)
outcome(r, "yield", "yield", [M_D3B], [E[4]], "96% yield", q(96, "%", raw_text="96% yield"), role="main product")
outcome(r, "dr", "selectivity", [M_D3B], [E[4]], ">20:1 dr", {"ratio": ">20:1", "metric": "dr"}, role="main product", scope={"metric": "dr"})
outcome(r, "D_incorporation", "other", [M_D3B], [E[4]], "D (62%)", q(62, "%", raw_text="D (62%)"),
        role="deuterium-labeled product", scope={"metric": "deuterium incorporation", "position": "alpha-prime"})
r["review"]["unresolved"] = ["Yield/dr/deuterium-incorporation measurement bases are not printed in Figure 2b.",
                              "Exact substrate charges, xylene volume, workup and purification are unreported without SI."]
reactions.append(r)
deut_id = r["reaction_id"]

# Figure 2c: six distinct control rows.
control_rows = [
    (1, "P10", M_K3PO4, M_XYL, 98, 99, "successful"),
    (2, "P10-1", M_K3PO4, M_XYL, 91, 8, "successful"),
    (3, None, M_K3PO4, M_XYL, 89, 0, "successful"),
    (4, "P10-2", M_K3PO4, M_XYL, 93, 89, "successful"),
    (5, "P10", M_K3PO4, M_MEOH, None, None, "messy"),
    (6, "P10", None, M_XYL, None, None, "no reaction"),
]
control_ids = []
for entry, catalyst, base, solvent, yld, ee, result in control_rows:
    ev = [E[4]]
    observed = result == "successful"
    r = base_reaction(f"figure2c_control_{entry}", ev, "specific_experiment", ["control", "mechanistic_study"],
                      f"Page 4, Figure 2c, entry {entry}",
                      classification=class_annulation(ev) if observed else class_empty(ev, "not_applicable", "No clean product-forming transformation was reported."))
    ps = [part(f"c{entry}_1b", substrate_1["1b"], "reactant", ev), part(f"c{entry}_2a", M_2A, "reactant", ev),
          part(f"c{entry}_solv", solvent, "solvent", ev)]
    if catalyst:
        ps.append(part(f"c{entry}_cat", cat_ids[catalyst], "catalyst", ev, mol_percent=10, reference=substrate_1["1b"]))
    if base:
        ps.append(part(f"c{entry}_base", base, "reagent", ev, equivalents=2.0, reference=substrate_1["1b"]))
    if observed:
        ps.append(part(f"c{entry}_prod", product_3["3b"], "product", ev, product_role="main product"))
    r["participants"] = ps
    set_solvent(r, solvent, ev)
    add_step(r, ev, [p["participant_id"] for p in ps], "cat. (10 mol%), K3PO4 (2.0 equiv), solvent, RT, 8 h", hours=8)
    if observed:
        outcome(r, "yield", "yield", [product_3["3b"]], ev, f"{yld}% yield", q(yld, "%", raw_text=f"{yld}% yield"), role="main product")
        outcome(r, "ee", "selectivity", [product_3["3b"]], ev, f"{ee}% ee", q(ee, "%", raw_text=f"{ee}% ee"),
                role="main product", scope={"metric": "ee", "major_configuration": "not named"})
        outcome(r, "dr", "selectivity", [product_3["3b"]], ev, ">20:1 dr", {"ratio": ">20:1", "metric": "dr"},
                role="main product", scope={"metric": "dr"})
        r["review"]["unresolved"] = ["Yield, ee and dr measurement bases are not printed in Figure 2c."]
    elif result == "messy":
        outcome(r, "messy", "other", [], ev, "messy", "messy reaction mixture", scope={"product_identity": "not assigned"})
        r["review"]["unresolved"] = ["The source reports only 'messy'; no product, conversion, yield, or detection method is assigned."]
    else:
        outcome(r, "no_reaction", "no_reaction", [], ev, "no reaction", "no reaction", scope={"detection_method": "unreported"})
        r["review"]["unresolved"] = ["No-reaction detection method/limit is unreported."]
    reactions.append(r); control_ids.append(r["reaction_id"])

# Figure 2a is an author-proposed mechanism, not experimental proof.
r = base_reaction("figure2a_proposed_mechanism", [E[4], T[3]], "proposed_mechanism", ["mechanistic_study"],
                  "Page 4, Figure 2a", reaction_name="plausible mechanism",
                  classification=class_empty([E[4]], "not_applicable", "Proposed elementary intermediates are not a separate observed synthesis record."))
r["participants"] = [part("mech_1", G1, "reactant", [E[4]]), part("mech_2", G2, "reactant", [E[4]]),
                     part("mech_3", G34, "product", [E[4]], product_role="proposed cycle product")] + [
    part(f"mech_{k}", v, "unknown", [E[4]]) for k, v in mech_mols.items()]
r["notes"] = [f("base-mediated deprotonation; alpha-addition; intramolecular beta-addition; H-shift/isomerization; protonation",
                 ["deprotonation", "alpha-addition", "intramolecular beta-addition", "H-shift/isomerization", "protonation"], [E[4]]),
              f("plausible mechanism", "author proposal, not directly observed proof", [E[4]])]
r["review"]["unresolved"] = ["Intermediates A–E are proposed; no direct observation, kinetic proof, or complete atom mapping is reported on pages 1–5."]
reactions.append(r)
mechanism_id = r["reaction_id"]

# Two generic literature examples from Figure 1a; no current-paper experiment is claimed.
literature_ids = []
for suffix, name, inputs, products in [
    ("figure1a_jorgensen_prior", "chiral ammonium salt-catalyzed conjugate addition (ref. 14a)", prior_mols["prior_Jorgensen_inputs"], prior_mols["prior_Jorgensen_products"]),
    ("figure1a_maruoka_prior", "chiral ammonium salt-catalyzed Mannich reaction and alkylation (ref. 14b)", prior_mols["prior_Maruoka_inputs"], prior_mols["prior_Maruoka_products"]),
]:
    r = base_reaction(suffix, [E[2], T[1]], "literature_example", [], f"Page 2, Figure 1a, {name}", reaction_name=name,
                      classification=class_empty([E[2]], "unresolved", "Generic prior-work structures lack a safe complete atom mapping."))
    r["participants"] = [part(suffix + "_in", inputs, "reactant", [E[2]]), part(suffix + "_out", products, "product", [E[2]], product_role="prior-literature product set")]
    r["review"]["unresolved"] = ["This is a cited prior-work schematic, not an experiment newly reported in this paper; exact conditions/outcomes are not extracted from the cited papers."]
    reactions.append(r); literature_ids.append(r["reaction_id"])

bundle["molecules"] = list(molecules.values())
bundle["reactions"] = reactions
bundle["source_context"] = [
    {"context_id": f"{DOC}:page1_overview", "document_id": DOC, "evidence_ids": [E[1], T[1]],
     "kind": "article_overview", "title": "Enantioselective Synthesis of Multifunctionalized 4H-Pyrans via Formal [4 + 2] Annulation Process by Bifunctional Phosphonium Salt Catalysis",
     "note": "Abstract/graphical abstract summarize up to 99% yield, >20:1 dr and >99.9% ee; these maxima are not duplicated as standalone experiments."},
    {"context_id": f"{DOC}:scheme1_bioactive_examples", "document_id": DOC, "evidence_ids": [E[2]],
     "kind": "non_reaction_context", "note": "Scheme 1 shows representative bioactive 4H-pyran-containing compounds, not experiments conducted in this paper."},
    {"context_id": f"{DOC}:xray_3l", "document_id": DOC, "evidence_ids": [E[3], T[3]],
     "kind": "stereochemical_characterization", "compound": "3l", "method": "single-crystal X-ray diffraction", "ccdc": "1947826"},
    {"context_id": f"{DOC}:supporting_information_notice", "document_id": DOC, "evidence_ids": [E[4], T[4]],
     "kind": "availability_notice", "note": "The article says SI contains experimental procedures and characterization, but no SI file was explicitly supplied or associated in this session."},
    {"context_id": f"{DOC}:page5_references", "document_id": DOC, "evidence_ids": [E[5], T[5]],
     "kind": "bibliography", "note": "Page 5 contains references only; it was inspected and contributes no current-paper synthesis record."},
]

st.write(RUN / "submitted.json", bundle)
st.write(RUN / "raw_model/stereo_checks.json", stereo_diagnostics)
st.write(RUN / "raw_model/source_mentions.json", st.source_mentions(
    DATASET, ["1a", "2a", "3a", "1b", "3b", "3l", "6l", "D2O", "no reaction", "messy"], context=320))

raw_inventory = {
    "document_id": DOC,
    "pages_attempted": [1, 2, 3, 4, 5],
    "table1_rows": [{"entry": e[0], "catalyst": e[1], "solvent": molecules[e[2]]["names"]["original"]["normalized_value"] if e[2] != M_XYL else "xylene",
                     "base": None if e[3] is None else molecules[e[3]]["names"]["original"]["normalized_value"],
                     "yield_percent": e[4], "ee_percent": e[5], "catalyst_mol_percent": e[6]} for e in table1],
    "table2_products": {**{k: {"yield_percent": v[0], "ee_percent": v[1], "dr": ">20:1", "description": v[2]} for k, v in scope_3.items()},
                        **{k: {"yield_percent": v[0], "ee_percent": v[1], "dr": ">20:1", "description": v[2]} for k, v in scope_4.items()}},
    "table3_products": {k: {"yield_percent": v[0], "dr": v[1], "ee_percent": v[2], "description": v[3]} for k, v in scope_6.items()},
    "figure2b": {"product": "d-3b", "yield_percent": 96, "dr": ">20:1", "D_incorporation_percent": 62},
    "figure2c": [{"entry": x[0], "catalyst": x[1], "base": None if x[2] is None else "K3PO4", "solvent": "MeOH" if x[3] == M_MEOH else "xylene",
                  "yield_percent": x[4], "ee_percent": x[5], "result": x[6]} for x in control_rows],
    "specific_experiment_count": sum(r["record_type"] == "specific_experiment" for r in reactions),
    "total_reaction_record_count": len(reactions),
    "limitations": ["No SI was supplied.", "Dense structures remain unresolved/not_encoded where exact stereochemical transcription was unsafe.",
                    "Machine records remain needs_review; visual inspection does not establish chemical accuracy."],
}
st.write(RUN / "raw_model/pages_0001_0005_extraction.json", raw_inventory)

page_records = {
    1: [],
    2: literature_ids + [SCHEME_T1] + table1_ids,
    3: [SCHEME_T2] + table2_ids + [SCHEME_T3] + table3_ids + prose_ids,
    4: [mechanism_id, deut_id] + control_ids,
    5: [],
}
coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {"page": 1, "state": "attempted", "regions": [
            {"region": "title, abstract, and graphical abstract", "finding": "article overview; reported maxima and generic transformation are not standalone experiments", "evidence_ids": [E[1], T[1]]},
            {"region": "two-column introduction", "finding": "background and cited prior work; Table 1 discussion begins at page bottom", "evidence_ids": [E[1], T[1]]}],
         "record_ids": page_records[1], "unresolved": ["Graphical abstract structures are generic and not independently atom-mapped."]},
        {"page": 2, "state": "attempted", "regions": [
            {"region": "Figure 1a", "finding": "two cited prior-work generic reaction schemes, classified as literature examples", "evidence_ids": [E[2], T[1]]},
            {"region": "Figure 1b and Scheme 1", "finding": "current-work overview plus bioactive context; no additional concrete experiment", "evidence_ids": [E[2]]},
            {"region": "Table 1 general arrow, entries 1–15, and footnotes", "finding": "15 optimization experiments including entry 15 no reaction; rows 14–15 retain printed catalyst loading conflict", "evidence_ids": [E[2], T[2]]},
            {"region": "lower-column scope discussion", "finding": "prose summary of Table 2; exact scope rows are recorded from page 3", "evidence_ids": [E[2], T[2]]}],
         "record_ids": page_records[2], "unresolved": ["Exact isomeric catalyst/product SMILES are not safely encoded from the compact drawings."]},
        {"page": 3, "state": "attempted", "regions": [
            {"region": "Table 2 general arrow, panels 3a–y and 4a–e, footnote", "finding": "30 distinct scope experiments with isolated yield, dr and ee", "evidence_ids": [E[3], T[3]]},
            {"region": "Table 3 general arrow, panels 6a–l, footnote", "finding": "12 distinct scope experiments with isolated yield, dr and ee", "evidence_ids": [E[3], T[3]]},
            {"region": "X-ray inset and prose", "finding": "3l stereochemical assignment by X-ray, CCDC 1947826; characterization context, not a synthesis record", "evidence_ids": [E[3], T[3]]},
            {"region": "lower prose referring to Figure S2", "finding": "two gram-scale products and two derivatizations are retained as partial prose-only experiments because SI/Figure S2 was not supplied", "evidence_ids": [E[3], T[3]]},
            {"region": "mechanism discussion at lower right", "finding": "author-proposed cycle continues to Figure 2a on page 4", "evidence_ids": [E[3], T[3]]}],
         "record_ids": page_records[3], "unresolved": ["Dense product stereochemical graphs are not encoded; most exact R/S descriptors are unprinted.", "Figure S2 details and exact gram-scale/derivatization conditions are unavailable."]},
        {"page": 4, "state": "attempted", "regions": [
            {"region": "Figure 2a", "finding": "author-proposed mechanism A–E; retained separately and not treated as experimental proof", "evidence_ids": [E[4], T[3]]},
            {"region": "Figure 2b", "finding": "one deuterium-labeling mechanistic experiment yielding d-3b (96%, >20:1 dr, 62% D)", "evidence_ids": [E[4], T[3]]},
            {"region": "Figure 2c entries 1–6", "finding": "six distinct control experiments, including messy MeOH and no-base no-reaction outcomes", "evidence_ids": [E[4], T[4]]},
            {"region": "Supporting Information and CCDC notices", "finding": "availability/context only; no SI association inferred", "evidence_ids": [E[4], T[4]]},
            {"region": "conclusion and references start", "finding": "summary/bibliography; no additional concrete synthesis arrow", "evidence_ids": [E[4], T[4]]}],
         "record_ids": page_records[4], "unresolved": ["Figure 2b/c analytical bases and full procedures are not printed."]},
        {"page": 5, "state": "attempted", "regions": [
            {"region": "both columns", "finding": "references only; no current-paper reaction arrow, table, or standalone characterization dataset", "evidence_ids": [E[5], T[5]]}],
         "record_ids": page_records[5], "unresolved": []},
    ],
}
st.write(RUN / "coverage/pages_0001_0005.json", coverage)

tasks = json.loads((RUN / "tasks.json").read_text())
for task in tasks:
    if task["document_id"] == DOC and task["evidence_id"] in {T[i] for i in range(1, 6)}:
        task["state"] = "attempted"
        task["reason"] = "Original page image and extracted text inspected; findings mapped in coverage/pages_0001_0005.json."
st.write(RUN / "tasks.json", tasks)

print(json.dumps({"molecules": len(bundle["molecules"]), "reactions": len(reactions),
                  "specific_experiments": raw_inventory["specific_experiment_count"]}, ensure_ascii=False))
