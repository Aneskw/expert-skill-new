#!/usr/bin/env python3
"""Session-local transcription for doc_e1c2e1b1afa31a18.

This program contains source-derived transcription only.  It does not alter the
extraction tools, schema, taxonomy, or quality rules.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

DOC = "doc_e1c2e1b1afa31a18"
RUN = ROOT / "runs/current"
bundle = json.loads((RUN / "submitted.json").read_text())
bundle["molecules"] = []
bundle["reactions"] = []
bundle["source_context"] = []


def ev(page: int, image: bool = True) -> str:
    return f"{DOC}:page{page}:{'image' if image else 'text'}"


def f(raw=None, normalized=None, evidence=(), method="original_extraction", derivation=None,
      status=None):
    return st.fact(raw, normalized, evidence, method=method, derivation=derivation, status=status)


def inferred(value, evidence, derivation):
    return f(None, value, evidence, method="inference", derivation=derivation)


def set_occurrence(m, label, evidence):
    m["occurrences"] = [{"document_id": DOC,
                         "compound_label": f(label, label, evidence),
                         "evidence_ids": list(evidence)}]


def molecule(mid, label, name=None, evidence=(ev(2),), smiles=None,
             chirality="unresolved", geometry="unresolved", representation="not_encoded",
             stereo_note=None):
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    set_occurrence(m, label, evidence)
    m["names"]["original"] = f(name or label, name or label, evidence)
    m["names"]["english"] = f(name or label, name or label, evidence)
    m["names"]["abbreviations"] = f(label, [label], evidence)
    m["names"]["aliases"] = f()
    m["cas_numbers"] = f()
    if smiles:
        m["representation_type"] = "smiles"
        m["smiles"] = f(f"source name/drawing: {name or label}", smiles, evidence,
                          method="structure_derivation",
                          derivation="Candidate SMILES derived from the unambiguous simple source identity; parsing is not source verification.")
        diag = st.stereo_inspect(smiles)
        m["stereochemistry"]["isomeric_smiles"] = f(
            f"candidate from {name or label}", diag["canonical_isomeric_smiles"], evidence,
            method="structure_derivation", derivation="RDKit canonical isomeric SMILES from the candidate structure.")
        m["stereochemistry"]["identity_key"] = f(
            None, diag["stereo_graph_key"], evidence, method="structure_derivation",
            derivation="RDKit stereo graph key from the candidate SMILES; not a sample identity.")
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = f()
        m["stereochemistry"]["isomeric_smiles"] = f()
        m["stereochemistry"]["identity_key"] = f()
    for key, value in (("chirality_status", chirality), ("geometry_status", geometry),
                       ("representation_status", representation)):
        m["stereochemistry"][key] = inferred(value, evidence,
            "Status assigned after inspection of the source drawing/name; unresolved means the supplied pages do not support a complete assignment.")
    if stereo_note:
        if chirality in {"specified", "partially_specified"}:
            m["stereochemistry"]["assignment"] = f(stereo_note, stereo_note, evidence)
        elif geometry in {"specified", "partially_specified"}:
            m["stereochemistry"]["double_bonds"] = f(stereo_note, stereo_note, evidence)
        else:
            m["stereochemistry"]["axial"] = f(stereo_note, stereo_note, evidence)
    m["identity_note"] = f(
        "Structure retained by source label; no SI was supplied and no unsupported complete structure was forced.",
        "Structure retained by source label; no SI supplied.", evidence)
    m["external_properties"] = []
    bundle["molecules"].append(m)
    return mid


def template_molecule(mid, label, fixed, variables, evidence):
    m = molecule(mid, label, label, evidence=evidence, chirality="unresolved",
                 geometry="unresolved", representation="not_encoded")
    m = bundle["molecules"].pop()
    m["representation_type"] = "generic_template"
    m["generic_template"] = {
        "fixed_structure": f(fixed, fixed, evidence),
        "variables": [{
            "label": v,
            "attachment_points": f(f"source drawing attachment for {v}", f"source drawing attachment for {v}", evidence),
            "definition": f(f"source variable {v}", f"source variable {v}", evidence),
            "allowed_values": f(),
            "constraints": f(),
        } for v in variables],
        "instance_mapping": [],
        "reason": f("Source uses R/Ar/PG variables.", "Generic source drawing with unexpanded variables.", evidence),
    }
    bundle["molecules"].append(m)
    return mid


def quantity(number, unit, raw=None, qualifier=None):
    out = {"number": number, "unit": unit}
    if raw is not None:
        out["raw_text"] = raw
    if qualifier:
        out["qualifier"] = qualifier
    return out


def quantity_raw(value):
    if not isinstance(value, dict):
        return str(value)
    if value.get("raw_text"):
        return value["raw_text"]
    if "number" in value:
        return f"{value['number']:g} {value['unit']}"
    return str(value)


def participant(r, pid, mid, role, evidence, *, amount=None, mass=None, volume=None,
                equiv=None, mol_percent=None, reference=None, stage="reaction", product_role=None):
    p = st.new_participant(pid, mid, role, evidence, stage=stage)
    if product_role:
        p["product_role"] = f(product_role, product_role, evidence)
    for key, val in (("amount", amount), ("mass", mass), ("volume", volume),
                     ("equivalents", equiv), ("mol_percent", mol_percent)):
        if val is not None:
            p["charges"][key] = f(quantity_raw(val), val, evidence)
    if reference:
        p["equivalence_reference"] = f(reference, reference, evidence)
    r["participants"].append(p)
    return pid


def classification(r, labels, evidence, note):
    changes = []
    label_objs = []
    for i, code in enumerate(labels, 1):
        if code == "C-C":
            atoms = [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}]
            order = 1
        elif code == "DB-C-O":
            atoms = [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "O"}]
            order = 2
        elif code == "DB-C-N":
            atoms = [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "N"}]
            order = 2
        elif code == "N-H":
            atoms = [{"atom_map": None, "element": "N"}, {"atom_map": None, "element": "H"}]
            order = 1
        else:
            raise ValueError(code)
        cid = f"b{i}"
        changes.append(f(note, {"change_id": cid, "atoms": atoms, "before_order": 0,
                                "after_order": order,
                                "change_type": "hydrogen_gain" if code == "N-H" else "formed"},
                         evidence, method="structure_derivation",
                         derivation=note))
        lab = st.bond_label(atoms[0]["element"], atoms[1]["element"], order)
        label_objs.append({**lab, "supporting_change_ids": [cid]})
    # Merge duplicate category labels while retaining all supporting changes.
    merged = {}
    for lab in label_objs:
        merged.setdefault(lab["code"], {"category_number": lab["category_number"],
                                        "code": lab["code"], "supporting_change_ids": []})
        merged[lab["code"]]["supporting_change_ids"].extend(lab["supporting_change_ids"])
    r["classification"]["mapped_reaction_smiles"] = f()
    r["classification"]["bond_changes"] = changes
    r["classification"]["labels"] = f(note, list(merged.values()), evidence,
                                          method="structure_derivation", derivation=note)
    r["classification"]["coverage"] = inferred("partial", evidence, note)
    r["classification"]["limitations"] = [
        "Partial class-level assignment from the displayed transformation; exact atom maps are unresolved.",
        "Hydrogen transfers, leaving groups, coproducts, and bond-order rearrangements are not claimed complete.",
    ]
    r["classification"]["model_votes"] = []


def base_reaction(rid, evidence, record_type="specific_experiment", categories=("other",),
                  location=None, name=None, transform=None):
    r = st.new_reaction(f"{DOC}:{rid}", DOC, evidence)
    r["record_type"] = record_type
    r["experiment_categories"] = f(list(categories), list(categories), evidence)
    r["location"] = f(location or rid, {"source": location or rid}, evidence)
    r["reaction_name"] = f(name, name, evidence) if name else f()
    r["transformation_type"] = inferred(transform, evidence,
        "Transformation type assigned from the displayed reactant/product classes.") if transform else f()
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                   "checked_fields": [], "unresolved": []}
    r["completeness"] = "partial"
    r["evidence_ids"] = list(evidence)
    return r


def set_solvent(r, raw, normalized, volume, evidence):
    r["solvent_conditions"]["composition"] = f(raw, [normalized], evidence)
    if volume is not None:
        r["solvent_conditions"]["total_volume"] = f(quantity_raw(volume), volume, evidence)


def reaction_step(r, sid, pids, evidence, operation, temperature=None, time=None):
    s = st.new_step(sid, evidence)
    s["participant_ids"] = list(pids)
    s["operation"] = f(operation, operation, evidence)
    if temperature is not None:
        s["temperature"] = f(quantity_raw(temperature), temperature, evidence)
    if time is not None:
        s["time"] = f(quantity_raw(time), time, evidence)
    r["steps"].append(s)


def outcome(r, oid, kind, targets, raw, normalized, evidence, *, basis=None,
            yield_type=None, scope=None, product_role=None):
    o = st.new_outcome(oid, kind, targets, evidence)
    o["value"] = f(raw, normalized, evidence)
    if basis is not None:
        o["measurement_basis"] = f(basis, basis, evidence)
    if yield_type is not None:
        o["yield_type"] = f(yield_type, yield_type, evidence)
    if scope is not None:
        o["scope"] = f(scope, scope, evidence)
    if product_role is not None:
        o["product_role"] = f(product_role, product_role, evidence)
    r["outcomes"].append(o)


# Common simple molecules. Candidate SMILES are intentionally limited to source-unambiguous identities.
simple = [
    ("mol_benzoic_acid", "PhCOOH", "benzoic acid", "O=C(O)c1ccccc1", ev(2)),
    ("mol_ea", "EA", "ethyl acetate", "CCOC(C)=O", ev(2)),
    ("mol_toluene", "toluene", "toluene", "Cc1ccccc1", ev(2)),
    ("mol_acetonitrile", "CH3CN", "acetonitrile", "CC#N", ev(2)),
    ("mol_chloroform", "CHCl3", "chloroform", "ClC(Cl)Cl", ev(2)),
    ("mol_thf", "THF", "tetrahydrofuran", "C1CCOC1", ev(2)),
    ("mol_dioxane", "dioxane", "1,4-dioxane", "O1CCOCC1", ev(2)),
    ("mol_dcm", "DCM", "dichloromethane", "ClCCl", ev(2)),
    ("mol_tfa", "CF3COOH", "trifluoroacetic acid", "O=C(O)C(F)(F)F", ev(2)),
]
for mid, label, name, smiles, evidence in simple:
    molecule(mid, label, name, evidence=(evidence,), smiles=smiles,
             chirality="achiral", geometry="not_applicable", representation="not_applicable")
molecule("mol_sulfur", "S8", "elemental sulfur (S8)", evidence=(ev(2),),
         chirality="unresolved", geometry="not_applicable", representation="not_encoded")
molecule("mol_table1_solvent", "solvent specified", "Table 1 row-specific solvent", evidence=(ev(2),),
         chirality="unresolved", geometry="unresolved", representation="not_encoded")

# Catalyst synthesis and generic route identities.
for mid, label, name, chiral, geom, rep, note in [
    ("mol_phosphole", "starting phosphole", "starting 1-phenylphosphole (6 mL, 32 mmol)", "unresolved", "unresolved", "not_encoded", None),
    ("mol_cinn_acetal", "trans-cinnamaldehyde dimethyl acetal", "trans-cinnamaldehyde dimethyl acetal", "unresolved", "specified", "not_encoded", "trans geometry stated by name"),
    ("mol_1", "1", "exo phosphine sulfide acetal 1", "partially_specified", "unresolved", "not_encoded", "source-drawn exo stereochemistry; absolute configuration not assigned"),
    ("mol_1p", "1′", "endo phosphine sulfide acetal 1′", "partially_specified", "unresolved", "not_encoded", "source-drawn endo stereochemistry; absolute configuration not assigned"),
    ("mol_2", "2", "aldehyde 2", "partially_specified", "unresolved", "not_encoded", "source-drawn bicyclic stereochemistry; absolute configuration not assigned"),
    ("mol_3", "3", "diastereomeric sulfinyl imine mixture 3", "mixture", "unresolved", "not_encoded", "two sets of diastereomeric chiral sulfinyl imines; composition not reported"),
    ("mol_4", "4", "JiaPhos precursor diastereomer 4", "specified", "unresolved", "not_encoded", "source-drawn configuration; 4d absolute structure established by X-ray for one member"),
    ("mol_4p", "4′", "JiaPhos precursor diastereomer 4′", "specified", "unresolved", "not_encoded", "source-drawn configuration; 4d′ absolute structure established by X-ray for one member"),
    ("mol_5", "5", "JiaPhos diastereomer 5", "specified", "unresolved", "not_encoded", "source supplies P/center stereodescriptors, but label-to-descriptor differs between Scheme 2 and prose"),
    ("mol_5p", "5′", "JiaPhos diastereomer 5′", "specified", "unresolved", "not_encoded", "source supplies P/center stereodescriptors, but label-to-descriptor differs between Scheme 2 and prose"),
    ("mol_tbut_sulfinamide", "(RS)-tert-butanesulfinamide", "(RS)-tert-butanesulfinamide", "specified", "not_applicable", "not_encoded", "source explicitly specifies (RS)"),
    ("mol_ti_oipr4", "Ti(OiPr)4", "titanium(IV) isopropoxide", "unresolved", "not_applicable", "not_encoded", None),
    ("mol_raney_ni", "Raney Ni", "Raney nickel", "unresolved", "not_applicable", "not_encoded", None),
]:
    molecule(mid, label, name, evidence=(ev(1), ev(2)), chirality=chiral, geometry=geom,
             representation=rep, stereo_note=note)
for mid, scheme_value, prose_value in [
    ("mol_5", "5 (Sp,R,RS)", "(Rp,S,RS)-5"),
    ("mol_5p", "5′ (Rp,S,RS)", "(Sp,R,RS)-5′"),
]:
    target = next(m for m in bundle["molecules"] if m["molecule_id"] == mid)
    conflict = f([scheme_value, prose_value], None, (ev(2), ev(2, False)), status="conflict")
    conflict["alternatives"] = [
        f(scheme_value, scheme_value, (ev(2),)),
        f(prose_value, prose_value, (ev(2, False),)),
    ]
    target["stereochemistry"]["assignment"] = conflict
    target["stereochemistry"]["atom_centers"] = f(
        "P-chirogenic and carbon-center stereochemistry is explicitly depicted, but the label-to-descriptor string conflicts.",
        "Specified stereogenic P/carbon centers; exact descriptor assignment unresolved by source conflict.",
        (ev(2), ev(2, False)))
template_molecule("mol_armgbr", "ArMgBr", "arylmagnesium bromide with source-defined Ar",
                  ["Ar"], (ev(2),))

# Annulation generic families and instance molecules.
template_molecule("mol_family6", "6", "3-methylene-2-oxindole family 6 with Ar1/R1/PG variables",
                  ["Ar1", "R1", "PG"], (ev(2), ev(3)))
template_molecule("mol_family7", "7", "γ-benzyl allenoate family 7 with Ar2/R2 variables",
                  ["Ar2", "R2"], (ev(2), ev(3)))
template_molecule("mol_family8", "8/8′", "3,3′-spirocyclic oxindole product family with source-drawn stereochemistry",
                  ["Ar1", "Ar2", "R1", "R2", "PG"], (ev(1), ev(3)))

six_labels = list("abcdefghijklmnopqrstuvw")
for lab in six_labels:
    geometry = "specified" if lab == "a" else "partially_specified"
    note = ("6a is named (E)-3-benzylidene-1-tosylindolin-2-one; the source drawing retains the exocyclic alkene geometry"
            if lab == "a" else "Source drawing specifies the exocyclic alkene geometry, but no text R/S-style descriptor is assigned.")
    molecule(f"mol_6{lab}", f"6{lab}", f"methyleneindolinone 6{lab}", evidence=(ev(2), ev(3)),
             chirality="unspecified", geometry=geometry, representation="not_encoded", stereo_note=note)
for lab in list("abcdefghi"):
    molecule(f"mol_7{lab}", f"7{lab}", f"γ-benzyl allenoate 7{lab}", evidence=(ev(2), ev(3)),
             chirality="unspecified", geometry="unresolved", representation="not_encoded",
             stereo_note="Potential allene axial stereochemistry is not assigned in the supplied pages; no racemic composition is inferred.")
for lab in ["a", "b", "c", "d", "e"]:
    molecule(f"mol_5{lab}", f"5{lab}", f"JiaPhos catalyst 5{lab}", evidence=(ev(2), ev(3)),
             chirality="specified", geometry="unresolved", representation="not_encoded",
             stereo_note="Source-drawn P-chirogenic JiaPhos configuration; exact complete isomeric SMILES not encoded.")
    molecule(f"mol_5{lab}p", f"5{lab}′", f"diastereomeric JiaPhos catalyst 5{lab}′", evidence=(ev(2), ev(3)),
             chirality="specified", geometry="unresolved", representation="not_encoded",
             stereo_note="Opposite P-center stereochemistry relative to the unprimed catalyst family is stated; complete isomeric SMILES not encoded.")
molecule("mol_9aa", "9aa", "minor (3 + 2) annulation product 9aa", evidence=(ev(2), ev(3)),
         chirality="unresolved", geometry="unresolved", representation="not_encoded")


def product_molecule(label, primed=False):
    mid = "mol_" + label.replace("′", "p")
    if any(x["molecule_id"] == mid for x in bundle["molecules"]):
        return mid
    molecule(mid, label, f"3,3′-spirocyclic oxindole {label}", evidence=(ev(1), ev(3)),
             chirality="specified", geometry="unresolved", representation="not_encoded",
             stereo_note=("Source-drawn stereobonds define the primed antipodal product series; absolute R/S labels are not printed for each member."
                          if primed else
                          "Source-drawn stereobonds define the unprimed product series; absolute R/S labels are not printed for each member."))
    return mid


# Scheme 2: six displayed reaction branches/arrows.
s2e = (ev(1, False), ev(2), ev(2, False))
r = base_reaction("scheme2_cycloaddition_sulfurization", s2e, categories=("scale_up",),
                  location="Scheme 2, first two operations; page 2", name="JiaPhos scaffold synthesis",
                  transform="Diels–Alder cycloaddition followed by sulfurization")
p0 = participant(r, "phosphole", "mol_phosphole", "reactant", s2e, volume=quantity(6, "mL"), amount=quantity(32, "mmol"))
p1 = participant(r, "acetal", "mol_cinn_acetal", "reactant", s2e)
ps = participant(r, "sulfur", "mol_sulfur", "reagent", s2e)
pt = participant(r, "toluene", "mol_toluene", "solvent", s2e)
pa = participant(r, "product_1", "mol_1", "product", s2e, product_role="major/exo isolated product")
pb = participant(r, "product_1p", "mol_1p", "product", s2e, product_role="minor/endo product")
r["reaction_scale"] = f("starting phosphole: 32 mmol", {"reference_molecule_id":"mol_phosphole", "amount":quantity(32,"mmol")}, s2e)
reaction_step(r, "cycloaddition", [p0, p1], s2e, "cycloaddition", temperature=quantity(140, "°C"))
reaction_step(r, "sulfurization", [ps, pt, pa, pb], s2e, "sulfurization in toluene")
outcome(r, "yield_1", "yield", ["mol_1"], "52%", quantity(52, "%"), s2e, basis="column chromatography", product_role="exo-1")
outcome(r, "amount_1", "product_amount", ["mol_1"], "6.63 g", quantity(6.63, "g"), (ev(2),), product_role="exo-1")
outcome(r, "yield_1p", "yield", ["mol_1p"], "10%", quantity(10, "%"), (ev(2),), product_role="minor 1′")
outcome(r, "combined_yield", "yield", ["mol_1", "mol_1p"], "total of 61% yield", quantity(61, "%"), (ev(1, False),), product_role="exo/endo mixture")
r["review"]["unresolved"] = ["Time, reagent charges, atmosphere, and workup are not reported on supplied pages.",
                                "Scheme component yields sum to 62%, whereas prose reports 61% total; retained as source rounding/statement difference."]
classification(r, ["C-C", "C-C"], s2e, "Two new C–C connections are visible for the cycloaddition; sulfurization changes are not completely mapped.")
bundle["reactions"].append(r)

r = base_reaction("scheme2_deprotection_1_to_2", (ev(1, False), ev(2)), categories=("other",),
                  location="Scheme 2, 1→2; page 2", transform="acetal deprotection to aldehyde")
p1 = participant(r, "substrate_1", "mol_1", "reactant", r["evidence_ids"])
ptfa = participant(r, "tfa", "mol_tfa", "reagent", r["evidence_ids"])
pdcm = participant(r, "dcm", "mol_dcm", "solvent", r["evidence_ids"])
p2 = participant(r, "product_2", "mol_2", "product", r["evidence_ids"], product_role="desired product")
set_solvent(r, "DCM", "dichloromethane", None, r["evidence_ids"])
reaction_step(r, "deprotection", [p1, ptfa, pdcm, p2], r["evidence_ids"], "deprotection", temperature=quantity(0, "°C"))
outcome(r, "yield_2", "yield", ["mol_2"], "94% yield", quantity(94, "%"), r["evidence_ids"], product_role="desired product")
outcome(r, "amount_2", "product_amount", ["mol_2"], "5.51 g", quantity(5.51, "g"), (ev(2),), product_role="desired product")
r["review"]["unresolved"] = ["Time, charges, atmosphere, workup, purification, and yield measurement basis are unreported on supplied pages."]
classification(r, ["DB-C-O"], r["evidence_ids"], "Aldehyde C=O formation is visible; hydrolysis partners and all broken bonds are unresolved.")
bundle["reactions"].append(r)

r = base_reaction("scheme2_condensation_2_to_3", (ev(1, False), ev(2)), categories=("other",),
                  location="Scheme 2, 2→3; page 2", transform="sulfinyl imine condensation")
p2 = participant(r, "substrate_2", "mol_2", "reactant", r["evidence_ids"])
psa = participant(r, "sulfinamide", "mol_tbut_sulfinamide", "reactant", r["evidence_ids"])
pti = participant(r, "ti_oipr4", "mol_ti_oipr4", "reagent", r["evidence_ids"])
pthf = participant(r, "thf", "mol_thf", "solvent", r["evidence_ids"])
p3 = participant(r, "product_3", "mol_3", "product", r["evidence_ids"], product_role="diastereomeric mixture")
set_solvent(r, "THF", "tetrahydrofuran", None, r["evidence_ids"])
reaction_step(r, "condensation", [p2, psa, pti, pthf, p3], r["evidence_ids"], "condensation")
outcome(r, "yield_3", "yield", ["mol_3"], "80% yield", quantity(80, "%"), r["evidence_ids"], product_role="mixture of diastereomers")
outcome(r, "amount_3", "product_amount", ["mol_3"], "5.70 g", quantity(5.70, "g"), (ev(2),), product_role="mixture of diastereomers")
r["review"]["unresolved"] = ["Temperature, time, charges, atmosphere, workup, purification, diastereomer composition, and yield basis are unreported."]
classification(r, ["DB-C-N"], r["evidence_ids"], "Imine C=N formation is visible; dehydration atoms and coproduct are unresolved.")
bundle["reactions"].append(r)

r = base_reaction("scheme2_grignard_3_to_4_4p", (ev(1, False), ev(2)), categories=("other", "scale_up"),
                  location="Scheme 2, 3→4/4′; page 2", transform="diastereoselective Grignard addition to sulfinyl imine")
p3 = participant(r, "substrate_3", "mol_3", "reactant", r["evidence_ids"])
pg = participant(r, "armgbr", "mol_armgbr", "reagent", r["evidence_ids"])
pthf = participant(r, "thf", "mol_thf", "solvent", r["evidence_ids"])
p4 = participant(r, "product_4", "mol_4", "product", r["evidence_ids"], product_role="separated diastereomer")
p4p = participant(r, "product_4p", "mol_4p", "product", r["evidence_ids"], product_role="separated diastereomer")
set_solvent(r, "THF", "tetrahydrofuran", None, r["evidence_ids"])
reaction_step(r, "grignard_addition", [p3, pg, pthf, p4, p4p], r["evidence_ids"], "ArMgBr addition", temperature=quantity(0, "°C"))
outcome(r, "yield_4", "yield", ["mol_4"], "42% yield", quantity(42, "%"), (ev(2),), basis="separated by column chromatography", product_role="diastereomer 4")
outcome(r, "amount_4", "product_amount", ["mol_4"], "2.8 g", quantity(2.8, "g"), (ev(2),), product_role="diastereomer 4")
outcome(r, "yield_4p", "yield", ["mol_4p"], "45% yield", quantity(45, "%"), (ev(2),), basis="separated by column chromatography", product_role="diastereomer 4′")
outcome(r, "amount_4p", "product_amount", ["mol_4p"], "3 g", quantity(3, "g"), (ev(2),), product_role="diastereomer 4′")
r["review"]["unresolved"] = ["Time, quantities, atmosphere, workup, and exact Ar identity for the displayed scale are not fully stated on supplied pages."]
classification(r, ["C-C"], r["evidence_ids"], "One new aryl C–C connection from addition to the imine carbon is visible; exact atom mapping is unresolved.")
bundle["reactions"].append(r)

for suffix, src, prod, yld, mass in [("4_to_5", "mol_4", "mol_5", 40, 2.16),
                                     ("4p_to_5p", "mol_4p", "mol_5p", 37, 2.0)]:
    r = base_reaction(f"scheme2_raney_ni_{suffix}", (ev(1, False), ev(2), ev(2, False)), categories=("other", "scale_up"),
                      location=f"Scheme 2, {suffix.replace('_','→')}; page 2", transform="Raney-Ni reduction/deprotection")
    ps = participant(r, "substrate", src, "reactant", r["evidence_ids"])
    prn = participant(r, "raney_ni", "mol_raney_ni", "reagent", r["evidence_ids"])
    pp = participant(r, "product", prod, "product", r["evidence_ids"], product_role="JiaPhos diastereomer")
    reaction_step(r, "raney_ni_reduction", [ps, prn, pp], r["evidence_ids"], "Raney Ni reduction")
    outcome(r, "yield", "yield", [prod], f"{yld}% yield", quantity(yld, "%"), (ev(2),), product_role="JiaPhos product")
    outcome(r, "amount", "product_amount", [prod], f"{mass:g} g", quantity(mass, "g"), (ev(2),), product_role="JiaPhos product")
    overall = 15.6 if prod == "mol_5" else 14.5
    outcome(r, "overall_five_step_yield", "yield", [prod], f"{overall}% yield by five steps", quantity(overall, "%"),
            (ev(2, False),), yield_type="overall five-step yield",
            scope={"starting_material":"starting phosphole (32.0 mmol)", "ending_product":prod})
    r["notes"].append(f(
        "Scheme 2 and the prose assign opposite Rp/Sp descriptor strings to labels 5 and 5′.",
        "Source conflict retained: label-to-stereodescriptor assignment differs between drawing and prose.",
        (ev(2), ev(2, False))))
    r["review"]["unresolved"] = ["Source conflict in 5/5′ stereodescriptor assignment requires expert review.",
                                  "Temperature, time, solvent, quantities, atmosphere, workup, and yield basis are unreported.",
                                  "Prose reports five-step overall yields 15.6%/14.5%, distinct from the displayed 40%/37% arrow yields."]
    classification(r, ["N-H"], r["evidence_ids"], "Net N–H formation accompanies sulfinyl-group removal; hydrogen source and full mapping are unresolved.")
    bundle["reactions"].append(r)


# Generic displayed annulations for the optimization and two catalyst scope schemes.
for rid, page, catmid, catlabel, solventmid, solventname, allene_amt, catalyst_pct in [
    ("generic_table1_annulation", 2, "mol_5a", "5a/5a′ series", "mol_table1_solvent", "solvent specified by row", 0.2, 20),
    ("generic_scheme3_annulation", 3, "mol_5a", "5a", "mol_ea", "ethyl acetate", 0.2, 20),
    ("generic_scheme4_annulation", 3, "mol_5ap", "5a′", "mol_toluene", "toluene", 0.3, 30),
]:
    evidence = (ev(page), ev(page, False))
    r = base_reaction(rid, evidence, record_type="generic_scheme", categories=("substrate_scope",),
                      location=f"page {page}, displayed general annulation", name="Enantioselective allene–olefin (4 + 2) annulation",
                      transform="formal (4 + 2) annulation")
    p6 = participant(r, "family_6", "mol_family6", "reactant", evidence, amount=quantity(0.1,"mmol"))
    p7 = participant(r, "family_7", "mol_family7", "reactant", evidence, amount=quantity(allene_amt,"mmol"))
    pc = participant(r, "catalyst", catmid, "catalyst", evidence, mol_percent=quantity(catalyst_pct, "mol%"), reference="mol_family6")
    pids = [p6, p7, pc]
    if rid != "generic_table1_annulation":
        pids.append(participant(r, "benzoic_acid", "mol_benzoic_acid", "additive", evidence, amount=quantity(0.02,"mmol")))
    psol = participant(r, "solvent", solventmid, "solvent", evidence, volume=quantity(1.0,"mL"))
    pids.append(psol)
    p8 = participant(r, "family_8", "mol_family8", "product", evidence, product_role="major (4 + 2) product family")
    pids.append(p8)
    set_solvent(r, solventname, solventname, quantity(1.0, "mL"), evidence)
    reaction_step(r, "annulation", pids, evidence, "annulation",
                  temperature=quantity(50, "°C"), time=quantity(12, "h"))
    r["review"]["unresolved"] = ["Generic scheme variables are not concrete experiments; workup, atmosphere, and purification are not reported on supplied pages."]
    if rid == "generic_scheme3_annulation":
        r["notes"].append(f("Caption says substrates ‘1 (0.1 mmol), 2 (0.2 mmol)’ while the scheme labels them 6 and 7.",
                            "Source label conflict retained; normalized participants follow the displayed 6/7 structures.", evidence))
        r["review"]["unresolved"].append("Scheme 3 caption substrate numbering 1/2 conflicts with displayed 6/7.")
    classification(r, ["C-C", "C-C"], evidence, "Formal (4 + 2) annulation visibly forms two new C–C ring connections; exact atom mapping and bond shifts are unresolved.")
    bundle["reactions"].append(r)


# Table 1 optimization rows.
table_rows = [
    (1,"5a",20,"toluene",93,"23:1",99,False,False),
    (2,"5a",20,"EA",99,"25:1",97,False,False),
    (3,"5a",20,"EA",95,"32:1",99,True,False),
    (4,"5b",20,"EA",98,">20:1",98,True,False),
    (5,"5c",20,"EA",95,">20:1",85,True,False),
    (6,"5d",20,"EA",93,">20:1",97,True,False),
    (7,"5e",20,"EA",86,"20:1",86,True,False),
    (8,"5a′",20,"toluene",51,"15:1",-98,False,False),
    (9,"5a′",25,"toluene",65,"14:1",-99,False,False),
    (10,"5a′",30,"toluene",72,"13:1",-98,False,False),
    (11,"5a′",30,"toluene",82,"18:1",-98,False,True),
    (12,"5a′",30,"EA",24,"9:1",-99,False,True),
    (13,"5a′",30,"CH3CN",25,"3:1",-89,False,True),
    (14,"5a′",30,"CHCl3",10,"5:1",-98,False,True),
    (15,"5a′",30,"THF",77,"7:1",-99,False,True),
    (16,"5a′",30,"dioxane",74,"8:1",-99,False,True),
    (17,"5a′",30,"toluene",95,">20:1",-98,True,True),
    (18,"5b′",30,"toluene",82,">20:1",-99,True,True),
    (19,"5c′",30,"toluene",93,">20:1",-99,True,True),
    (20,"5d′",30,"toluene",93,">20:1",-99,True,True),
    (21,"5e′",30,"toluene",62,">20:1",-96,True,True),
]
solvent_ids = {"EA":"mol_ea", "toluene":"mol_toluene", "CH3CN":"mol_acetonitrile",
               "CHCl3":"mol_chloroform", "THF":"mol_thf", "dioxane":"mol_dioxane"}
for entry, cat, pct, solv, yld, ratio, ee, acid, high_allene in table_rows:
    evidence = (ev(2), ev(2, False), ev(3, False)) if entry in (1,2,3,8,11,17) else (ev(2),)
    primed = "′" in cat
    prod = product_molecule("8aa′" if primed else "8aa", primed)
    catmid = "mol_" + cat.replace("′", "p")
    r = base_reaction(f"table1_entry_{entry:02d}", evidence, categories=("optimization",),
                      location=f"Table 1, entry {entry}; page 2", name="Phosphine-catalyzed (4 + 2) annulation of γ-benzyl allenoate 7a",
                      transform="formal (4 + 2) annulation")
    r["general_scheme_ids"] = [f"{DOC}:generic_table1_annulation"]
    p6 = participant(r, "substrate_6a", "mol_6a", "reactant", evidence, amount=quantity(0.1,"mmol"))
    p7 = participant(r, "allenoate_7a", "mol_7a", "reactant", evidence, amount=quantity(0.3 if high_allene else 0.2,"mmol"))
    pc = participant(r, "catalyst", catmid, "catalyst", evidence, mol_percent=quantity(pct,"mol%"), reference="mol_6a")
    psol = participant(r, "solvent", solvent_ids[solv], "solvent", evidence, volume=quantity(1.0,"mL"))
    pids = [p6,p7,pc,psol]
    if acid:
        pids.append(participant(r, "benzoic_acid", "mol_benzoic_acid", "additive", evidence, amount=quantity(0.02,"mmol")))
    pp = participant(r, "main_product", prod, "product", evidence, product_role="major (4 + 2) product")
    pm = participant(r, "minor_product", "mol_9aa", "product", evidence, product_role="minor (3 + 2) product")
    pids += [pp,pm]
    r["reaction_scale"] = f("6a (0.1 mmol)", {"reference_molecule_id":"mol_6a", "amount":quantity(0.1,"mmol")}, evidence)
    set_solvent(r, solv, {"EA":"ethyl acetate"}.get(solv,solv), quantity(1.0,"mL"), evidence)
    reaction_step(r, "annulation", pids, evidence, "annulation", temperature=quantity(50,"°C"), time=quantity(12,"h"))
    outcome(r, "isolated_yield", "yield", [prod], f"{yld}%", quantity(yld,"%"), evidence, basis="isolated yield", yield_type="isolated", product_role="major (4 + 2) product")
    outcome(r, "chemoselectivity", "selectivity", [prod,"mol_9aa"], ratio, ratio, evidence, basis="1H NMR", scope={"metric":"product ratio","order":"8aa/8aa′:9aa"})
    outcome(r, "ee", "selectivity", [prod], f"{ee}% ee", quantity(ee,"%") if ee >= 0 else f"{ee}% ee (reported signed value)", evidence,
            basis="HPLC analysis on a chiral stationary phase", scope={"metric":"ee","sign":"reported signed ee; absolute configuration not inferred"})
    r["review"]["unresolved"] = ["Workup, purification details, atmosphere, and absolute product configuration are not reported on supplied pages.",
                                  "Minor-product absolute amount/yield is not separately reported."]
    classification(r, ["C-C","C-C"], evidence, "Formal (4 + 2) product formation visibly creates two C–C ring connections; exact atom mapping is unresolved.")
    bundle["reactions"].append(r)


# Scope rows transcribed from Scheme 3 and Scheme 4.
scheme3 = [
 ("8aa","6a","7a",95,">20:1",99,None), ("8ba","6b","7a",85,">20:1",99,None),
 ("8ca","6c","7a",89,">20:1",98,None), ("8da","6d","7a",98,">20:1",99,None),
 ("8ea","6e","7a",88,">20:1",97,None), ("8fa","6f","7a",94,">20:1",99,None),
 ("8ga","6g","7a",97,">20:1",99,None), ("8ha","6h","7a",86,">20:1",99,None),
 ("8ia","6i","7a",94,">20:1",99,None), ("8ja","6j","7a",81,">20:1",99,None),
 ("8ka","6k","7a",98,">20:1",99,None), ("8la","6l","7a",95,">20:1",99,None),
 ("8ma","6m","7a",78,">20:1",99,None), ("8na","6n","7a",79,">20:1",99,None),
 ("8oa","6o","7a",93,">20:1",99,None), ("8pa","6p","7a",89,">20:1",99,None),
 ("8qa","6q","7a",93,">20:1",97,None), ("8ra","6r","7a",97,">20:1",99,None),
 ("8sa","6s","7a",74,"8:1",98,None), ("8ta","6t","7a",82,"8:1",80,None),
 ("8ua","6u","7a",76,"16:1",82,None), ("8ab","6a","7b",99,">20:1",99,None),
 ("8ac","6a","7c",96,">20:1",99,None), ("8ad","6a","7d",92,">20:1",99,None),
 ("8ae","6a","7e",93,">20:1",99,None), ("8af","6a","7f",92,">20:1",99,None),
 ("8ag","6a","7g",90,">20:1",99,"80 °C"), ("8ah","6a","7h",99,">20:1",98,"80 °C"),
 ("8ai","6a","7i",66,">20:1",99,"80 °C"),
]
scheme4 = [
 ("8aa′","6a","7a",95,">20:1",-98,None), ("8ba′","6b","7a",99,">20:1",-98,None),
 ("8ca′","6c","7a",79,">20:1",-98,None), ("8da′","6d","7a",94,">20:1",-99,None),
 ("8ea′","6e","7a",99,">20:1",-98,None), ("8fa′","6f","7a",93,">20:1",-91,None),
 ("8ga′","6g","7a",96,">20:1",-98,None), ("8ha′","6h","7a",72,">20:1",-99,None),
 ("8ia′","6i","7a",71,">20:1",-99,None), ("8ka′","6k","7a",98,">20:1",-98,None),
 ("8la′","6l","7a",88,">20:1",-98,None), ("8ma′","6m","7a",71,">20:1",-99,None),
 ("8va′","6v","7a",72,">20:1",-98,None), ("8oa′","6o","7a",91,">20:1",-98,None),
 ("8pa′","6p","7a",80,">20:1",-98,None), ("8qa′","6q","7a",93,">20:1",-98,None),
 ("8ra′","6r","7a",91,">20:1",-98,None), ("8wa′","6w","7a",94,">20:1",-99,None),
 ("8sa′","6s","7a",75,">20:1",-94,None), ("8ta′","6t","7a",79,"10:1",-99,None),
 ("8ua′","6u","7a",88,"15:1",-99,None), ("8ab′","6a","7b",88,">20:1",-99,None),
 ("8ac′","6a","7c",69,">20:1",-99,None), ("8ad′","6a","7d",78,">20:1",-99,None),
 ("8ae′","6a","7e",72,">20:1",-98,"THF"), ("8ag′","6a","7g",88,">20:1",-98,None),
 ("8ah′","6a","7h",90,">20:1",-98,None), ("8ai′","6a","7i",78,">20:1",-97,None),
]


def scope_record(label, s6, s7, yld, dr, ee_value, exception, primed):
    evidence = (ev(3), ev(3, False), ev(4, False)) if primed else (ev(3), ev(3, False))
    prod = product_molecule(label, primed)
    rid = "scheme4_" if primed else "scheme3_"
    r = base_reaction(rid + label.replace("′","p"), evidence, categories=("substrate_scope",),
                      location=("Scheme 4" if primed else "Scheme 3") + f", product {label}; page 3",
                      name="Enantioselective allene–olefin (4 + 2) annulation",
                      transform="formal (4 + 2) annulation")
    r["general_scheme_ids"] = [f"{DOC}:generic_scheme{'4' if primed else '3'}_annulation"]
    p6 = participant(r, "methyleneindolinone", "mol_"+s6, "reactant", evidence, amount=quantity(0.1,"mmol"))
    p7 = participant(r, "allenoate", "mol_"+s7, "reactant", evidence, amount=quantity(0.3 if primed else 0.2,"mmol"))
    catmid = "mol_5ap" if primed else "mol_5a"
    pc = participant(r, "catalyst", catmid, "catalyst", evidence,
                     mol_percent=quantity(30 if primed else 20,"mol%"), reference="mol_"+s6)
    pac = participant(r, "benzoic_acid", "mol_benzoic_acid", "additive", evidence, amount=quantity(0.02,"mmol"))
    solvent = "mol_toluene" if primed else "mol_ea"
    solvent_name = "toluene" if primed else "ethyl acetate"
    if primed and exception == "THF":
        solvent, solvent_name = "mol_thf", "THF"
    psol = participant(r, "solvent", solvent, "solvent", evidence, volume=quantity(1.0,"mL"))
    pp = participant(r, "product", prod, "product", evidence, product_role="major isolated product")
    r["reaction_scale"] = f(f"{s6} (0.1 mmol)", {"reference_molecule_id":"mol_"+s6, "amount":quantity(0.1,"mmol")}, evidence)
    set_solvent(r, solvent_name, solvent_name, quantity(1.0,"mL"), evidence)
    temp = 80 if (not primed and exception == "80 °C") else 50
    reaction_step(r, "annulation", [p6,p7,pc,pac,psol,pp], evidence, "annulation",
                  temperature=quantity(temp,"°C"), time=quantity(12,"h"))
    outcome(r, "isolated_yield", "yield", [prod], f"{yld}% yield", quantity(yld,"%"), evidence,
            basis="isolated yield", yield_type="isolated", product_role="major product")
    outcome(r, "dr", "selectivity", [prod], dr, dr, evidence, basis="1H NMR",
            scope={"metric":"dr","major/minor configurations":"not assigned"})
    outcome(r, "ee", "selectivity", [prod], f"{ee_value}% ee", quantity(ee_value,"%") if ee_value >= 0 else f"{ee_value}% ee (reported signed value)", evidence,
            basis="HPLC analysis", scope={"metric":"ee","signed value retained":"absolute configuration not inferred"})
    if not primed:
        r["notes"].append(f("Scheme 3 caption calls substrates 1 and 2, while the drawing labels them 6 and 7.",
                            "Source label conflict retained; molecule links follow displayed 6/7 labels.", evidence))
    r["review"]["unresolved"] = ["Absolute product configuration, detailed HPLC stationary phase, atmosphere, workup, and purification are not reported on supplied pages.",
                                  "Complex source drawings are retained by label; complete SMILES and atom mapping require expert review."]
    if not primed:
        r["review"]["unresolved"].append("Scheme 3 caption substrate numbering 1/2 conflicts with displayed 6/7.")
    classification(r, ["C-C","C-C"], evidence, "Formal (4 + 2) product formation visibly creates two C–C ring connections; exact atom mapping and hydrogen/bond shifts are unresolved.")
    bundle["reactions"].append(r)


for row in scheme3:
    scope_record(*row, primed=False)

# Separate the parenthesized gram-scale 8aa result from the standard-scale row.
evidence = (ev(3), ev(3, False))
prod = "mol_8aa"
r = base_reaction("scheme3_8aa_scaleup", evidence, categories=("substrate_scope","scale_up"),
                  location="Scheme 3, parenthesized 8aa scale-up result; page 3",
                  name="Enantioselective allene–olefin (4 + 2) annulation",
                  transform="formal (4 + 2) annulation")
r["general_scheme_ids"] = [f"{DOC}:generic_scheme3_annulation"]
p6 = participant(r, "methyleneindolinone", "mol_6a", "reactant", evidence, amount=quantity(1,"mmol"))
p7 = participant(r, "allenoate", "mol_7a", "reactant", evidence, amount=quantity(2,"mmol"))
pc = participant(r, "catalyst", "mol_5a", "catalyst", evidence, mol_percent=quantity(20,"mol%"), reference="mol_6a")
pac = participant(r, "benzoic_acid", "mol_benzoic_acid", "additive", evidence, amount=quantity(0.02,"mmol"))
ps = participant(r, "solvent", "mol_ea", "solvent", evidence, volume=quantity(1,"mL"))
pp = participant(r, "product", prod, "product", evidence, product_role="major isolated product")
r["reaction_scale"] = f("6a (1 mmol)", {"reference_molecule_id":"mol_6a", "amount":quantity(1,"mmol")}, evidence)
set_solvent(r, "EA (caption; footnote f modifies substrate amounts only)", "ethyl acetate", quantity(1,"mL"), evidence)
reaction_step(r, "annulation", [p6,p7,pc,pac,ps,pp], evidence, "annulation", temperature=quantity(50,"°C"), time=quantity(12,"h"))
outcome(r, "isolated_yield", "yield", [prod], "73% yield", quantity(73,"%"), evidence, basis="isolated yield", yield_type="isolated")
outcome(r, "dr", "selectivity", [prod], ">20:1 dr", ">20:1", evidence, basis="1H NMR", scope={"metric":"dr","major/minor":"unassigned"})
outcome(r, "ee", "selectivity", [prod], "99% ee", quantity(99,"%"), evidence, basis="HPLC analysis", scope={"metric":"ee"})
r["notes"].append(f("Caption says 1/2 but drawing and footnote f say 6/7.", "Source numbering conflict retained.", evidence))
r["review"]["unresolved"] = ["Footnote f changes substrate amounts but does not explicitly restate additive/catalyst/solvent scaling; inherited unchanged values require expert review.",
                              "Scheme 3 caption substrate numbering conflict and missing workup/purification remain unresolved."]
classification(r, ["C-C","C-C"], evidence, "Formal (4 + 2) product formation visibly creates two C–C ring connections; exact atom mapping is unresolved.")
bundle["reactions"].append(r)

for row in scheme4:
    scope_record(*row, primed=True)


# Non-reaction/context regions retained separately from formal syntheses.
bundle["source_context"] = [
    {"context_id": f"{DOC}:metadata", "kind": "article_metadata", "evidence_ids": [ev(1), ev(1,False)],
     "value": {"title": "Design of 1-Phosphanorbornene Derivatives as Chiral Organocatalysts for Enantioselective (4 + 2) Annulation Reactions of γ-Benzyl Allenoates",
               "doi": "10.1021/acs.orglett.1c00833"}},
    {"context_id": f"{DOC}:graphical_abstract", "kind": "duplicate_summary", "evidence_ids": [ev(1)],
     "value": "Graphical abstract summarizes the 5a/5a′ scope and antipodal products; it is not an extra experiment."},
    {"context_id": f"{DOC}:scheme1", "kind": "background_and_design", "evidence_ids": [ev(2)],
     "value": "Scheme 1A depicts cited prior catalysts; Scheme 1B is design rationale for the current catalyst family, not a preparative experiment."},
    {"context_id": f"{DOC}:xray", "kind": "characterization", "evidence_ids": [ev(2), ev(4,False)],
     "value": "Single-crystal X-ray structures of 4d and 4d′ and CCDC accessions are characterization evidence, not separate syntheses."},
    {"context_id": f"{DOC}:alkyl_no_product", "kind": "scope_limitation", "evidence_ids": [ev(4), ev(4,False)],
     "value": "Alkyl-substituted methyleneindolinones reacted with γ-benzyl allenoates gave no product; individual labels, quantities, and row-level conditions are not supplied, so no concrete control row is fabricated."},
    {"context_id": f"{DOC}:si_notice", "kind": "supplement_notice", "evidence_ids": [ev(4), ev(4,False)],
     "value": "The article states SI exists, but no SI PDF was provided or explicitly associated in this session."},
    {"context_id": f"{DOC}:references", "kind": "bibliography", "evidence_ids": [ev(4), ev(5), ev(6)],
     "value": "Pages 4–6 bibliography/prior-work citations are source context, not experiments of this paper."},
]

# Raw transcription snapshot for review before merge.
raw = {"document_id": DOC, "table1_rows": table_rows, "scheme3_rows": scheme3,
       "scheme4_rows": scheme4, "scheme2_formal_records": 6,
       "notes": ["Signed negative ee values preserved exactly as reported.",
                 "Scope product configurations are source-drawn but not converted to unsupported R/S SMILES."]}
(RUN / "raw_model/transcription.json").write_text(json.dumps(raw, ensure_ascii=False, indent=2) + "\n")
(RUN / "submitted.json").write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")
print(json.dumps({"molecules": len(bundle["molecules"]), "reactions": len(bundle["reactions"]),
                  "source_context": len(bundle["source_context"])}, ensure_ascii=False))
