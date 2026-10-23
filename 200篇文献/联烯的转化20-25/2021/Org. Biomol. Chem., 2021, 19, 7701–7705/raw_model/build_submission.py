#!/usr/bin/env python3
"""Build the source-derived v3.1 submission for PDF pages 1--5.

The script deliberately keeps unencoded source drawings unresolved.  It does
not perform name-to-structure lookup or promote machine observations to source
review.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st


DOC = "doc_cb508a9ba1bdcee5"
RUN = ROOT / "runs/current"
DATASET = json.loads((RUN / "dataset.json").read_text())
EV = {n: f"{DOC}:page{n}:image" for n in range(1, 6)}
TXT = {n: f"{DOC}:page{n}:text" for n in range(1, 6)}


def f(raw=None, normalized=None, pages=(), method="original_extraction", derivation=None,
      status=None):
    evidence = []
    for p in pages:
        evidence.extend([EV[p], TXT[p]])
    return st.fact(raw, normalized, evidence, method=method, derivation=derivation,
                   status=status)


def inferred(value, pages, derivation):
    return f(None, value, pages, method="inference", derivation=derivation)


def derived(value, pages, derivation):
    return f(None, value, pages, method="structure_derivation", derivation=derivation)


MOLECULES = {}
STEREO_DIAGNOSTICS = {}


def add_molecule(label, pages, name=None, smiles=None, *, chirality="unresolved",
                 geometry="unresolved", representation=None, note=None, mixture=None):
    mid = f"{DOC}:mol:{label}"
    evidence = [item for p in pages for item in (EV[p], TXT[p])]
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f(label, label, pages)
    m["names"]["original"] = f(name or label, name or label, pages)
    m["names"]["english"] = f(name, name, pages) if name else st.fact()
    m["names"]["abbreviations"] = f(label, [label], pages)
    m["names"]["aliases"] = st.fact()
    m["cas_numbers"] = st.fact()
    m["identity_note"] = f(note, note, pages) if note else st.fact()
    if smiles:
        diag = st.stereo_inspect(smiles)
        STEREO_DIAGNOSTICS[mid] = diag
        m["representation_type"] = "smiles"
        m["smiles"] = f(smiles, smiles, pages)
        m["stereochemistry"]["atom_centers"] = derived(
            diag["atom_centers"], pages, "Computed with standard stereo from the source-derived candidate SMILES.")
        m["stereochemistry"]["double_bonds"] = derived(
            diag["double_bonds"], pages, "Computed with standard stereo from the source-derived candidate SMILES.")
        m["stereochemistry"]["axial"] = st.fact()
        m["stereochemistry"]["mixture"] = st.fact()
        m["stereochemistry"]["assignment"] = st.fact()
        m["stereochemistry"]["chirality_status"] = derived(
            chirality, pages, "Status based on the candidate graph; absence of a source assignment is retained.")
        m["stereochemistry"]["geometry_status"] = derived(
            geometry, pages, "Status based on the candidate graph and source drawing.")
        rep = representation or ("partial" if any(x["specified"] == "Unspecified" for x in diag["potential_stereo"]) else "encoded")
        m["stereochemistry"]["representation_status"] = derived(
            rep, pages, "Compared candidate SMILES with standard stereo diagnostics; source identity remains machine-checked only.")
        m["stereochemistry"]["isomeric_smiles"] = derived(
            diag["canonical_isomeric_smiles"], pages, "Canonical isomeric SMILES from standard stereo.")
        m["stereochemistry"]["identity_key"] = derived(
            diag["stereo_graph_key"], pages, "Stereo graph key from standard stereo; not a sample identity.")
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = st.fact()
        m["stereochemistry"]["atom_centers"] = st.fact()
        m["stereochemistry"]["double_bonds"] = st.fact()
        m["stereochemistry"]["axial"] = st.fact()
        m["stereochemistry"]["mixture"] = f(mixture, mixture, pages) if mixture else st.fact()
        m["stereochemistry"]["assignment"] = st.fact()
        m["stereochemistry"]["chirality_status"] = f(
            "No absolute/relative configuration assigned in supplied pages", chirality, pages)
        m["stereochemistry"]["geometry_status"] = f(
            mixture or "Geometry not fully assigned in supplied pages", geometry, pages)
        m["stereochemistry"]["representation_status"] = inferred(
            representation or "not_encoded", pages,
            "The source drawing/label was retained without forcing an uncertain isomeric SMILES.")
        m["stereochemistry"]["isomeric_smiles"] = st.fact()
        m["stereochemistry"]["identity_key"] = st.fact()
    MOLECULES[label] = m
    return mid


def add_template(label, pages, fixed, variables, note):
    mid = f"{DOC}:mol:{label}"
    evidence = [item for p in pages for item in (EV[p], TXT[p])]
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f(label, label, pages)
    m["names"]["original"] = f(label, label, pages)
    m["names"]["english"] = st.fact()
    m["names"]["abbreviations"] = f(label, [label], pages)
    m["names"]["aliases"] = st.fact()
    m["cas_numbers"] = st.fact()
    m["representation_type"] = "generic_template"
    m["smiles"] = st.fact()
    m["generic_template"] = {
        "fixed_structure": f(fixed, fixed, pages),
        "variables": [{
            "label": var,
            "attachment_points": f(desc, desc, pages),
            "definition": f(desc, desc, pages),
            "allowed_values": st.fact(),
            "constraints": st.fact(),
        } for var, desc in variables],
        "instance_mapping": [],
        "reason": f(note, note, pages),
    }
    for key, value in (("chirality_status", "unresolved"),
                       ("geometry_status", "unresolved"),
                       ("representation_status", "not_encoded")):
        m["stereochemistry"][key] = inferred(
            value, pages, "Generic R-group drawing cannot establish a unique stereochemical identity.")
    for key in ("atom_centers", "double_bonds", "axial", "mixture", "assignment",
                "isomeric_smiles", "identity_key"):
        m["stereochemistry"][key] = st.fact()
    m["identity_note"] = f(note, note, pages)
    MOLECULES[label] = m
    return mid


# Source-labelled substrates and products.  Only 1a and 2a are transcribed to
# candidate SMILES because their names and drawings are explicit in the main text.
add_molecule("1a", (1, 2, 3, 4), "2-phenyl-2,3-dihydrophthalazine-1,4-dione",
             "O=C1NN(c2ccccc2)C(=O)c2ccccc21", chirality="achiral",
             geometry="not_applicable", representation="encoded",
             note="Candidate connectivity is derived from the printed name and drawing; not independently structure-verified.")
add_molecule("2a", (1, 2, 3, 4), "buta-2,3-dien-1-ylbenzene", "C=C=CCc1ccccc1",
             chirality="unspecified", geometry="unspecified", representation="partial",
             note="Terminal allene candidate; standard stereo reports an unspecified potential element and cannot certify axial assignment.")
add_molecule("3a", (1, 2, 3, 4), note="Product drawing contains a potentially stereogenic quaternary center and an alkene; no absolute configuration is assigned in supplied pages.")
add_molecule("1e", (2, 3, 4), "2-(p-tolyl)-2,3-dihydrophthalazine-1,4-dione",
             "O=C1NN(c2ccc(C)cc2)C(=O)c2ccccc21", chirality="achiral",
             geometry="not_applicable", representation="encoded",
             note="Candidate connectivity derived from the para-methyl source drawing; not independently structure-verified.")
add_molecule("1f", (2, 3, 4), "2-(4-fluorophenyl)-2,3-dihydrophthalazine-1,4-dione",
             "O=C1NN(c2ccc(F)cc2)C(=O)c2ccccc21", chirality="achiral",
             geometry="not_applicable", representation="encoded",
             note="Candidate connectivity derived from the para-fluoro source drawing; not independently structure-verified.")

for label in [f"1{c}" for c in "bcdghijklmnopqrstuvw"]:
    add_molecule(label, (2, 3, 4) if label in {"1e", "1f"} else (2,),
                 note="Concrete source drawing retained by compound label; complete stereochemical SMILES not asserted.")

for label in [f"3{c}" for c in "bcdefghijklmnopqrstuvw"]:
    add_molecule(label, (2, 3), note="Source product drawing retained; wedge/alkene depiction was not converted to an absolute stereochemical assignment.")
for label in ("3c_prime", "3j_prime", "3k_prime", "3t_prime", "3u_prime"):
    add_molecule(label, (2,), note="Separately drawn regioisomer; retained as a distinct molecule without invented absolute configuration.")

for c in "bcdefghijklmnopq":
    pages = (3,)
    add_molecule(f"2{c}", pages, note="Allene source drawing retained by label; no unique stereochemical SMILES asserted.")
for c in "abcdefghijklmn":
    if c == "m":
        add_molecule("4m", (3,), geometry="mixture", representation="not_encoded",
                     mixture="E/Z = 9:1", note="Reported E/Z mixture; individual E and Z structures are not separately assigned in supplied pages.")
    else:
        add_molecule(f"4{c}", (3,), note="Source product drawing retained; configuration is not fully assigned in supplied pages.")

add_molecule("5", (3,), note="Derivatization product source drawing retained; exact stereochemical SMILES not asserted.")
for label in ("1e_d2", "deuterio_3a", "1a_d5", "3a_d4"):
    add_molecule(label, (3, 4), mixture="Partial deuterium incorporation as reported" if label in {"1e_d2", "deuterio_3a"} else None,
                 note="Isotopologue/partially labelled material retained from Scheme 5 without inventing a unique isotopic microstate.")


# Common reagents/solvents. Candidate SMILES are checked below and remain
# source-identified rather than structure-verified.
simple = {
    "MeCN": ("acetonitrile", "CC#N"),
    "AgOAc": ("silver acetate", "CC(=O)[O-].[Ag+]"),
    "NaOAc": ("sodium acetate", "CC(=O)[O-].[Na+]"),
    "CuOAc2": ("copper(II) acetate", "CC(=O)[O-].CC(=O)[O-].[Cu+2]"),
    "Ag2CO3": ("silver carbonate", "O=C([O-])[O-].[Ag+].[Ag+]"),
    "CsOAc": ("cesium acetate", "CC(=O)[O-].[Cs+]"),
    "K2HPO4": ("dipotassium hydrogen phosphate", "OP(=O)([O-])[O-].[K+].[K+]"),
    "K3PO4": ("tripotassium phosphate", "O=P([O-])([O-])[O-].[K+].[K+].[K+]"),
    "Na2CO3": ("sodium carbonate", "O=C([O-])[O-].[Na+].[Na+]"),
    "NaOPiv": ("sodium pivalate", "CC(C)(C)C(=O)[O-].[Na+]"),
    "DCE": ("1,2-dichloroethane", "ClCCCl"),
    "THF": ("tetrahydrofuran", "C1CCOC1"),
    "toluene": ("toluene", "Cc1ccccc1"),
    "dioxane": ("1,4-dioxane", "O1CCOCC1"),
    "HFIP": ("hexafluoroisopropanol", "OC(C(F)(F)F)C(F)(F)F"),
    "CD3OD": ("methanol-d4", "[2H]OC([2H])([2H])[2H]"),
    "NaIO4": ("sodium periodate", "[Na+].[O-]I(=O)(=O)=O"),
    "CH2Cl2": ("dichloromethane", "ClCCl"),
    "water": ("water", "O"),
    "hexane": ("hexane", "CCCCCC"),
    "EtOAc": ("ethyl acetate", "CCOC(C)=O"),
    "BHT": ("2,6-di-tert-butyl-4-methylphenol", "Cc1cc(C(C)(C)C)c(O)c(C(C)(C)C)c1"),
}
for label, (name, smi) in simple.items():
    pages = (4,) if label in {"CD3OD", "BHT", "hexane", "EtOAc"} else ((3, 4) if label in {"NaIO4", "CH2Cl2", "water"} else (2, 3, 4))
    add_molecule(label, pages, name, smi, chirality="achiral", geometry="not_applicable",
                 representation="encoded", note="Common reagent/solvent candidate from an explicit source label.")

for label, name, pages in (
    ("CpRh_dimer", "[Cp*RhCl2]2", (1, 2, 3, 4)),
    ("CpIr_dimer", "[Cp*IrCl2]2", (2,)),
    ("Ru_pcymene_dimer", "[Ru(p-cymene)Cl2]2", (2,)),
    ("RuCl3", "ruthenium(III) chloride", (3,)),
    ("TEMPO", "TEMPO", (4,)),
    ("silica", "silica gel (200-300 mesh)", (4,)),
):
    add_molecule(label, pages, name=name, note="Explicit source reagent; unique molecular representation not asserted.")


# Generic R-group structures used only by general/prior-work schemes.
add_template("generic_1", (2,), "N-aryl phthalazinone/pyridazinone scaffold (X = C or N)",
             [("R1", "substituent shown on the heterocyclic scaffold"), ("R2", "substituent shown on the N-aryl ring"), ("X", "C or N")],
             "Scheme 1/2 generic R-group substrate; not a concrete molecule.")
add_template("generic_2", (2, 3), "substituted allene scaffold", [("R", "allene substituent")],
             "Scheme 1/3 generic allene; not a concrete molecule.")
add_template("generic_3", (2,), "annulated indazole product scaffold (X = C or N)",
             [("R1", "source-defined scaffold substituent"), ("R2", "source-defined aryl substituent"), ("X", "C or N")],
             "Scheme 1/2 generic product; configuration not assigned.")
add_template("generic_4", (3,), "annulated product bearing allene-derived substituent R1", [("R1", "aryl/aliphatic substituent varied in Scheme 3")],
             "Scheme 3 generic product; not a concrete molecule.")

for i in range(1, 4):
    add_template(f"prior_reactant_{i}", (2,), f"Scheme 1 prior-work generic reactant set {i}", [("R1", "source R-group"), ("R2", "source R-group")], "Prior-work drawing is generic and is not a current-paper experiment.")
    add_template(f"prior_product_{i}", (2,), f"Scheme 1 prior-work generic product set {i}", [("R1", "source R-group"), ("R2", "source R-group")], "Prior-work drawing is generic and is not a current-paper experiment.")


def mid(label):
    return MOLECULES[label]["molecule_id"]


def participant(pid, label, role, pages, *, amount=None, mass=None, volume=None,
                equiv=None, molpct=None, reference=None, stage="reaction", product_role=None):
    evidence = [item for p in pages for item in (EV[p], TXT[p])]
    p = st.new_participant(pid, mid(label), role, evidence, stage=stage)
    p["product_role"] = f(product_role, product_role, pages) if product_role else st.fact()
    vals = {"amount": amount, "mass": mass, "volume": volume, "equivalents": equiv, "mol_percent": molpct}
    for key, value in vals.items():
        if value is not None:
            p["charges"][key] = f(value.get("raw_text", value), value, pages)
    if reference:
        p["equivalence_reference"] = f(reference, mid(reference), pages)
    return p


def step(step_id, pages, pids, operation, temperature=None, time=None, atmosphere=None,
         stage="reaction"):
    evidence = [item for p in pages for item in (EV[p], TXT[p])]
    s = st.new_step(step_id, evidence, stage=stage)
    s["participant_ids"] = list(pids)
    s["operation"] = f(operation, operation, pages)
    if temperature is not None:
        s["temperature"] = f(temperature.get("raw_text", temperature) if isinstance(temperature, dict) else temperature, temperature, pages)
    if time is not None:
        s["time"] = f(time.get("raw_text", time) if isinstance(time, dict) else time, time, pages)
    if atmosphere is not None:
        s["atmosphere"] = f(atmosphere, atmosphere, pages)
    return s


def outcome(oid, kind, targets, pages, value=None, yield_type=None, basis=None,
            role="product", scope=None):
    evidence = [item for p in pages for item in (EV[p], TXT[p])]
    o = st.new_outcome(oid, kind, [mid(t) for t in targets], evidence)
    o["product_role"] = f(role, role, pages)
    if value is not None:
        o["value"] = f(value.get("raw_text", value) if isinstance(value, dict) else value, value, pages)
    if yield_type is not None:
        o["yield_type"] = f(yield_type, yield_type, pages)
    if basis is not None:
        o["measurement_basis"] = f(basis, basis, pages)
    if scope is not None:
        o["scope"] = f(scope, scope, pages)
    return o


def classification(kind, pages):
    evidence = [item for p in pages for item in (EV[p], TXT[p])]
    if kind in {"annulation", "annulation_deuteration"}:
        changes = [
            {"change_id": "b1", "atoms": [{"element": "C", "source_role": "N-aryl ortho carbon"}, {"element": "C", "source_role": "allene-derived carbon"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
            {"change_id": "b2", "atoms": [{"element": "C", "source_role": "allene-derived carbon"}, {"element": "N", "source_role": "substrate N-H nitrogen"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
        ]
        labels = [
            {"category_number": 2, "code": "C-C", "supporting_change_ids": ["b1"]},
            {"category_number": 5, "code": "C-N", "supporting_change_ids": ["b2"]},
        ]
        limitation = "Partial net classification from source drawings: two ring-closing sigma bonds are represented; exact atom maps, all allene bond-order shifts, implicit-H losses, and coproducts remain unresolved."
        if kind == "annulation_deuteration":
            changes += [
                {"change_id": "b3", "atoms": [{"element": "C", "source_role": "product methyl carbon"}, {"element": "H", "isotope": 2}], "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain"},
                {"change_id": "b4", "atoms": [{"element": "C", "source_role": "product N-aryl ortho carbon"}, {"element": "H", "isotope": 2}], "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain"},
            ]
            labels.append({"category_number": 1, "code": "C-H", "supporting_change_ids": ["b3", "b4"]})
            limitation += " Two fractional C-D incorporation sites are represented under the C-H taxonomy; individual isotopic microstates remain unresolved."
        coverage = "partial"
    elif kind == "deuteration":
        changes = [{"change_id": "b1", "atoms": [{"element": "C"}, {"element": "H", "isotope": 2}], "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain"}]
        labels = [{"category_number": 1, "code": "C-H", "supporting_change_ids": ["b1"]}]
        limitation = "Partial isotope-exchange classification; individual isotopic microstates and broken C-H/C-D bonds are not fully mapped."
        coverage = "partial"
    elif kind == "oxidation":
        changes = [{"change_id": "b1", "atoms": [{"element": "C"}, {"element": "O"}], "before_order": 0, "after_order": 2, "change_type": "formed"}]
        labels = [{"category_number": 14, "code": "DB-C-O", "supporting_change_ids": ["b1"]}]
        limitation = "Partial drawing-derived classification of the new carbonyl; cleavage partners and exact atom mapping are unresolved."
        coverage = "partial"
    elif kind == "not_applicable":
        changes, labels, limitation, coverage = [], [], "No observed product transformation was reported.", "not_applicable"
    else:
        changes, labels, limitation, coverage = [], None, "Source does not support a complete net bond-change assignment.", "unresolved"
    return {
        "standard_id": st.STANDARD,
        "mapped_reaction_smiles": st.fact(),
        "bond_changes": [st.fact(None, c, evidence, method="structure_derivation", derivation=limitation) for c in changes],
        "labels": st.fact(None, labels, evidence, method="structure_derivation", derivation=limitation) if labels is not None else st.fact(),
        "coverage": st.fact(None, coverage, evidence, method="inference", derivation=limitation),
        "limitations": [limitation],
        "model_votes": [],
    }


REACTIONS = []


def reaction(suffix, pages, record_type="specific_experiment", categories=(), *, location,
             name="Rh(III)-catalyzed [4 + 1] annulation", transformation="C-H oxidative cyclization/annulation",
             class_kind="annulation", unresolved=()):
    rid = f"{DOC}:{suffix}"
    evidence = [item for p in pages for item in (EV[p], TXT[p])]
    r = st.new_reaction(rid, DOC, evidence)
    r["record_type"] = record_type
    r["experiment_categories"] = f(list(categories), list(categories), pages) if categories else st.fact()
    r["location"] = f(location, {"pages": list(pages), "location": location}, pages)
    r["reaction_name"] = f(name, name, pages) if name else st.fact()
    r["transformation_type"] = inferred(transformation, pages, "Normalized from the source reaction description and drawings.") if transformation else st.fact()
    r["classification"] = classification(class_kind, pages)
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                   "checked_fields": [], "unresolved": list(unresolved)}
    r["completeness"] = "partial"
    REACTIONS.append(r)
    return r


def set_solvent(r, label, pages, volume=None, mixture_ratio=None):
    r["solvent_conditions"]["composition"] = f(label, [label], pages)
    if volume is not None:
        r["solvent_conditions"]["total_volume"] = f(volume.get("raw_text", volume), volume, pages)
    if mixture_ratio is not None:
        r["solvent_conditions"]["mixture_ratio"] = f(mixture_ratio, mixture_ratio, pages)


def set_scale(r, reference, amount, pages):
    r["reaction_scale"] = f(amount.get("raw_text", amount), {"reference_molecule_id": mid(reference), "amount": amount}, pages)


# Scheme 1: three literature examples and the current-work general concept.
for i, citation in ((1, "Wu and Ji, J. Org. Chem. 2018"), (2, "Karishma et al., J. Org. Chem. 2018"), (3, "Karishma et al., Adv. Synth. Catal. 2021")):
    r = reaction(f"scheme1_prior_{i}", (2,), "literature_example", location=f"Scheme 1, eqn ({i}), previous work",
                 name="Prior-work transformation reproduced in Scheme 1", transformation=None,
                 class_kind="unresolved", unresolved=("Generic prior-work drawing; conditions and concrete identities belong to the cited paper and were not reconstructed here.",))
    r["participants"] = [participant("reactants", f"prior_reactant_{i}", "reactant", (2,)), participant("product", f"prior_product_{i}", "product", (2,), product_role="reported prior-work product")]
    r["notes"] = [f(citation, citation, (2,))]

r_scheme1 = reaction("scheme1_this_work", (2,), "generic_scheme", location="Scheme 1, This work",
                     unresolved=("Generic R-group transformation; concrete experiments are recorded separately.",))
r_scheme1["participants"] = [participant("substrate", "generic_1", "reactant", (2,)), participant("allene", "generic_2", "reactant", (2,)), participant("product", "generic_3", "product", (2,), product_role="generic annulation product")]

r_scheme2 = reaction("scheme2_general", (2,), "generic_scheme", location="Scheme 2 substrate-scope general arrow",
                     unresolved=("Generic R-group transformation; each displayed product is recorded separately.",))
r_scheme2["participants"] = [participant("substrate", "generic_1", "reactant", (2,), amount={"number": 0.2, "unit": "mmol", "raw_text": "0.2 mmol"}), participant("allene", "2a", "reactant", (2,), amount={"number": 0.3, "unit": "mmol", "raw_text": "0.3 mmol"}), participant("product", "generic_3", "product", (2,), product_role="generic scope product"), participant("cat", "CpRh_dimer", "catalyst", (2,), molpct={"number": 2.5, "unit": "mol%", "raw_text": "2.5 mol%"}, reference="generic_1"), participant("ox", "AgOAc", "reagent", (2,), equiv={"number": 1.0, "unit": "equiv", "raw_text": "1.0 equiv."}, reference="generic_1"), participant("add", "NaOAc", "additive", (2,), equiv={"number": 1.0, "unit": "equiv", "raw_text": "1.0 equiv."}, reference="generic_1"), participant("solv", "MeCN", "solvent", (2,), volume={"number": 1.0, "unit": "mL", "raw_text": "1.0 mL"})]
r_scheme2["steps"] = [step("reaction", (2,), [p["participant_id"] for p in r_scheme2["participants"]], "stirred under the Scheme 2 conditions", {"number": 120, "unit": "°C", "raw_text": "120 °C"}, {"number": 12, "unit": "h", "raw_text": "12 h"}, "air")]
set_solvent(r_scheme2, "MeCN", (2,), {"number": 1.0, "unit": "mL", "raw_text": "1.0 mL"})

r_scheme3 = reaction("scheme3_general", (3,), "generic_scheme", location="Scheme 3 allene-scope general arrow",
                     unresolved=("Generic allene substituent R1; concrete experiments are recorded separately.",))
r_scheme3["participants"] = [participant("substrate", "1a", "reactant", (3,), amount={"number": 0.2, "unit": "mmol", "raw_text": "0.2 mmol"}), participant("allene", "generic_2", "reactant", (3,), amount={"number": 0.3, "unit": "mmol", "raw_text": "0.3 mmol"}), participant("product", "generic_4", "product", (3,), product_role="generic scope product"), participant("cat", "CpRh_dimer", "catalyst", (3,), molpct={"number": 2.5, "unit": "mol%", "raw_text": "2.5 mol%"}, reference="1a"), participant("ox", "AgOAc", "reagent", (3,), equiv={"number": 1.0, "unit": "equiv", "raw_text": "1.0 equiv."}, reference="1a"), participant("add", "NaOAc", "additive", (3,), equiv={"number": 1.0, "unit": "equiv", "raw_text": "1.0 equiv."}, reference="1a"), participant("solv", "MeCN", "solvent", (3,), volume={"number": 1.0, "unit": "mL", "raw_text": "1.0 mL"})]
r_scheme3["steps"] = [step("reaction", (3,), [p["participant_id"] for p in r_scheme3["participants"]], "stirred under the Scheme 3 conditions", {"number": 120, "unit": "°C", "raw_text": "120 °C"}, {"number": 12, "unit": "h", "raw_text": "12 h"}, "air")]
set_solvent(r_scheme3, "MeCN", (3,), {"number": 1.0, "unit": "mL", "raw_text": "1.0 mL"})


# Table 1 optimization rows. Entry 19 is also Scheme 2 product 3a and the
# detailed general procedure on page 4; those evidences are merged here.
table_rows = [
    (1, "CpRh_dimer", "AgOAc", None, "MeCN", 83, 100),
    (2, "CpIr_dimer", "AgOAc", None, "MeCN", 35, 100),
    (3, "Ru_pcymene_dimer", "AgOAc", None, "MeCN", "Trace", 100),
    (4, "CpRh_dimer", "CuOAc2", None, "MeCN", 52, 100),
    (5, "CpRh_dimer", "Ag2CO3", None, "MeCN", 67, 100),
    (6, "CpRh_dimer", "AgOAc", "NaOAc", "MeCN", 85, 100),
    (7, "CpRh_dimer", "AgOAc", "CsOAc", "MeCN", 55, 100),
    (8, "CpRh_dimer", "AgOAc", "K2HPO4", "MeCN", 82, 100),
    (9, "CpRh_dimer", "AgOAc", "K3PO4", "MeCN", 45, 100),
    (10, "CpRh_dimer", "AgOAc", "Na2CO3", "MeCN", 65, 100),
    (11, "CpRh_dimer", "AgOAc", "NaOPiv", "MeCN", 74, 100),
    (12, "CpRh_dimer", "AgOAc", "NaOAc", "DCE", 46, 100),
    (13, "CpRh_dimer", "AgOAc", "NaOAc", "THF", 35, 100),
    (14, "CpRh_dimer", "AgOAc", "NaOAc", "toluene", 42, 100),
    (15, "CpRh_dimer", "AgOAc", "NaOAc", "dioxane", 27, 100),
    (16, "CpRh_dimer", "AgOAc", "NaOAc", "HFIP", 24, 100),
    (17, "CpRh_dimer", "AgOAc", "NaOAc", "MeCN", 84, 100),
    (18, "CpRh_dimer", "AgOAc", "NaOAc", "MeCN", 85, 100),
    (19, "CpRh_dimer", "AgOAc", "NaOAc", "MeCN", 88, 120),
]

for entry, cat, oxidant, additive, solvent, yld, temp in table_rows:
    pages = (2, 4) if entry == 19 else ((1, 2) if entry == 1 else (2,))
    categories = ("optimization", "substrate_scope") if entry == 19 else ("optimization",)
    r = reaction(f"table1_entry_{entry:02d}", pages, categories=categories,
                 location=f"Table 1 entry {entry}" + ("; Scheme 2 product 3a; Scheme 5 standard row; Experimental section" if entry == 19 else ""),
                 unresolved=("No explicitly supplied SI; compound-specific characterization is unavailable.",))
    cat_pct = 5.0 if entry == 18 else 2.5
    ox_eq = 2.0 if entry == 17 else (1.0 if oxidant == "AgOAc" else None)
    r["participants"] = [
        participant("substrate", "1a", "reactant", pages, amount={"number": 0.20, "unit": "mmol", "raw_text": "0.20 mmol"}, mass={"number": 47.6, "unit": "mg", "raw_text": "47.6 mg"} if entry == 19 else None),
        participant("allene", "2a", "reactant", pages, amount={"number": 0.30, "unit": "mmol", "raw_text": "0.30 mmol"}, mass={"number": 39.1, "unit": "mg", "raw_text": "39.1 mg"} if entry == 19 else None),
        participant("cat", cat, "catalyst", pages, molpct={"number": cat_pct, "unit": "mol%", "raw_text": f"{cat_pct:g} mol%"}, reference="1a", mass={"number": 3.1, "unit": "mg", "raw_text": "3.1 mg"} if entry == 19 else None),
        participant("oxidant", oxidant, "reagent", pages, equiv={"number": ox_eq, "unit": "equiv", "raw_text": f"{ox_eq:g} equiv."} if ox_eq is not None else None, reference="1a" if ox_eq is not None else None, amount={"number": 0.20, "unit": "mmol", "raw_text": "0.20 mmol"} if entry == 19 else None, mass={"number": 33.4, "unit": "mg", "raw_text": "33.4 mg"} if entry == 19 else None),
        participant("solvent", solvent, "solvent", pages, volume={"number": 1.0, "unit": "mL", "raw_text": "1.0 mL"}),
        participant("product", "3a", "product", pages, product_role="desired product"),
    ]
    if additive:
        add_eq = 1.0 if additive == "NaOAc" else None
        r["participants"].insert(-1, participant("additive", additive, "additive", pages,
            equiv={"number": add_eq, "unit": "equiv", "raw_text": "1.0 equiv."} if add_eq else None,
            reference="1a" if add_eq else None,
            amount={"number": 0.20, "unit": "mmol", "raw_text": "0.20 mmol"} if entry == 19 else None,
            mass={"number": 16.4, "unit": "mg", "raw_text": "16.4 mg"} if entry == 19 else None))
    if entry == 19:
        r["participants"].extend([
            participant("silica", "silica", "purification_agent", (4,), stage="purification"),
            participant("hexane", "hexane", "purification_agent", (4,), stage="purification"),
            participant("EtOAc", "EtOAc", "purification_agent", (4,), stage="purification"),
        ])
    set_scale(r, "1a", {"number": 0.20, "unit": "mmol", "raw_text": "0.20 mmol"}, pages)
    set_solvent(r, solvent, pages, {"number": 1.0, "unit": "mL", "raw_text": "1.0 mL"})
    all_reaction_pids = [p["participant_id"] for p in r["participants"] if p["stage"] == "reaction"]
    r["steps"] = [step("reaction", pages, all_reaction_pids, "stirred", {"number": temp, "unit": "°C", "raw_text": f"{temp} °C"}, {"number": 12, "unit": "h", "raw_text": "12 h"}, "air")]
    if entry == 19:
        r["steps"] += [
            step("cool", (4,), all_reaction_pids, "cooled to room temperature", "room temperature"),
            step("evaporate", (4,), ["solvent"], "solvent removed under vacuum", stage="workup"),
            step("chromatography", (4,), ["product", "silica", "hexane", "EtOAc"], "silica gel (200-300 mesh) column chromatography using hexane/ethyl acetate (6:1)", stage="purification"),
        ]
        r["procedure"]["original_text"] = f("A mixture of 1a (47.6 mg, 0.20 mmol), 2a (39.1 mg, 0.30 mmol), [Rh*CpCl2]2 (3.1 mg, 2.5 mol%), NaOAc (16.4 mg, 0.20 mmol) and AgOAc (33.4 mg, 0.20 mmol) in MeCN (1.0 mL) was stirred at 120 °C under open air for 12 hours; cooled, evaporated, and purified on silica gel with hexane/ethyl acetate (6:1).", "See raw transcription", (4,))
        r["procedure"]["addition_order"] = f("A mixture of all listed reaction components", ["1a", "2a", "[Rh*CpCl2]2", "NaOAc", "AgOAc", "MeCN"], (4,))
        r["notes"] = [f("The same 88% standard-condition outcome appears in Table 1 entry 19, Scheme 2 as 3a, and Scheme 5's uninhibited standard row; these were merged as one experiment.", "duplicate source mentions merged", (2, 4))]
    yield_pages = (2, 4) if entry == 19 else (2,)
    if yld == "Trace":
        r["outcomes"] = [outcome("trace", "trace", ["3a"], yield_pages, "Trace", role="desired product", scope="reported product 3a")]
    else:
        r["outcomes"] = [outcome("yield", "yield", ["3a"], yield_pages, {"number": yld, "unit": "%", "raw_text": f"{yld}%"}, "isolated yield", "isolated", role="desired product")]
    r["general_scheme_ids"] = [r_scheme1["reaction_id"], r_scheme2["reaction_id"]] if entry == 19 else [r_scheme1["reaction_id"]]


# Scheme 2 scope: 3a is already the standard Table 1 entry 19.
scope2 = {
    "3b": (26, None), "3c": (70, ("3c", "3c_prime", [9, 1])), "3d": (79, None),
    "3e": (76, None), "3f": (71, None), "3g": (72, None), "3h": (47, None),
    "3i": (62, None), "3j": (87, ("3j", "3j_prime", [1, 1])),
    "3k": (77, ("3k", "3k_prime", [1, 1])), "3l": (74, None),
    "3m": (72, None), "3n": (62, None), "3o": (61, None), "3p": (65, None),
    "3q": (63, None), "3r": (59, None), "3s": (59, None),
    "3t": (55, ("3t", "3t_prime", [9, 1])),
    "3u": (58, ("3u", "3u_prime", [5, 1])), "3v": (68, None), "3w": (52, None),
}
for prod, (yld, ratio) in scope2.items():
    substrate = "1" + prod[1]
    targets = [prod] if ratio is None else [ratio[0], ratio[1]]
    r = reaction(f"scheme2_{prod}", (2, 3) if prod in {"3b", "3c", "3d", "3e", "3f", "3g", "3h", "3i", "3j", "3k", "3l", "3m", "3n", "3o", "3p", "3q", "3r", "3s", "3t", "3u", "3v", "3w"} else (2,), categories=("substrate_scope",), location=f"Scheme 2 product {prod}", unresolved=("No explicitly associated SI; detailed workup, characterization and full configuration assignment are unavailable.",))
    r["participants"] = [
        participant("substrate", substrate, "reactant", (2,), amount={"number": 0.2, "unit": "mmol", "raw_text": "0.2 mmol"}),
        participant("allene", "2a", "reactant", (2,), amount={"number": 0.3, "unit": "mmol", "raw_text": "0.3 mmol"}),
        participant("cat", "CpRh_dimer", "catalyst", (2,), molpct={"number": 2.5, "unit": "mol%", "raw_text": "2.5 mol%"}, reference=substrate),
        participant("oxidant", "AgOAc", "reagent", (2,), equiv={"number": 1.0, "unit": "equiv", "raw_text": "1.0 equiv."}, reference=substrate),
        participant("additive", "NaOAc", "additive", (2,), equiv={"number": 1.0, "unit": "equiv", "raw_text": "1.0 equiv."}, reference=substrate),
        participant("solvent", "MeCN", "solvent", (2,), volume={"number": 1.0, "unit": "mL", "raw_text": "1.0 mL"}),
    ] + [participant("product" + str(i + 1), t, "product", (2,), product_role="desired product" if i == 0 else "regioisomer") for i, t in enumerate(targets)]
    set_scale(r, substrate, {"number": 0.2, "unit": "mmol", "raw_text": "0.2 mmol"}, (2,))
    set_solvent(r, "MeCN", (2,), {"number": 1.0, "unit": "mL", "raw_text": "1.0 mL"})
    r["steps"] = [step("reaction", (2,), [p["participant_id"] for p in r["participants"] if p["stage"] == "reaction"], "stirred under Scheme 2 conditions", {"number": 120, "unit": "°C", "raw_text": "120 °C"}, {"number": 12, "unit": "h", "raw_text": "12 h"}, "air")]
    r["outcomes"] = [outcome("combined_yield" if ratio else "yield", "yield", targets, (2,), {"number": yld, "unit": "%", "raw_text": f"{yld}%"}, "isolated yield", "isolated", role="desired product(s)")]
    if ratio:
        a, b, nums = ratio
        r["outcomes"].append(outcome("regioisomer_ratio", "selectivity", [a, b], (2,), {"ratio": nums, "order": [a.replace("_prime", "′"), b.replace("_prime", "′")], "raw_text": f"{nums[0]}:{nums[1]}"}, role="regioisomer mixture", scope={"metric": "regioisomer ratio", "order": [a, b]}))
    r["general_scheme_ids"] = [r_scheme2["reaction_id"]]


# Scheme 3 scope. 4a--4n correspond to 2b--2o; 2p and 2q are N.R.
scope3_yields = dict(zip([f"4{c}" for c in "abcdefghijklmn"], [90, 75, 83, 66, 77, 64, 80, 74, 69, 61, 52, 73, 83, 46]))
for idx, (prod, yld) in enumerate(scope3_yields.items(), start=0):
    allene = "2" + chr(ord("b") + idx)
    r = reaction(f"scheme3_{prod}", (3,), categories=("substrate_scope",), location=f"Scheme 3 product {prod}", unresolved=("Scheme 3 does not state the yield measurement basis; retained as unreported.", "No explicitly associated SI; detailed workup and full configuration assignment are unavailable."))
    r["participants"] = [
        participant("substrate", "1a", "reactant", (3,), amount={"number": 0.2, "unit": "mmol", "raw_text": "0.2 mmol"}),
        participant("allene", allene, "reactant", (3,), amount={"number": 0.3, "unit": "mmol", "raw_text": "0.3 mmol"}),
        participant("cat", "CpRh_dimer", "catalyst", (3,), molpct={"number": 2.5, "unit": "mol%", "raw_text": "2.5 mol%"}, reference="1a"),
        participant("oxidant", "AgOAc", "reagent", (3,), equiv={"number": 1.0, "unit": "equiv", "raw_text": "1.0 equiv."}, reference="1a"),
        participant("additive", "NaOAc", "additive", (3,), equiv={"number": 1.0, "unit": "equiv", "raw_text": "1.0 equiv."}, reference="1a"),
        participant("solvent", "MeCN", "solvent", (3,), volume={"number": 1.0, "unit": "mL", "raw_text": "1.0 mL"}),
        participant("product", prod, "product", (3,), product_role="desired product"),
    ]
    set_scale(r, "1a", {"number": 0.2, "unit": "mmol", "raw_text": "0.2 mmol"}, (3,))
    set_solvent(r, "MeCN", (3,), {"number": 1.0, "unit": "mL", "raw_text": "1.0 mL"})
    r["steps"] = [step("reaction", (3,), [p["participant_id"] for p in r["participants"]], "stirred under Scheme 3 conditions", {"number": 120, "unit": "°C", "raw_text": "120 °C"}, {"number": 12, "unit": "h", "raw_text": "12 h"}, "air")]
    r["outcomes"] = [outcome("yield", "yield", [prod], (3,), {"number": yld, "unit": "%", "raw_text": f"{yld}%"}, role="desired product")]
    if prod == "4m":
        r["outcomes"].append(outcome("E_Z_ratio", "selectivity", [prod], (3,), {"ratio": [9, 1], "order": ["E", "Z"], "raw_text": "E/Z = 9:1"}, basis="1H NMR", role="E/Z mixture", scope={"metric": "E/Z ratio", "order": ["E", "Z"]}))
    r["general_scheme_ids"] = [r_scheme3["reaction_id"]]

for allene in ("2p", "2q"):
    r = reaction(f"scheme3_{allene}_no_reaction", (3,), categories=("substrate_scope", "control"), location=f"Scheme 3 {allene}, N.R.", class_kind="not_applicable", unresolved=("No product identity or detection limit is reported; no expected product was invented.",))
    r["participants"] = [participant("substrate", "1a", "reactant", (3,), amount={"number": 0.2, "unit": "mmol", "raw_text": "0.2 mmol"}), participant("allene", allene, "reactant", (3,), amount={"number": 0.3, "unit": "mmol", "raw_text": "0.3 mmol"}), participant("cat", "CpRh_dimer", "catalyst", (3,), molpct={"number": 2.5, "unit": "mol%", "raw_text": "2.5 mol%"}, reference="1a"), participant("oxidant", "AgOAc", "reagent", (3,), equiv={"number": 1.0, "unit": "equiv", "raw_text": "1.0 equiv."}, reference="1a"), participant("additive", "NaOAc", "additive", (3,), equiv={"number": 1.0, "unit": "equiv", "raw_text": "1.0 equiv."}, reference="1a"), participant("solvent", "MeCN", "solvent", (3,), volume={"number": 1.0, "unit": "mL", "raw_text": "1.0 mL"})]
    set_scale(r, "1a", {"number": 0.2, "unit": "mmol", "raw_text": "0.2 mmol"}, (3,))
    set_solvent(r, "MeCN", (3,), {"number": 1.0, "unit": "mL", "raw_text": "1.0 mL"})
    r["steps"] = [step("reaction", (3,), [p["participant_id"] for p in r["participants"]], "stirred under Scheme 3 conditions", {"number": 120, "unit": "°C", "raw_text": "120 °C"}, {"number": 12, "unit": "h", "raw_text": "12 h"}, "air")]
    r["outcomes"] = [outcome("no_reaction", "no_reaction", [], (3,), "N.R.", role="no desired annulation product")]
    r["general_scheme_ids"] = [r_scheme3["reaction_id"]]


# Scheme 4a: separate 5 mmol scale-up.
r = reaction("scheme4a_scale_up", (3,), categories=("scale_up",), location="Scheme 4a gram-scale synthesis", unresolved=("MeCN volume and yield measurement basis are not reported in the supplied page.",))
r["participants"] = [participant("substrate", "1a", "reactant", (3,), amount={"number": 5.0, "unit": "mmol", "raw_text": "5.0 mmol"}), participant("allene", "2a", "reactant", (3,), amount={"number": 7.5, "unit": "mmol", "raw_text": "7.5 mmol"}), participant("cat", "CpRh_dimer", "catalyst", (3,), molpct={"number": 2.5, "unit": "mol%", "raw_text": "2.5 mol%"}, reference="1a"), participant("oxidant", "AgOAc", "reagent", (3,), equiv={"number": 1.0, "unit": "equiv", "raw_text": "1.0 equiv."}, reference="1a"), participant("additive", "NaOAc", "additive", (3,), equiv={"number": 1.0, "unit": "equiv", "raw_text": "1.0 equiv."}, reference="1a"), participant("solvent", "MeCN", "solvent", (3,)), participant("product", "3a", "product", (3,), product_role="desired product")]
set_scale(r, "1a", {"number": 5.0, "unit": "mmol", "raw_text": "5.0 mmol"}, (3,))
set_solvent(r, "MeCN", (3,))
r["steps"] = [step("reaction", (3,), [p["participant_id"] for p in r["participants"]], "stirred under Scheme 4a conditions", {"number": 120, "unit": "°C", "raw_text": "120 °C"}, {"number": 12, "unit": "h", "raw_text": "12 h"})]
r["outcomes"] = [outcome("yield", "yield", ["3a"], (3,), {"number": 77, "unit": "%", "raw_text": "77%"}, role="desired product")]
r["general_scheme_ids"] = [r_scheme1["reaction_id"]]

# Scheme 4b: oxidation/derivatization.
r = reaction("scheme4b_derivatization", (3,), categories=("derivatization",), location="Scheme 4b derivatization", name="Oxidation of 3a", transformation="oxidative derivatization", class_kind="oxidation", unresolved=("Reaction temperature, time and yield measurement basis are not reported in the supplied page.", "Exact atom mapping/cleavage coproducts remain unresolved."))
r["participants"] = [participant("substrate", "3a", "reactant", (3,)), participant("cat", "RuCl3", "catalyst", (3,), molpct={"number": 3.0, "unit": "mol%", "raw_text": "3.0 mol%"}, reference="3a"), participant("oxidant", "NaIO4", "reagent", (3,), equiv={"number": 5.0, "unit": "equiv", "raw_text": "5.0 equiv."}, reference="3a"), participant("solvent1", "CH2Cl2", "solvent", (3,)), participant("solvent2", "MeCN", "solvent", (3,)), participant("solvent3", "water", "solvent", (3,)), participant("product", "5", "product", (3,), product_role="derivatization product")]
set_solvent(r, "CH2Cl2/MeCN/H2O", (3,), mixture_ratio="1:1:2")
r["steps"] = [step("reaction", (3,), [p["participant_id"] for p in r["participants"]], "reacted under Scheme 4b conditions")]
r["outcomes"] = [outcome("yield", "yield", ["5"], (3,), {"number": 44, "unit": "%", "raw_text": "44%"}, role="derivatization product")]


def standard_mech_participants(substrates, product_labels, pages=(4,), extra=()):
    items = [participant(f"substrate{i+1}", label, "reactant", pages) for i, label in enumerate(substrates)]
    items += [participant("allene", "2a", "reactant", pages)] if "2a" not in substrates and product_labels else []
    items += [participant("cat", "CpRh_dimer", "catalyst", pages, molpct={"number": 2.5, "unit": "mol%", "raw_text": "2.5 mol%"}, reference=substrates[0]), participant("oxidant", "AgOAc", "reagent", pages, equiv={"number": 1.0, "unit": "equiv", "raw_text": "1.0 equiv."}, reference=substrates[0]), participant("additive", "NaOAc", "additive", pages, equiv={"number": 1.0, "unit": "equiv", "raw_text": "1.0 equiv."}, reference=substrates[0]), participant("solvent", "MeCN", "solvent", pages, volume={"number": 1.0, "unit": "mL", "raw_text": "1.0 mL"})]
    items += list(extra)
    items += [participant(f"product{i+1}", label, "product", pages, product_role="observed product") for i, label in enumerate(product_labels)]
    return items


# Scheme 5a H/D exchange of 1e without allene.
r = reaction("scheme5a_hd_exchange_1e", (3, 4), categories=("mechanistic_study",), location="Scheme 5a, H/D exchange of 1e", name="H/D exchange experiment", transformation="H/D exchange", class_kind="deuteration", unresolved=("The phrase 'standard conditions' does not state whether every standard reagent/solvent amount is unchanged; only explicitly shown CD3OD dose and outcome are certain.", "Measurement method for deuterium incorporation is not stated on supplied pages."))
r["participants"] = [participant("substrate", "1e", "reactant", (4,)), participant("CD3OD", "CD3OD", "reagent", (4,), equiv={"number": 10.0, "unit": "equiv", "raw_text": "10.0 equiv."}, reference="1e"), participant("product", "1e_d2", "product", (4,), product_role="partially deuterated recovered substrate")]
r["steps"] = [step("reaction", (4,), [p["participant_id"] for p in r["participants"]], "treated with CD3OD under standard conditions")]
r["outcomes"] = [outcome("D_incorporation", "selectivity", ["1e_d2"], (4,), {"number": 55, "unit": "%", "raw_text": "55% D; H/D = 0.9:1.1"}, role="partially deuterated material", scope={"metric": "D incorporation", "positions": "ortho-C-H bonds"})]

# Scheme 5a productive deuteration experiment.
r = reaction("scheme5a_deuterio_3a", (3, 4), categories=("mechanistic_study",), location="Scheme 5a, formation of deuterio-3a", name="Deuterium-labelling annulation", class_kind="annulation_deuteration", unresolved=("Measurement method for positional D incorporation is not stated on supplied pages.",))
r["participants"] = standard_mech_participants(["1a"], ["deuterio_3a"], extra=(participant("CD3OD", "CD3OD", "reagent", (4,), equiv={"number": 10.0, "unit": "equiv", "raw_text": "10.0 equiv."}, reference="1a"),))
r["steps"] = [step("reaction", (4,), [p["participant_id"] for p in r["participants"]], "reacted with 2a and CD3OD under standard conditions", {"number": 120, "unit": "°C", "raw_text": "120 °C"}, {"number": 12, "unit": "h", "raw_text": "12 h"}, "air")]
r["outcomes"] = [outcome("methyl_D", "selectivity", ["deuterio_3a"], (4,), {"number": 20, "unit": "%", "raw_text": "H/D (20% D)"}, role="deuterated product", scope={"metric": "D incorporation", "position": "methyl group"}), outcome("ortho_D", "selectivity", ["deuterio_3a"], (4,), {"number": 10, "unit": "%", "raw_text": "D/H (10% D)"}, role="deuterated product", scope={"metric": "D incorporation", "position": "ortho-position of N-aryl ring"})]

# Scheme 5b competition and parallel KIE determinations.
r = reaction("scheme5b_kie_competition", (3, 4), categories=("mechanistic_study", "control"), location="Scheme 5b, competition KIE experiment", name="Kinetic isotope effect competition experiment", unresolved=("How kH/kD was quantified is not stated on supplied pages.",))
r["participants"] = standard_mech_participants(["1a", "1a_d5"], ["3a", "3a_d4"])
r["steps"] = [step("reaction", (4,), [p["participant_id"] for p in r["participants"]], "competition reaction under standard conditions", {"number": 120, "unit": "°C", "raw_text": "standard conditions"}, {"number": 15, "unit": "min", "raw_text": "15 min"}, "air")]
r["outcomes"] = [outcome("KIE", "selectivity", ["3a", "3a_d4"], (4,), {"ratio": [2.3, 1], "order": ["kH", "kD"], "raw_text": "kH/kD = 2.3:1"}, role="comparative rate outcome", scope={"metric": "kinetic isotope effect", "order": ["kH", "kD"]})]

r = reaction("scheme5b_kie_parallel", (3, 4), categories=("mechanistic_study", "control"), location="Scheme 5b, two parallel KIE runs", name="Kinetic isotope effect parallel runs", unresolved=("Two parallel reaction branches are represented in one comparative experiment record; individual conversions/yields and quantification method are not reported.",))
r["participants"] = standard_mech_participants(["1a", "1a_d5"], ["3a", "3a_d4"])
r["steps"] = [step("parallel_H_run", (4,), ["substrate1", "allene", "cat", "oxidant", "additive", "solvent", "product1"], "1a run under standard conditions", {"number": 120, "unit": "°C", "raw_text": "standard conditions"}, {"number": 15, "unit": "min", "raw_text": "15 min"}, "air"), step("parallel_D_run", (4,), ["substrate2", "allene", "cat", "oxidant", "additive", "solvent", "product2"], "1a-d5 run under standard conditions", {"number": 120, "unit": "°C", "raw_text": "standard conditions"}, {"number": 15, "unit": "min", "raw_text": "15 min"}, "air")]
r["outcomes"] = [outcome("KIE", "selectivity", ["3a", "3a_d4"], (4,), {"ratio": [2.3, 1], "order": ["kH", "kD"], "raw_text": "kH/kD = 2.3:1"}, role="comparative rate outcome", scope={"metric": "kinetic isotope effect", "order": ["kH", "kD"]})]

# Scheme 5c electronic competition.
r = reaction("scheme5c_competition", (3, 4), categories=("mechanistic_study", "control"), location="Scheme 5c competition experiment", name="Electronic competition experiment", unresolved=("Measurement basis for the product ratio is not stated.",))
r["participants"] = standard_mech_participants(["1e", "1f"], ["3e", "3f"])
r["steps"] = [step("reaction", (4,), [p["participant_id"] for p in r["participants"]], "equimolar 1e and 1f reacted with 2a under standard conditions", {"number": 120, "unit": "°C", "raw_text": "standard conditions"}, {"number": 12, "unit": "h", "raw_text": "standard conditions"}, "air")]
r["outcomes"] = [outcome("product_ratio", "selectivity", ["3e", "3f"], (4,), {"ratio": [1.3, 1], "order": ["3e", "3f"], "raw_text": "3e:3f = 1.3:1"}, role="competition products", scope={"metric": "product ratio", "order": ["3e", "3f"]})]

# Scheme 5d radical scavengers. The 88% standard row is merged into Table 1 entry 19.
for scavenger, yld in (("TEMPO", 79), ("BHT", 70)):
    r = reaction(f"scheme5d_{scavenger.lower()}", (3, 4), categories=("mechanistic_study", "control"), location=f"Scheme 5d radical inhibition, {scavenger}", unresolved=("Yield measurement basis is not stated in Scheme 5d; it was not copied from Table 1.",))
    r["participants"] = standard_mech_participants(["1a"], ["3a"], extra=(participant("scavenger", scavenger, "additive", (4,), equiv={"number": 3.0, "unit": "equiv", "raw_text": "3.0 equiv."}, reference="1a"),))
    r["steps"] = [step("reaction", (4,), [p["participant_id"] for p in r["participants"]], "reacted under standard conditions with radical scavenger", {"number": 120, "unit": "°C", "raw_text": "120 °C"}, {"number": 12, "unit": "h", "raw_text": "12 h"}, "air")]
    r["outcomes"] = [outcome("yield", "yield", ["3a"], (4,), {"number": yld, "unit": "%", "raw_text": f"{yld}%"}, role="desired product")]


# Scheme 6 is an author-proposed mechanism, not an observed concrete experiment.
r = reaction("scheme6_proposed_mechanism", (4,), "proposed_mechanism", location="Scheme 6 plausible mechanism", name="Plausible Rh(III)-catalyzed annulation mechanism", transformation=None, class_kind="unresolved", unresolved=("Arrow-pushing intermediates I-IV and catalyst regeneration are proposed, not independently observed experiments.",))
r["participants"] = [participant("substrate", "1a", "reactant", (4,)), participant("allene", "2a", "reactant", (4,)), participant("product", "3a", "product", (4,), product_role="mechanistic-cycle product")]
r["notes"] = [f("Authors propose C-H activation, allene insertion, beta-hydride elimination, Rh-H insertion and reductive elimination; AgOAc/NaOAc reoxidize Rh(I).", "author-proposed mechanism", (4,))]


# Link general schemes where applicable.
for r in REACTIONS:
    if r["record_type"] == "specific_experiment" and not r["general_scheme_ids"]:
        if "scheme2_" in r["reaction_id"]:
            r["general_scheme_ids"] = [r_scheme2["reaction_id"]]
        elif "scheme3_" in r["reaction_id"]:
            r["general_scheme_ids"] = [r_scheme3["reaction_id"]]


bundle = st.empty(DATASET["documents"])
bundle["molecules"] = list(MOLECULES.values())
bundle["reactions"] = REACTIONS
bundle["source_context"] = [
    {"context_id": f"{DOC}:title", "kind": "bibliographic", "evidence_ids": [EV[1], TXT[1]], "title": "Direct synthesis of indazole derivatives via Rh(III)-catalyzed C-H activation of phthalazinones and allenes", "doi": "10.1039/d1ob01458g"},
    {"context_id": f"{DOC}:figure1", "kind": "non_synthetic_illustration", "evidence_ids": [EV[1], TXT[1]], "note": "Figure 1 shows pharmaceutical scaffold motifs and is not a synthesis experiment in this paper."},
    {"context_id": f"{DOC}:global_stereochemistry_claim", "kind": "author_summary", "evidence_ids": [EV[1], TXT[1]], "note": "Abstract states high Z-selectivity, but supplied pages give a numeric E/Z order only for 4m; no global per-product configuration was invented."},
    {"context_id": f"{DOC}:esi_notice", "kind": "supplementary_information_notice", "evidence_ids": [EV[1], TXT[1]], "note": "The paper states ESI is available, but no SI file was explicitly supplied or associated in this session."},
    {"context_id": f"{DOC}:page5", "kind": "non_experimental_page", "evidence_ids": [EV[5], TXT[5]], "note": "Page 5 contains conflicts, acknowledgements and references only; no current-paper experiment was promoted from it."},
]

st.write(RUN / "submitted.json", bundle)
st.write(RUN / "raw_model/stereo_diagnostics.json", STEREO_DIAGNOSTICS)

inventory = {
    "document_id": DOC,
    "pages": {
        "1": {"regions": ["title/abstract", "introduction", "Figure 1", "Results opening"], "reaction_ids": []},
        "2": {"regions": ["Scheme 1", "Table 1 entries 1-19", "Scheme 2 products 3a-3w"], "reaction_ids": [r["reaction_id"] for r in REACTIONS if 2 in (st.value(r["location"]).get("pages", []) if isinstance(st.value(r["location"]), dict) else [])]},
        "3": {"regions": ["Scheme 3 products 4a-4n and 2p/2q N.R.", "Scheme 4a-b", "mechanistic discussion opening"], "reaction_ids": [r["reaction_id"] for r in REACTIONS if 3 in (st.value(r["location"]).get("pages", []) if isinstance(st.value(r["location"]), dict) else [])]},
        "4": {"regions": ["Scheme 5a-d", "Scheme 6", "general procedure"], "reaction_ids": [r["reaction_id"] for r in REACTIONS if 4 in (st.value(r["location"]).get("pages", []) if isinstance(st.value(r["location"]), dict) else [])]},
        "5": {"regions": ["conflicts/acknowledgements/references"], "reaction_ids": []},
    },
    "counts": {"molecules": len(MOLECULES), "reactions": len(REACTIONS)},
}
st.write(RUN / "raw_model/page_inventory.json", inventory)

def records_on(page):
    result = []
    for rxn in REACTIONS:
        loc = st.value(rxn["location"])
        if isinstance(loc, dict) and page in loc.get("pages", []):
            result.append(rxn["reaction_id"])
    return result


coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {
            "page": 1,
            "state": "attempted",
            "regions": [
                {"region": "title, bibliographic header and abstract", "finding": "Current-paper title/DOI and author summary; abstract claims moderate-to-good yields and high Z-selectivity but is not a row-level experiment.", "evidence_ids": [EV[1], TXT[1]]},
                {"region": "Introduction and Figure 1", "finding": "Background and pharmaceutical scaffold illustrations; no current-paper synthesis experiment.", "evidence_ids": [EV[1], TXT[1]]},
                {"region": "Results and discussion opening at lower right", "finding": "Model substrates and initial conditions for Table 1 entry 1 begin here and continue on page 2; linked to the single entry-1 record.", "evidence_ids": [EV[1], TXT[1]]},
            ],
            "record_ids": records_on(1),
            "unresolved": ["The ESI availability notice does not constitute an explicitly supplied/associated SI file.", "The entry-1 result sentence continues onto page 2."],
        },
        {
            "page": 2,
            "state": "attempted",
            "regions": [
                {"region": "Scheme 1 at upper left", "finding": "Three cited prior-work arrows were classified as literature examples; the 'This work' arrow was retained as a generic scheme, not a concrete experiment.", "evidence_ids": [EV[2], TXT[2]]},
                {"region": "Table 1 at left", "finding": "All 19 optimization rows were transcribed separately; row 19 is merged with the repeated 3a 88% standard experiment elsewhere.", "evidence_ids": [EV[2], TXT[2]]},
                {"region": "Scheme 2 at lower right", "finding": "General scope arrow plus 3a-3w were inventoried; 3a duplicates Table 1 entry 19, while the other 22 substrate-scope experiments are separate. Regioisomer pairs and ratios were retained.", "evidence_ids": [EV[2], TXT[2]]},
            ],
            "record_ids": records_on(2),
            "unresolved": ["Complete product structures and absolute configuration assignments were not forced into SMILES without the SI.", "Alternative Table 1 oxidant/additive quantities not explicitly stated in the caption remain unreported rather than inferred."],
        },
        {
            "page": 3,
            "state": "attempted",
            "regions": [
                {"region": "upper-left prose continuing Scheme 2", "finding": "Narrative cross-check for products 3b-3w, including regioisomer mixtures; evidence was merged into the corresponding Scheme 2 records rather than duplicated.", "evidence_ids": [EV[3], TXT[3]]},
                {"region": "Scheme 3 at lower left", "finding": "General allene-scope arrow, 14 successful products 4a-4n, and two N.R. substrates 2p/2q were recorded. 4m retains printed E/Z = 9:1 and 1H NMR basis.", "evidence_ids": [EV[3], TXT[3]]},
                {"region": "Scheme 4 at lower right", "finding": "The 5.0 mmol scale-up to 3a and the RuCl3/NaIO4 derivatization to 5 were recorded as distinct experiments.", "evidence_ids": [EV[3], TXT[3]]},
                {"region": "right-column mechanistic-study prose", "finding": "Narrative values for H/D exchange, KIE, electronic competition and radical inhibition were linked to Scheme 5 records on page 4.", "evidence_ids": [EV[3], TXT[3]]},
            ],
            "record_ids": records_on(3),
            "unresolved": ["Scheme 3 does not state the yield measurement basis; those fields remain missing.", "Scheme 4 omits several operational details; no SI details were inferred.", "The printed 4m ratio is preserved in the reported E/Z order even though the abstract globally says high Z-selectivity."],
        },
        {
            "page": 4,
            "state": "attempted",
            "regions": [
                {"region": "Scheme 5a-d at left", "finding": "Two H/D experiments, competition and parallel KIE determinations, electronic competition, and TEMPO/BHT inhibition rows were captured. The uninhibited 88% row is the existing standard 3a experiment and was not duplicated.", "evidence_ids": [EV[4], TXT[4]]},
                {"region": "Scheme 6 at upper right", "finding": "Arrow-pushing cycle classified as a proposed mechanism, not as an experimentally observed sequence.", "evidence_ids": [EV[4], TXT[4]]},
                {"region": "Experimental section at lower right", "finding": "Detailed 3a general procedure, absolute reagent masses/amounts, open-air condition, cooling, evaporation and silica chromatography were merged into Table 1 entry 19.", "evidence_ids": [EV[4], TXT[4]]},
            ],
            "record_ids": records_on(4),
            "unresolved": ["KIE/product-ratio quantification methods are not reported in the supplied pages.", "Scheme 5 radical-inhibition yield bases are unreported.", "Proposed intermediates I-IV were not promoted to observed compounds."],
        },
        {
            "page": 5,
            "state": "attempted",
            "regions": [
                {"region": "full page", "finding": "Conflicts of interest, acknowledgements, and references only; no current-paper experimental arrow, table row, spectrum, mass peak list or synthesis record.", "evidence_ids": [EV[5], TXT[5]]},
            ],
            "record_ids": [],
            "unresolved": [],
        },
    ],
}
st.write(RUN / "coverage/pages_0001_0005.json", coverage)
print(json.dumps(inventory["counts"], ensure_ascii=False))
