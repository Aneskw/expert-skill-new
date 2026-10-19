#!/usr/bin/env python3
"""Build the source-derived v3.1 submission for doc_31c3846e9d6f4c16.

This file is a machine transcription aid. It does not certify source review or
chemical correctness; standard merge will quarantine all records as
needs_review.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

DOC = "doc_31c3846e9d6f4c16"
RUN = ROOT / "runs/current"


def ev(page: int, image: bool = True) -> str:
    return f"{DOC}:page{page}:{'image' if image else 'text'}"


def fx(raw=None, normalized=None, evidence=(), method="original_extraction", derivation=None,
       status=None):
    return st.fact(raw, normalized, evidence, method=method, derivation=derivation,
                   status=status)


def inferred(normalized, evidence, derivation, raw=None):
    return fx(raw, normalized, evidence, method="inference", derivation=derivation)


def derived(normalized, evidence, derivation, raw=None):
    return fx(raw, normalized, evidence, method="structure_derivation", derivation=derivation)


def missing(evidence=()):
    return fx(evidence=evidence)


def molecule(mid, label, evidence, name=None, smiles=None, chirality="unresolved",
             geometry="unresolved", representation="unresolved", assignment=None,
             mixture=None, note=None):
    evidence = list(dict.fromkeys(evidence))
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"] = [{
        "document_id": DOC,
        "compound_label": fx(label, label, evidence),
        "evidence_ids": evidence,
    }]
    m["names"]["original"] = fx(name or label, name or label, evidence)
    m["names"]["english"] = fx(name, name, evidence) if name else missing()
    m["names"]["abbreviations"] = missing()
    m["names"]["aliases"] = missing()
    m["cas_numbers"] = missing()
    if smiles:
        m["representation_type"] = "smiles"
        m["smiles"] = derived(
            smiles, evidence,
            "Candidate connectivity was transcribed from the source drawing; RDKit parsing and stereo diagnostics do not verify the drawing identity.",
            raw=f"candidate connectivity derived from source drawing for {label}",
        )
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = missing()
    m["generic_template"] = None
    stereo = m["stereochemistry"]
    stereo["chirality_status"] = inferred(
        chirality, evidence,
        "Status follows explicit source stereochemical marks/descriptors when present; otherwise uncertainty is preserved rather than interpreted as achiral or racemic.",
        raw=label,
    )
    stereo["geometry_status"] = inferred(
        geometry, evidence,
        "Alkene/imine geometry was retained only when explicitly assigned in the source; otherwise it remains unspecified or unresolved.",
        raw=label,
    )
    stereo["representation_status"] = inferred(
        representation, evidence,
        "The representation-status judgment compares the source stereochemical information with the stored candidate representation.",
        raw=label,
    )
    stereo["atom_centers"] = missing()
    stereo["double_bonds"] = missing()
    stereo["axial"] = missing()
    stereo["mixture"] = fx(mixture, mixture, evidence) if mixture is not None else missing()
    stereo["assignment"] = fx(assignment, assignment, evidence) if assignment is not None else missing()
    stereo["isomeric_smiles"] = missing()
    stereo["identity_key"] = missing()
    if smiles and chirality == "achiral":
        stereo["assignment"] = derived(
            {"assessment": "no tetrahedral or allene-axis stereogenic element in the stored connectivity"},
            evidence,
            "Candidate graph and source drawing were inspected for stereogenic elements; this is a structural inference, not a source descriptor.",
        )
    m["identity_note"] = inferred(
        note or "Identity retained at the source-label/drawing level; no CAS number was reported and no external identity lookup was performed.",
        evidence,
        "Registry note records the limits of the supplied main-text source.",
        raw=label,
    )
    m["external_properties"] = []
    return m


def participant(pid, mid, role, evidence, stage="reaction", amount=None, mass=None,
                volume=None, equiv=None, mol_percent=None, reference=None,
                product_role=None, method="original_extraction"):
    p = st.new_participant(pid, mid, role, evidence, stage=stage)
    p["product_role"] = fx(product_role, product_role, evidence) if product_role else missing()
    values = {
        "amount": amount,
        "mass": mass,
        "volume": volume,
        "equivalents": equiv,
        "mol_percent": mol_percent,
    }
    for key, value in values.items():
        if value is not None:
            p["charges"][key] = fx(value.get("raw_text") or value, value, evidence,
                                      method=method)
    if reference:
        p["equivalence_reference"] = inferred(
            reference, evidence,
            "The loading/equivalent denominator is the hydrazone substrate used as the standard-condition reference.",
            raw="standard-condition substrate basis",
        )
    return p


def step(step_id, evidence, participant_ids, operation, temperature=None, time=None,
         atmosphere=None, pressure=None):
    s = st.new_step(step_id, evidence)
    s["participant_ids"] = list(participant_ids)
    s["operation"] = fx(operation, operation, evidence)
    if temperature is not None:
        s["temperature"] = fx(temperature.get("raw_text") if isinstance(temperature, dict) else temperature,
                                temperature, evidence)
    if time is not None:
        s["time"] = fx(time.get("raw_text") if isinstance(time, dict) else time,
                       time, evidence)
    if atmosphere is not None:
        s["atmosphere"] = fx(atmosphere, atmosphere, evidence)
    if pressure is not None:
        s["pressure"] = fx(pressure.get("raw_text") if isinstance(pressure, dict) else pressure,
                           pressure, evidence)
    return s


def outcome(oid, kind, targets, evidence, value=None, yield_type=None,
            basis=None, scope=None):
    o = st.new_outcome(oid, kind, targets, evidence)
    o["product_role"] = fx("main product", "main product", evidence) if targets else missing()
    if value is not None:
        raw = value.get("raw_text") if isinstance(value, dict) else value
        o["value"] = fx(raw, value, evidence)
    if yield_type is not None:
        o["yield_type"] = fx(yield_type, yield_type, evidence)
    if basis is not None:
        o["measurement_basis"] = fx(basis, basis, evidence)
    if scope is not None:
        o["scope"] = fx(scope, scope, evidence)
    return o


def base_reaction(suffix, evidence, record_type="specific_experiment", categories=None,
                  location=None, name=None, transformation=None):
    evidence = list(dict.fromkeys(evidence))
    r = st.new_reaction(f"{DOC}:{suffix}", DOC, evidence)
    r["record_type"] = record_type
    r["experiment_categories"] = fx(categories, categories, evidence) if categories is not None else missing()
    r["location"] = fx(location, location, evidence) if location is not None else missing()
    r["reaction_name"] = fx(name, name, evidence) if name else missing()
    if transformation:
        r["transformation_type"] = inferred(
            transformation, evidence,
            "Transformation wording is normalized from the source reaction title/scheme and is not a database assignment.",
            raw=name or transformation,
        )
    r["review"] = {
        "status": "needs_review",
        "method": "model",
        "reviewer": None,
        "checked_fields": [],
        "unresolved": [
            "No Supporting Information file was supplied or explicitly associated; detailed procedure, workup, purification, and characterization may be missing.",
        ],
    }
    r["completeness"] = "partial"
    return r


def classify_cc(r, evidence, count=2, note="Two net C-C connections are implied by the depicted [3+2] annulation product."):
    changes = []
    for i in range(count):
        changes.append({
            "change_id": f"cc{i+1}",
            "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}],
            "before_order": 0,
            "after_order": 1,
            "change_type": "formed",
        })
    r["classification"]["bond_changes"] = [derived(c, evidence, note) for c in changes]
    r["classification"]["labels"] = derived(
        [{"category_number": 2, "code": "C-C", "supporting_change_ids": [c["change_id"] for c in changes]}],
        evidence, note,
    )
    r["classification"]["coverage"] = inferred(
        "partial", evidence,
        "Category-level net-bond classification only; dense source drawings were not atom-mapped and hydrogen/fragment provenance is incomplete.",
    )
    r["classification"]["limitations"] = [
        "No atom-complete mapped reaction was asserted.",
        "The two C-C connections are category-level structural deductions from the depicted annulation.",
    ]
    r["classification"]["model_votes"] = []


def classify_hydrogen_gain(r, evidence, element="C", count=2, order_decrease=None):
    changes = []
    if order_decrease:
        a, b = order_decrease
        changes.append({"change_id": "order1", "atoms": [{"atom_map": None, "element": a}, {"atom_map": None, "element": b}],
                        "before_order": 2, "after_order": 1, "change_type": "order_decreased"})
    for i in range(count):
        changes.append({"change_id": f"h{i+1}", "atoms": [{"atom_map": None, "element": element}, {"atom_map": None, "element": "H"}],
                        "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain", "hydrogen_count": 1})
    r["classification"]["bond_changes"] = [derived(c, evidence, "Net reduction shown in Scheme 2; atom maps and hydrogen source are unresolved.") for c in changes]
    code = "C-H" if element == "C" else "N-H"
    number = 1 if element == "C" else 29
    r["classification"]["labels"] = derived(
        [{"category_number": number, "code": code,
          "supporting_change_ids": [c["change_id"] for c in changes if c["change_type"] == "hydrogen_gain"]}],
        evidence, "Only newly gained X-H bonds are classified; reduction of the original double bond is retained as an order decrease, not a new sigma bond.")
    r["classification"]["coverage"] = inferred("partial", evidence, "Hydrogen source/atom mapping is incomplete.")
    r["classification"]["limitations"] = ["Category-level reduction classification; no complete atom mapping."]
    r["classification"]["model_votes"] = []


def classify_none(r, evidence, coverage="not_applicable", reason="No observed synthetic transformation is asserted for this record."):
    r["classification"]["bond_changes"] = []
    r["classification"]["labels"] = inferred([] if coverage == "not_applicable" else None, evidence, reason)
    r["classification"]["coverage"] = inferred(coverage, evidence, reason)
    r["classification"]["limitations"] = [reason]
    r["classification"]["model_votes"] = []


bundle = json.loads((RUN / "submitted.json").read_text())
bundle["molecules"] = []
bundle["reactions"] = []
bundle["source_context"] = []

# Molecule registry: simple candidate graphs are limited to two clearly drawn substrates/solvents.
molecules = {}


def add_mol(*args, **kwargs):
    m = molecule(*args, **kwargs)
    molecules[m["molecule_id"]] = m
    return m["molecule_id"]


M_1A = add_mol(f"{DOC}:mol:1a", "1a", [ev(2), ev(3), ev(4)],
               "acetophenone 2-pyridylhydrazone 1a",
               smiles="CC(=NNc1ccccn1)c1ccccc1", chirality="unspecified",
               geometry="unspecified", representation="partial",
               note="Candidate connectivity for the clearly drawn model hydrazone; C=N geometry is not assigned in the supplied main text.")
M_2A = add_mol(f"{DOC}:mol:2a", "2a", [ev(2), ev(3), ev(4)], "phenylallene 2a",
               smiles="C=C=Cc1ccccc1", chirality="achiral", geometry="not_applicable",
               representation="not_applicable",
               note="Terminal CH2 end prevents allene axial stereogenicity in the stored candidate connectivity.")
M_3A = add_mol(f"{DOC}:mol:3a", "(±)-3a", [ev(2), ev(3), ev(4)], "racemic syn indane 3a",
               chirality="racemic", geometry="not_applicable", representation="not_encoded",
               assignment={"relative_configuration": "syn phenyl/azo", "methods": ["NOE", "single-crystal X-ray diffraction"], "absolute_configuration": None},
               mixture={"type": "racemate", "diastereomer_observation": "single diastereomer in crude 1H NMR"},
               note="Source assigns syn relative configuration and labels the product (±); no absolute descriptor is printed in the supplied main text.")

agents = {
    "co_oac": ("Co(OAc)2·4H2O", "cobalt(II) acetate tetrahydrate"),
    "mn_oac": ("Mn(OAc)2", "manganese(II) acetate"),
    "etoh": ("EtOH", "ethanol"),
    "na_oac": ("NaOAc", "sodium acetate"),
    "na_opiv": ("NaOPiv", "sodium pivalate"),
    "cu_oac": ("Cu(OAc)2", "copper(II) acetate"),
    "pivoh": ("PivOH", "pivalic acid"),
    "co_cl2": ("CoCl2·6H2O", "cobalt(II) chloride hexahydrate"),
    "co_acac": ("Co(acac)2", "cobalt(II) acetylacetonate"),
    "co_no3": ("Co(NO3)2·6H2O", "cobalt(II) nitrate hexahydrate"),
    "co_so4": ("CoSO4·7H2O", "cobalt(II) sulfate heptahydrate"),
    "hfip": ("HFIP", "hexafluoroisopropanol"),
    "tfe": ("TFE", "2,2,2-trifluoroethanol"),
    "tbuoh": ("t-BuOH", "tert-butanol"),
    "meoh": ("MeOH", "methanol"),
    "dce": ("DCE", "1,2-dichloroethane"),
    "acn": ("CH3CN", "acetonitrile"),
    "acetone": ("acetone", "acetone"),
    "fe": ("Fe", "iron powder"),
    "nh4cl": ("NH4Cl", "ammonium chloride"),
    "water": ("H2O", "water"),
    "pdc": ("Pd/C", "palladium on carbon"),
    "h2": ("H2", "hydrogen"),
}
M = {}
for key, (label, name) in agents.items():
    pages = [ev(2), ev(3)] if key in {"co_oac", "mn_oac", "etoh"} else [ev(2)]
    if key in {"fe", "nh4cl", "water", "pdc", "h2", "meoh"}:
        pages = list(dict.fromkeys(pages + [ev(3)]))
    M[key] = add_mol(f"{DOC}:mol:{key}", label, pages, name,
                     chirality="unresolved", geometry="unresolved", representation="unresolved")

# Dense scope identities remain drawing-level records with explicit stereo status.
hydrazone_inputs = {}
for label in [f"1{x}" for x in "bcdefghijklmno"]:
    hydrazone_inputs[label] = add_mol(
        f"{DOC}:mol:{label}", label, [ev(3)], f"hydrazone precursor {label}",
        chirality="unspecified", geometry="unspecified", representation="not_encoded",
        note="Dense Table 2 source drawing retained without a guessed complete SMILES; imine geometry is not explicitly assigned.")

products = {"3a": M_3A}
product_labels = [f"3{x}" for x in list("bcdefghijklmnopqrstuvwxyz")]
product_labels += ["3aa"] + [f"3a{x}" for x in "bcdefghi"]
for label in product_labels:
    geom = "specified" if label in {"3ac", "3ad"} else "not_applicable"
    assignment = {"relative_configuration": "syn aryl/azo", "absolute_configuration": None}
    if label in {"3ac", "3ad"}:
        assignment["alkene_geometry"] = "E"
        assignment["assignment_method"] = "NOE"
    products[label] = add_mol(
        f"{DOC}:mol:{label}", f"(±)-{label}", [ev(3)], f"racemic syn indane {label}",
        chirality="racemic", geometry=geom, representation="not_encoded",
        assignment=assignment, mixture={"type": "racemate"},
        note="Individual dense product drawing retained without a guessed complete SMILES; relative syn stereochemistry follows the table's generic product depiction.")

allene_inputs = {}
for label in [f"for_{p}" for p in [f"3{x}" for x in "pqrstuvwxyz"] + ["3aa"]]:
    prod = label.removeprefix("for_")
    allene_inputs[prod] = add_mol(
        f"{DOC}:mol:allene_{prod}", f"allene precursor to {prod}", [ev(3)],
        f"allene precursor to {prod}", chirality="unresolved", geometry="unresolved",
        representation="not_encoded",
        note="The supplied main-text Table 3 drawing identifies this allene graphically but does not print an unambiguous individual 2-letter label beside every structure.")
allene_inputs["3ab"] = add_mol(f"{DOC}:mol:2m", "2m", [ev(2), ev(3)], "2-naphthyl-substituted allene 2m",
                                chirality="unspecified", geometry="not_applicable", representation="not_encoded")
allene_inputs["3ac"] = add_mol(f"{DOC}:mol:2n", "2n", [ev(2), ev(3)], "disubstituted allene 2n",
                                chirality="unspecified", geometry="not_applicable", representation="not_encoded",
                                note="Potential allene axial stereochemistry is not assigned in the supplied source.")
allene_inputs["3ad"] = add_mol(f"{DOC}:mol:2o", "2o", [ev(2), ev(3)], "disubstituted allene 2o",
                                chirality="unspecified", geometry="not_applicable", representation="not_encoded",
                                note="Potential allene axial stereochemistry is not assigned in the supplied source.")

late_inputs = {}
for prod, bio in [("3ae", "apocynin"), ("3af", "aspirin"), ("3ag", "L-menthol"),
                  ("3ah", "indole"), ("3ai", "morpholine")]:
    late_inputs[prod] = add_mol(
        f"{DOC}:mol:hydrazone_{prod}", f"{bio}-derived hydrazone precursor to {prod}", [ev(3)],
        f"{bio}-derived hydrazone", chirality="unresolved", geometry="unspecified",
        representation="not_encoded",
        note="Bioactive-derived hydrazone retained from the Table 3 drawing; complete connectivity/stereochemistry was not guessed.")

M_4 = add_mol(f"{DOC}:mol:4", "(±)-4", [ev(3)], "hydrazine derivative 4",
              chirality="racemic", geometry="not_applicable", representation="not_encoded",
              assignment={"relative_configuration": "retained indane relationship", "absolute_configuration": None},
              mixture={"type": "racemate"})
M_5 = add_mol(f"{DOC}:mol:5", "(±)-5", [ev(3)], "hydrogenated indane 5",
              chirality="racemic", geometry="not_applicable", representation="not_encoded",
              assignment={"diastereomer_observation": "single diastereomer", "absolute_configuration": None},
              mixture={"type": "racemate"})
M_6 = add_mol(f"{DOC}:mol:6", "6", [ev(3), ev(4)], "N-phenyl hydrazone 6",
              chirality="unspecified", geometry="unspecified", representation="not_encoded")
M_7 = add_mol(f"{DOC}:mol:7", "7", [ev(3), ev(4)], "N-methyl hydrazone 7",
              chirality="unspecified", geometry="unspecified", representation="not_encoded")
M_1F = hydrazone_inputs["1f"]
M_1J = hydrazone_inputs["1j"]
M_D5_1A = add_mol(f"{DOC}:mol:D5_1a", "[D]5-1a", [ev(3), ev(4)], "pentadeuterated hydrazone [D]5-1a",
                  chirality="unspecified", geometry="unspecified", representation="not_encoded")
M_D4_3A = add_mol(f"{DOC}:mol:D4_3a", "[D]4-3a", [ev(4)], "tetradeuterated product [D]4-3a",
                  chirality="racemic", geometry="not_applicable", representation="not_encoded",
                  assignment={"relative_configuration": "syn aryl/azo", "absolute_configuration": None}, mixture={"type": "racemate"})

bundle["molecules"] = list(molecules.values())

# Scheme 1: distinguish cited work from the present generic transformation.
r = base_reaction("literature_scheme1a", [ev(2)], record_type="literature_example",
                  location={"page": 2, "region": "Scheme 1a"},
                  name="Previous work: cobalt-catalyzed annulation with allenes",
                  transformation="cited cobalt/allene annulation overview")
r["notes"] = [fx("Scheme 1a groups prior 4+2, 5+1, 3+3, and 4+1 outcomes around a Co-pi-allyl intermediate.",
                 "Prior-work overview only; not a present-paper experiment.", [ev(2)])]
classify_none(r, [ev(2)], coverage="unresolved", reason="Multiple cited generic transformations are grouped; no single observed present-work bond change is assigned.")
bundle["reactions"].append(r)

r = base_reaction("literature_scheme1b", [ev(2)], record_type="literature_example",
                  location={"page": 2, "region": "Scheme 1b"},
                  name="Previous work: Co(II)-catalyzed [4+2] annulation with 1,3-diynes",
                  transformation="cited [4+2] annulation")
r["notes"] = [fx("Previous work", "Prior-work example only; not counted as a present-paper experiment.", [ev(2)])]
classify_none(r, [ev(2)], coverage="unresolved", reason="Cited generic scheme lacks an atom-complete current-paper experiment.")
bundle["reactions"].append(r)

GEN = f"{DOC}:general_scheme"
r = base_reaction("general_scheme", [ev(1), ev(2)], record_type="generic_scheme",
                  location={"pages": [1, 2], "regions": ["graphical abstract", "Scheme 1c"]},
                  name="Co-catalyzed [3+2] annulation reaction with allenes",
                  transformation="oxidative C-H activation/[3+2] annulation")
r["notes"] = [fx("2,3-migratory insertion followed by nucleophilic addition",
                 "Present-work generic scheme; alternative [4+2] path is depicted but not observed as the reported product class.", [ev(1), ev(2)])]
classify_cc(r, [ev(1), ev(2)])
bundle["reactions"].append(r)


def add_standard_participants(r, evidence, hydrazone=M_1A, allene=M_2A, product=None,
                              h_amount=True, a_amount=True, cobalt="co_oac", cobalt_mol=20,
                              additive="mn_oac", additive_mol=40, solvent="etoh",
                              solvent_volume=2, inherited=False):
    method = "inherited_general_condition" if inherited else "original_extraction"
    hpid = "hydrazone"
    apid = "allene"
    r["participants"].append(participant(hpid, hydrazone, "reactant", evidence,
        amount={"number": 0.2, "unit": "mmol", "raw_text": "0.2 mmol"} if h_amount else None,
        method=method))
    r["participants"].append(participant(apid, allene, "reactant", evidence,
        amount={"number": 0.6, "unit": "mmol", "raw_text": "0.6 mmol"} if a_amount else None,
        method=method))
    if cobalt:
        r["participants"].append(participant("cobalt", M[cobalt], "catalyst", evidence,
            mol_percent={"number": cobalt_mol, "unit": "mol%", "raw_text": f"{cobalt_mol} mol %"}, reference=hydrazone, method=method))
    if additive:
        r["participants"].append(participant("additive", M[additive], "additive", evidence,
            mol_percent={"number": additive_mol, "unit": "mol%", "raw_text": f"{additive_mol} mol %"}, reference=hydrazone, method=method))
    r["participants"].append(participant("solvent", M[solvent], "solvent", evidence,
        volume={"number": solvent_volume, "unit": "mL", "raw_text": f"{solvent_volume} mL"} if solvent_volume else None,
        method=method))
    if product:
        r["participants"].append(participant("product", product, "product", evidence, product_role="main product"))


def standard_step(r, evidence):
    pids = [p["participant_id"] for p in r["participants"]]
    r["steps"] = [step("reaction", evidence, pids, "heat/stir under the indicated optimization conditions",
                       {"number": 70, "unit": "°C", "raw_text": "70 °C"},
                       {"number": 6, "unit": "h", "raw_text": "6 h"}, "air")]


# Table 1. Row 12 is split into seven concrete solvent experiments to avoid representing a false solvent mixture.
opt_specs = [
    ("table1_entry01", 1, {}, [(83, "NMR"), (78, "isolated")], "none"),
    ("table1_entry02", 2, {"additive_mol": 20}, [(50, "NMR")], "20 mol % Mn(OAc)2"),
    ("table1_entry03", 3, {"additive": None}, [(10, "NMR")], "without Mn(OAc)2"),
    ("table1_entry04", 4, {"additive": "na_oac", "additive_mol": 40}, [(25, "NMR")], "NaOAc instead of Mn(OAc)2"),
    ("table1_entry05", 5, {"additive": "na_opiv", "additive_mol": 20}, [(40, "isolated")], "NaOPiv instead of Mn(OAc)2"),
    ("table1_entry06", 6, {"additive": "cu_oac", "additive_mol": 20}, [], "Cu(OAc)2 instead of Mn(OAc)2"),
    ("table1_entry07", 7, {"additive": "pivoh", "additive_mol": 20}, [], "PivOH instead of Mn(OAc)2"),
    ("table1_entry08", 8, {"cobalt": "co_cl2", "additive_mol": 20}, [(20, "NMR")], "with CoCl2·6H2O"),
    ("table1_entry09", 9, {"cobalt": "co_acac", "additive": "na_opiv", "additive_mol": 20}, [(35, "NMR")], "with Co(acac)2"),
    ("table1_entry10", 10, {"cobalt": "co_no3", "additive": "na_opiv", "additive_mol": 20}, [(18, "NMR")], "with Co(NO3)2·6H2O"),
    ("table1_entry11", 11, {"cobalt": "co_so4", "additive": "na_opiv", "additive_mol": 20}, [(19, "NMR")], "with CoSO4·7H2O"),
    ("table1_entry13", 13, {}, [], "at room temperature"),
    ("table1_entry14", 14, {}, [(72, "NMR")], "at 100 °C"),
    ("table1_entry15", 15, {"cobalt": None}, [], "without Co(OAc)2·4H2O"),
]
for suffix, entry, changes, yields, variation in opt_specs:
    evidence = [ev(2)]
    no_product = entry in {6, 7, 13, 15}
    r = base_reaction(suffix, evidence, categories=["optimization"] + (["control"] if entry in {13, 15} else []),
                      location={"page": 2, "region": "Table 1", "entry": entry},
                      name="Optimization of reaction conditions", transformation="[3+2] annulation")
    r["general_scheme_ids"] = [GEN]
    opts = {"product": None if no_product else M_3A}
    opts.update(changes)
    add_standard_participants(r, evidence, **opts)
    standard_step(r, evidence)
    if entry == 13:
        r["steps"][0]["temperature"] = fx("room temperature", "room temperature", evidence)
    if entry == 14:
        r["steps"][0]["temperature"] = fx("100 °C", {"number": 100, "unit": "°C", "raw_text": "100 °C"}, evidence)
    if yields:
        for i, (val, basis) in enumerate(yields, 1):
            r["outcomes"].append(outcome(f"yield{i}", "yield", [M_3A], evidence,
                {"number": val, "unit": "%", "raw_text": f"{val}%"},
                "isolated yield" if basis == "isolated" else "reported yield",
                "isolated after column chromatography" if basis == "isolated" else "1H NMR of crude mixture using 1,3,5-trimethoxybenzene as internal standard"))
        classify_cc(r, evidence)
    elif entry in {13, 15}:
        r["outcomes"] = [outcome("no_reaction", "no_reaction", [], evidence, "no reaction" if entry == 13 else "no product detected")]
        classify_none(r, evidence)
    else:
        r["outcomes"] = [outcome("table_marker", "other", [], evidence, "table dash (–)")]
        classify_none(r, evidence, coverage="unresolved", reason="The table reports a dash and prose only calls the condition deleterious; no quantitative product observation is assigned.")
    r["procedure"]["original_text"] = fx(variation, variation, evidence)
    if entry == 8:
        r["review"]["unresolved"].append("Table 1 entry 8 bears footnote d ('20 mol % additive was used'); the exact intended additive-loading interpretation should be checked against SI.")
    bundle["reactions"].append(r)

for solvent_key, solvent_label in [("hfip", "HFIP"), ("tfe", "TFE"), ("tbuoh", "t-BuOH"),
                                    ("meoh", "MeOH"), ("dce", "DCE"), ("acn", "CH3CN"),
                                    ("acetone", "acetone")]:
    suffix = f"table1_entry12_{solvent_key}"
    r = base_reaction(suffix, [ev(2)], categories=["optimization"],
                      location={"page": 2, "region": "Table 1", "entry": 12, "solvent": solvent_label},
                      name="Optimization of reaction solvent", transformation="[3+2] annulation")
    r["general_scheme_ids"] = [GEN]
    add_standard_participants(r, [ev(2)], product=M_3A, additive="na_opiv", additive_mol=20,
                              solvent=solvent_key)
    standard_step(r, [ev(2)])
    r["outcomes"] = [outcome("yield", "yield", [M_3A], [ev(2)],
                              {"number": 25, "unit": "%", "qualifier": "lt", "raw_text": "<25%"},
                              "reported yield", "1H NMR of crude mixture using 1,3,5-trimethoxybenzene as internal standard")]
    r["notes"] = [fx("HFIP, TFE, t-BuOH, MeOH, DCE, CH3CN, acetone <25%",
                     "Table row 12 aggregates seven separate solvent screens; this record represents one named solvent, not a solvent mixture.", [ev(2)])]
    classify_cc(r, [ev(2)])
    bundle["reactions"].append(r)

# Table 2: 3a is the same standard experiment as Table 1 entry 1 and is not duplicated.
table2 = {
    "3b": 73, "3c": 68, "3d": 71, "3e": 77, "3f": 81, "3g": 80, "3h": 70,
    "3i": 66, "3j": 79, "3k": 81, "3l": 73, "3m": 76, "3n": 75, "3o": 67,
}
for prod, yld in table2.items():
    input_label = "1" + prod[1:]
    hmid = hydrazone_inputs[input_label]
    r = base_reaction(f"table2_{prod}", [ev(3)], categories=["substrate_scope"],
                      location={"page": 3, "region": "Table 2", "product": prod},
                      name="Substrate scope with different hydrazone derivatives",
                      transformation="syn-diastereoselective [3+2] annulation")
    r["general_scheme_ids"] = [GEN]
    add_standard_participants(r, [ev(2), ev(3)], hydrazone=hmid, allene=M_2A, product=products[prod], inherited=True)
    standard_step(r, [ev(3)])
    r["reaction_scale"] = inferred({"reference_molecule_id": hmid, "amount": {"number": 0.2, "unit": "mmol"}},
                                   [ev(2), ev(3)], "Inherited from the optimized standard conditions explicitly invoked for the substrate scope.")
    r["solvent_conditions"]["composition"] = fx("EtOH", [M["etoh"]], [ev(3)])
    r["solvent_conditions"]["total_volume"] = fx("2 mL", {"number": 2, "unit": "mL"}, [ev(3)])
    r["outcomes"] = [outcome("yield", "yield", [products[prod]], [ev(3)],
                              {"number": yld, "unit": "%", "raw_text": f"{yld}%"}, "reported yield")]
    r["review"]["unresolved"].append("The supplied main-text Table 2 does not state the yield measurement basis; it is not assumed to be isolated.")
    classify_cc(r, [ev(3)])
    bundle["reactions"].append(r)

# Add the Table 2 location as evidence to the already represented standard experiment.
entry1 = next(x for x in bundle["reactions"] if x["reaction_id"].endswith("table1_entry01"))
entry1["evidence_ids"].append(ev(3))
entry1["notes"].append(fx("3a, R'=Me, 78%", "Table 2 repeats the isolated 78% standard experiment represented by Table 1 entry 1; no duplicate experiment was created.", [ev(3)]))

# Table 3 allene scope, including bioactive-derived hydrazones at 3ae-3ai.
table3 = {
    "3p": 90, "3q": 79, "3r": 66, "3s": 79, "3t": 78, "3u": 88, "3v": 70,
    "3w": 74, "3x": 81, "3y": 67, "3z": 69, "3aa": 81, "3ab": 95, "3ac": 53,
    "3ad": 59, "3ae": 66, "3af": 52, "3ag": 50, "3ah": 55, "3ai": 58,
}
for prod, yld in table3.items():
    late = prod in late_inputs
    hmid = late_inputs[prod] if late else M_1A
    amid = M_2A if late else allene_inputs[prod]
    r = base_reaction(f"table3_{prod}", [ev(3)], categories=["substrate_scope"],
                      location={"page": 3, "region": "Table 3", "product": prod},
                      name="Substrate scope with different allenes" if not late else "Late-stage modification of bioactive-derived hydrazones",
                      transformation="syn-diastereoselective [3+2] annulation")
    r["general_scheme_ids"] = [GEN]
    add_standard_participants(r, [ev(2), ev(3)], hydrazone=hmid, allene=amid, product=products[prod], inherited=True)
    standard_step(r, [ev(3)])
    r["reaction_scale"] = inferred({"reference_molecule_id": hmid, "amount": {"number": 0.2, "unit": "mmol"}},
                                   [ev(2), ev(3)], "Inherited from the optimized standard conditions explicitly invoked for the scope.")
    basis = "isolated" if prod == "3ae" else None
    r["outcomes"] = [outcome("yield", "yield", [products[prod]], [ev(3)],
                              {"number": yld, "unit": "%", "raw_text": f"{yld}%"},
                              "isolated yield" if basis else "reported yield",
                              "isolated" if basis else None)]
    if not basis:
        r["review"]["unresolved"].append("The supplied main-text Table 3 does not state the yield measurement basis; it is not assumed to be isolated.")
    if prod in {"3ac", "3ad"}:
        r["outcomes"].append(outcome("alkene_geometry", "selectivity", [products[prod]], [ev(2), ev(3)],
            "selective formation of the E-isomer", basis="NOE analysis",
            scope={"metric": "alkene geometry", "major_isomer": "E", "ratio": None}))
    classify_cc(r, [ev(3)])
    bundle["reactions"].append(r)

# Scheme 2 scale-up and two product transformations.
r = base_reaction("scheme2_scaleup_3a", [ev(3)], categories=["scale_up"],
                  location={"page": 3, "region": "Scheme 2a"}, name="Gram-scale synthesis",
                  transformation="syn-diastereoselective [3+2] annulation")
r["general_scheme_ids"] = [GEN]
add_standard_participants(r, [ev(3)], product=M_3A, h_amount=False, a_amount=False)
r["participants"][0]["charges"]["amount"] = fx("4.7 mmol", {"number": 4.7, "unit": "mmol"}, [ev(3)])
r["participants"][0]["charges"]["mass"] = fx("1 g", {"number": 1, "unit": "g"}, [ev(3)])
r["participants"][1]["charges"]["amount"] = fx("14.1 mmol", {"number": 14.1, "unit": "mmol"}, [ev(3)])
r["participants"][1]["charges"]["mass"] = fx("1.6 g", {"number": 1.6, "unit": "g"}, [ev(3)])
standard_step(r, [ev(3)])
r["reaction_scale"] = fx("1a: 4.7 mmol", {"reference_molecule_id": M_1A, "amount": {"number": 4.7, "unit": "mmol"}}, [ev(3)])
r["outcomes"] = [
    outcome("product_amount", "product_amount", [M_3A], [ev(3)], {"number": 1.16, "unit": "g", "raw_text": "1.16 g"}),
    outcome("yield", "yield", [M_3A], [ev(3)], {"number": 75, "unit": "%", "raw_text": "75%"}, "reported yield"),
]
r["review"]["unresolved"].append("The main text does not explicitly state the 75% scale-up yield measurement basis.")
classify_cc(r, [ev(3)])
bundle["reactions"].append(r)

r = base_reaction("scheme2_reduction_to_4", [ev(3)], categories=["derivatization"],
                  location={"page": 3, "region": "Scheme 2b, upper arrow"},
                  name="Selective reduction of the azo group", transformation="azo-to-hydrazine reduction")
r["participants"] = [
    participant("substrate", M_3A, "reactant", [ev(3)]),
    participant("iron", M["fe"], "reagent", [ev(3)], mol_percent={"number": 20, "unit": "mol%", "raw_text": "20 mol %"}, reference=M_3A),
    participant("nh4cl", M["nh4cl"], "reagent", [ev(3)], equiv={"number": 1, "unit": "equiv", "raw_text": "1.0 equiv"}, reference=M_3A),
    participant("water", M["water"], "solvent", [ev(3)]),
    participant("ethanol", M["etoh"], "solvent", [ev(3)]),
    participant("product", M_4, "product", [ev(3)], product_role="main product"),
]
r["solvent_conditions"]["composition"] = fx("H2O/EtOH (1:1)", [M["water"], M["etoh"]], [ev(3)])
r["solvent_conditions"]["mixture_ratio"] = fx("1:1", {"components": [M["water"], M["etoh"]], "ratio": [1, 1]}, [ev(3)])
r["steps"] = [step("reaction", [ev(3)], [p["participant_id"] for p in r["participants"]],
                   "heat from room temperature to reflux", "rt to reflux",
                   {"number": 15, "unit": "min", "raw_text": "15 min"})]
r["outcomes"] = [outcome("yield", "yield", [M_4], [ev(3)], {"number": 98, "unit": "%", "raw_text": "98%"},
                              "isolated yield", "isolated")]
r["review"]["unresolved"].append("The transition/ramp time from room temperature to reflux is not reported separately.")
classify_hydrogen_gain(r, [ev(3)], element="N", count=2, order_decrease=("N", "N"))
bundle["reactions"].append(r)

r = base_reaction("scheme2_hydrogenation_to_5", [ev(3)], categories=["derivatization"],
                  location={"page": 3, "region": "Scheme 2b, lower arrow"},
                  name="Selective hydrogenation of the exo-methylene moiety", transformation="alkene hydrogenation")
r["participants"] = [
    participant("substrate", M_3A, "reactant", [ev(3)]),
    participant("catalyst", M["pdc"], "catalyst", [ev(3)]),
    participant("hydrogen", M["h2"], "reagent", [ev(3)]),
    participant("solvent", M["meoh"], "solvent", [ev(3)]),
    participant("product", M_5, "product", [ev(3)], product_role="main product"),
]
r["solvent_conditions"]["composition"] = fx("MeOH", [M["meoh"]], [ev(3)])
r["steps"] = [step("reaction", [ev(3)], [p["participant_id"] for p in r["participants"]],
                   "Pd/C-catalyzed hydrogenation", "rt", {"number": 17, "unit": "h", "raw_text": "17 h"},
                   "H2", "1 atm")]
r["outcomes"] = [
    outcome("yield", "yield", [M_5], [ev(3)], {"number": 99, "unit": "%", "raw_text": "99%"}, "reported yield"),
    outcome("diastereoselectivity", "selectivity", [M_5], [ev(3)], "single diastereomer",
            scope={"metric": "diastereomeric outcome", "major_isomer": "unassigned", "ratio": None}),
]
r["notes"] = [fx("quantitative yield", "Prose calls the result quantitative while Scheme 2 prints 99%; both source expressions are retained without silently changing the numerical scheme value.", [ev(3)])]
r["review"]["unresolved"].append("Yield measurement basis and selectivity measurement method are not reported in the supplied main text; pressure is preserved as source text '1 atm' because atm is not a supported structured quantity unit in the unchanged schema/tool rules.")
classify_hydrogen_gain(r, [ev(3)], element="C", count=2, order_decrease=("C", "C"))
bundle["reactions"].append(r)

# Scheme 3 mechanistic studies.
for suffix, substrate, statement in [
    ("scheme3_control_6", M_6, "indane product not observed"),
    ("scheme3_control_7", M_7, "no reaction"),
]:
    r = base_reaction(suffix, [ev(3), ev(4)], categories=["control", "mechanistic_study"],
                      location={"page": 4, "region": "Scheme 3a"},
                      name="Directing-group control experiment", transformation="control under optimized annulation conditions")
    add_standard_participants(r, [ev(4)], hydrazone=substrate, allene=M_2A, product=None,
                              h_amount=False, a_amount=False)
    standard_step(r, [ev(4)])
    r["outcomes"] = [outcome("no_reaction", "no_reaction", [], [ev(3), ev(4)], statement)]
    classify_none(r, [ev(4)])
    bundle["reactions"].append(r)

r = base_reaction("scheme3_competition_1f_1j", [ev(4)], categories=["mechanistic_study"],
                  location={"page": 4, "region": "Scheme 3b"}, name="Competitive experiment",
                  transformation="competitive [3+2] annulation")
r["participants"] = [
    participant("hydrazone_1f", M_1F, "reactant", [ev(4)]),
    participant("hydrazone_1j", M_1J, "reactant", [ev(4)]),
    participant("allene", M_2A, "reactant", [ev(4)]),
    participant("cobalt", M["co_oac"], "catalyst", [ev(4)], mol_percent={"number": 20, "unit": "mol%", "raw_text": "20 mol %"}, reference=M_1F),
    participant("additive", M["mn_oac"], "additive", [ev(4)], mol_percent={"number": 40, "unit": "mol%", "raw_text": "40 mol %"}, reference=M_1F),
    participant("solvent", M["etoh"], "solvent", [ev(4)]),
    participant("product_3f", products["3f"], "product", [ev(4)], product_role="competition product"),
    participant("product_3j", products["3j"], "product", [ev(4)], product_role="competition product"),
]
r["notes"] = [fx("1:1 mixture of 4-Me 1f and 4-F 1j", "Reactant mixture ratio 1f:1j = 1:1; absolute charges are not printed in Scheme 3b.", [ev(4)])]
r["solvent_conditions"]["composition"] = fx("EtOH", [M["etoh"]], [ev(4)])
r["solvent_conditions"]["total_volume"] = fx("2 mL", {"number": 2, "unit": "mL"}, [ev(4)])
r["steps"] = [step("reaction", [ev(4)], [p["participant_id"] for p in r["participants"]],
                   "competitive reaction", {"number": 70, "unit": "°C", "raw_text": "70 °C"},
                   {"number": 6, "unit": "h", "raw_text": "6 h"}, "air")]
r["outcomes"] = [outcome("product_ratio", "selectivity", [products["3f"], products["3j"]], [ev(4)],
                              "1.1:1", scope={"metric": "product distribution", "order": ["3f", "3j"], "ratio": [1.1, 1]})]
r["review"]["unresolved"].append("The product-ratio measurement method is not stated in the supplied main text.")
classify_cc(r, [ev(4)], count=4, note="Two [3+2] products are observed in one competition experiment; two category-level C-C connections are assigned per product.")
bundle["reactions"].append(r)

r = base_reaction("scheme3_parallel_kie", [ev(4)], categories=["mechanistic_study"],
                  location={"page": 4, "region": "Scheme 3c"}, name="Parallel kinetic isotope effect",
                  transformation="parallel isotope-effect annulation study")
r["participants"] = [
    participant("h_substrate", M_1A, "reactant", [ev(4)]),
    participant("d_substrate", M_D5_1A, "reactant", [ev(4)]),
    participant("allene", M_2A, "reactant", [ev(4)]),
    participant("h_product", M_3A, "product", [ev(4)], product_role="H-series product"),
    participant("d_product", M_D4_3A, "product", [ev(4)], product_role="D-series product"),
    participant("cobalt", M["co_oac"], "catalyst", [ev(4)], mol_percent={"number": 20, "unit": "mol%", "raw_text": "20 mol %"}, reference=M_1A),
    participant("additive", M["mn_oac"], "additive", [ev(4)], mol_percent={"number": 40, "unit": "mol%", "raw_text": "40 mol %"}, reference=M_1A),
    participant("solvent", M["etoh"], "solvent", [ev(4)]),
]
r["steps"] = [step("parallel_reactions", [ev(4)], [p["participant_id"] for p in r["participants"]],
                   "parallel H/D kinetic experiments", {"number": 70, "unit": "°C", "raw_text": "70 °C"},
                   None, "air")]
r["outcomes"] = [outcome("kie", "other", [M_3A, M_D4_3A], [ev(4)], 1.23,
                              scope={"metric": "kH/kD", "isotope_pair": ["1a", "[D]5-1a"]})]
r["review"]["unresolved"].append("Individual rates, reaction duration, and kinetic measurement method are not reported in the supplied main text.")
classify_cc(r, [ev(4)], count=4, note="Two parallel [3+2] annulations are summarized; category-level C-C formation is retained without atom mapping.")
bundle["reactions"].append(r)

r = base_reaction("scheme3_competitive_kie", [ev(4)], categories=["mechanistic_study"],
                  location={"page": 4, "region": "Scheme 3d"}, name="Competitive kinetic isotope effect",
                  transformation="competitive isotope-effect annulation study")
r["participants"] = [
    participant("h_substrate", M_1A, "reactant", [ev(4)]),
    participant("d_substrate", M_D5_1A, "reactant", [ev(4)]),
    participant("allene", M_2A, "reactant", [ev(4)]),
    participant("h_product", M_3A, "product", [ev(4)], product_role="H-series product"),
    participant("d_product", M_D4_3A, "product", [ev(4)], product_role="D-series product"),
    participant("cobalt", M["co_oac"], "catalyst", [ev(4)], mol_percent={"number": 20, "unit": "mol%", "raw_text": "20 mol %"}, reference=M_1A),
    participant("additive", M["mn_oac"], "additive", [ev(4)], mol_percent={"number": 40, "unit": "mol%", "raw_text": "40 mol %"}, reference=M_1A),
    participant("solvent", M["etoh"], "solvent", [ev(4)]),
]
r["steps"] = [step("competitive_reaction", [ev(4)], [p["participant_id"] for p in r["participants"]],
                   "competitive H/D reaction", {"number": 70, "unit": "°C", "raw_text": "70 °C"},
                   {"number": 1, "unit": "h", "raw_text": "1 h"}, "air")]
r["outcomes"] = [outcome("product_isotope_ratio", "selectivity", [M_3A, M_D4_3A], [ev(4)],
                              "1.56", scope={"metric": "PH/PD", "order": ["3a", "[D]4-3a"], "ratio": [1.56, 1]})]
r["review"]["unresolved"].append("Product-isotope ratio measurement method is not reported in the supplied main text.")
classify_cc(r, [ev(4)], count=4, note="Two isotopologue annulation products are summarized; category-level C-C formation is retained without atom mapping.")
bundle["reactions"].append(r)

# Four prose-reported exchange variants; SI-only procedural details remain absent.
for substrate_mid, substrate_label, direction in [(M_1A, "1a", "H/D"), (M_D5_1A, "[D]5-1a", "D/H")]:
    for with_allene in (False, True):
        suffix = f"exchange_{'d5' if substrate_mid == M_D5_1A else 'h'}_{'with_allene' if with_allene else 'without_allene'}"
        r = base_reaction(suffix, [ev(3), ev(4)], categories=["mechanistic_study"],
                          location={"pages": [3, 4], "region": "prose preceding Scheme 3", "with_allene": with_allene},
                          name=f"{direction} exchange experiment", transformation="isotope-exchange control")
        r["participants"] = [participant("substrate", substrate_mid, "reactant", [ev(3), ev(4)])]
        if with_allene:
            r["participants"].append(participant("allene", M_2A, "reactant", [ev(3), ev(4)]))
        r["steps"] = [step("exchange_test", [ev(3), ev(4)], [p["participant_id"] for p in r["participants"]],
                           f"{direction} exchange experiment; procedural details referred to unavailable SI")]
        r["outcomes"] = [outcome("exchange_result", "other", [substrate_mid], [ev(3), ev(4)], f"no {direction} exchange")]
        r["review"]["unresolved"].append("Exchange reagent/medium, charges, temperature, duration, and analytical method are only referred to SI and are unavailable.")
        classify_none(r, [ev(3), ev(4)])
        bundle["reactions"].append(r)

r = base_reaction("proposed_mechanism_scheme4", [ev(4)], record_type="proposed_mechanism",
                  location={"page": 4, "region": "Scheme 4"}, name="Proposed reaction mechanism",
                  transformation="author-proposed cobalt catalytic cycle")
r["notes"] = [fx("Scheme 4. Proposed Reaction Mechanism",
                 "Author-proposed sequence: coordination, aerial oxidation/CMD, allene coordination and 2,3-migratory insertion, intramolecular nucleophilic addition, elimination/decomplexation, and catalyst reoxidation.", [ev(4)])]
classify_none(r, [ev(4)], reason="Arrow-pushing proposal is not an independently observed net experiment and receives no observed-transformation label.")
bundle["reactions"].append(r)

bundle["source_context"] = [
    {
        "context_id": f"{DOC}:article_metadata",
        "document_id": DOC,
        "evidence_ids": [ev(1), ev(1, False)],
        "title": "Cobalt-Catalyzed C-H Activation and [3 + 2] Annulation with Allenes: Diastereoselective Synthesis of Indane Derivatives",
        "authors": ["Arnab Dey", "Chandra M. R. Volla"],
        "journal": "Organic Letters",
        "year": 2021,
        "volume": 23,
        "pages": "5018-5023",
        "doi": "10.1021/acs.orglett.1c01521",
    },
    {
        "context_id": f"{DOC}:stereo_support",
        "document_id": DOC,
        "evidence_ids": [ev(2), ev(3)],
        "observation": "Crude 1H NMR indicated a single diastereomer for 3a; syn phenyl/azo configuration was assigned by NOE and single-crystal X-ray diffraction. 3w was also analyzed crystallographically; E geometry for 3ac/3ad was assigned by NOE.",
        "ccdc": ["2023190", "2023191", "2023192"],
    },
    {
        "context_id": f"{DOC}:si_notice",
        "document_id": DOC,
        "evidence_ids": [ev(4), ev(6)],
        "observation": "The article states that Supporting Information exists and refers experimental details to it, but no SI file was supplied or explicitly associated in this workspace.",
        "association_status": "unestablished",
    },
    {
        "context_id": f"{DOC}:spectra_triage",
        "document_id": DOC,
        "evidence_ids": [ev(1), ev(2), ev(3), ev(4), ev(5), ev(6)],
        "observation": "No standalone mass-spectrum peak table, molecular-weight table, or fragment-relative-abundance table appears on the six supplied main-text pages.",
    },
]

st.write(RUN / "submitted.json", bundle)

# Raw transcription inventory and deterministic stereo diagnostics.
inventory = {
    "document_id": DOC,
    "source_sha256": bundle["documents"][0]["sha256"],
    "title": bundle["documents"][0]["title"],
    "si_association": "unestablished",
    "table1": {
        "rows": list(range(1, 16)),
        "row12_split_solvents": ["HFIP", "TFE", "t-BuOH", "MeOH", "DCE", "CH3CN", "acetone"],
        "yield_cells": {"1": "83 (78 isolated)", "2": "50", "3": "10", "4": "25", "5": "40 isolated", "6": "–", "7": "–", "8": "20", "9": "35", "10": "18", "11": "19", "12": "<25%", "13": "– / prose no reaction", "14": "72", "15": "– / prose no product"},
    },
    "table2": table2,
    "table2_duplicate": {"3a": "same standard run as Table 1 entry 1; isolated 78%"},
    "table3": table3,
    "scheme2": ["gram-scale 3a", "3a to 4", "3a to 5"],
    "scheme3": ["control 6", "control 7", "competition 1f/1j", "parallel KIE", "competitive KIE"],
    "prose_exchange_variants": 4,
    "proposed_mechanism": "Scheme 4",
    "pages_visually_attempted": [1, 2, 3, 4, 5, 6],
    "record_count": len(bundle["reactions"]),
    "molecule_count": len(bundle["molecules"]),
}
st.write(RUN / "raw_model/extraction_inventory.json", inventory)

stereo = {}
for label, smiles in {
    "1a_candidate": "CC(=NNc1ccccn1)c1ccccc1",
    "2a_candidate": "C=C=Cc1ccccc1",
    "disubstituted_allene_connectivity_example_not_identity_assignment": "CC=C=Cc1ccccc1",
}.items():
    stereo[label] = st.stereo_inspect(smiles)
st.write(RUN / "raw_model/stereo_diagnostics.json", stereo)

# Coverage maps every visually inspected page and does not claim semantic verification.
ids = [r["reaction_id"] for r in bundle["reactions"]]
def ids_matching(*terms):
    return [rid for rid in ids if any(t in rid for t in terms)]

coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {"page": 1, "state": "attempted", "regions": [
            {"region": "title/abstract and graphical abstract", "finding": "Article metadata and present-work generic hydrazone/allene [3+2] annulation; graphical summary is duplicate evidence for the generic scheme, not a separate experiment.", "evidence_ids": [ev(1), ev(1, False)]},
            {"region": "body text columns", "finding": "Background/cited work and rationale; cited reactions were not imported as present-paper experiments.", "evidence_ids": [ev(1), ev(1, False)]},
        ], "record_ids": [GEN], "unresolved": ["Dense graphical-abstract structures were not converted to atom-complete verified SMILES."]},
        {"page": 2, "state": "attempted", "regions": [
            {"region": "Scheme 1a", "finding": "Previous-work overview grouped as a literature example.", "evidence_ids": [ev(2)]},
            {"region": "Scheme 1b", "finding": "Previous [4+2] diene work grouped as a literature example.", "evidence_ids": [ev(2)]},
            {"region": "Scheme 1c", "finding": "Present generic [3+2] scheme; alternative [4+2] path is proposed, not reported as the product outcome.", "evidence_ids": [ev(2)]},
            {"region": "Table 1 entries 1-15 and footnotes", "finding": "All rows represented; aggregated row 12 split into seven solvent experiments. Entry 1 includes NMR 83% and isolated 78%; entry 5 is isolated 40%; entries 13/15 retain explicit no-reaction/no-product prose.", "evidence_ids": [ev(2)]},
            {"region": "scope discussion", "finding": "Cross-page cues for 3n/3o, 3ab, and E-3ac/E-3ad mapped to page-3 scope records.", "evidence_ids": [ev(2), ev(2, False)]},
        ], "record_ids": ids_matching("literature_scheme1", "general_scheme", "table1_", "table3_3ab", "table3_3ac", "table3_3ad"), "unresolved": ["Entry 8 footnote-d additive interpretation requires unavailable SI.", "Table dashes for entries 6/7 are preserved rather than upgraded to no reaction."]},
        {"page": 3, "state": "attempted", "regions": [
            {"region": "Table 2", "finding": "Products 3a-3o inventoried. 3a is duplicate evidence for Table 1 entry 1; 3b-3o are independent scope records.", "evidence_ids": [ev(3)]},
            {"region": "Table 3", "finding": "Products 3p-3ai inventoried as independent scope records; E assignments for 3ac/3ad retained with NOE basis.", "evidence_ids": [ev(3)]},
            {"region": "Scheme 2a", "finding": "Independent gram-scale experiment to 3a.", "evidence_ids": [ev(3)]},
            {"region": "Scheme 2b upper/lower arrows", "finding": "Two independent derivatizations to 4 and 5; rt-to-reflux ramp remains un-timed; scheme 99% and prose 'quantitative' for 5 are both retained.", "evidence_ids": [ev(3)]},
            {"region": "mechanistic prose below Scheme 2", "finding": "Two directing-group controls and four SI-referred isotope-exchange variants identified; control schemes continue on page 4.", "evidence_ids": [ev(3), ev(3, False)]},
        ], "record_ids": ids_matching("table1_entry01", "table2_", "table3_", "scheme2_", "scheme3_control", "exchange_"), "unresolved": ["Most scope yield measurement bases are not stated.", "Complete dense substrate/product SMILES and absolute configurations are not safely recoverable from the main-text page alone.", "SI-only exchange conditions are unavailable."]},
        {"page": 4, "state": "attempted", "regions": [
            {"region": "Scheme 3a", "finding": "Two real no-product/no-reaction directing-group control experiments.", "evidence_ids": [ev(4)]},
            {"region": "Scheme 3b", "finding": "One competition experiment; 3f:3j = 1.1:1.", "evidence_ids": [ev(4)]},
            {"region": "Scheme 3c", "finding": "Parallel KIE experiment, kH/kD = 1.23.", "evidence_ids": [ev(4)]},
            {"region": "Scheme 3d", "finding": "Competitive isotope experiment, PH/PD = 1.56 after 1 h.", "evidence_ids": [ev(4)]},
            {"region": "Scheme 4 and explanatory prose", "finding": "Author-proposed catalytic cycle kept separate as proposed_mechanism, not experimental proof.", "evidence_ids": [ev(4), ev(4, False)]},
            {"region": "Associated Content", "finding": "SI is stated to exist, but no SI file is supplied or associated; CCDC 2023190-2023192 noted as source context.", "evidence_ids": [ev(4)]},
        ], "record_ids": ids_matching("scheme3_", "exchange_", "proposed_mechanism"), "unresolved": ["KIE raw rates and analytical methods are unreported in the supplied main text.", "Mechanistic intermediates were not treated as isolated experimental products."]},
        {"page": 5, "state": "attempted", "regions": [
            {"region": "notes/acknowledgments", "finding": "No synthetic experiment.", "evidence_ids": [ev(5), ev(5, False)]},
            {"region": "references 1-13", "finding": "Cited literature only; not imported as present-paper experiments.", "evidence_ids": [ev(5), ev(5, False)]},
        ], "record_ids": [], "unresolved": []},
        {"page": 6, "state": "attempted", "regions": [
            {"region": "references 13-16 and SI cross-reference", "finding": "Cited literature and a statement to see SI; no present-paper experiment or independent spectral table.", "evidence_ids": [ev(6), ev(6, False)]},
        ], "record_ids": [], "unresolved": ["The cited SI schematic/optimization details are unavailable because no SI file was supplied."]},
    ],
}
st.write(RUN / "coverage/pages_0001_0006.json", coverage)

tasks = json.loads((RUN / "tasks.json").read_text())
for task in tasks:
    task["state"] = "attempted"
    task["reason"] = "Original page image and extracted text were inspected; semantic accuracy remains unverified. See coverage/pages_0001_0006.json."
st.write(RUN / "tasks.json", tasks)

print(json.dumps({"molecules": len(bundle["molecules"]), "reactions": len(bundle["reactions"])}, ensure_ascii=False))
