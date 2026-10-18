#!/usr/bin/env python3
"""Build the source-specific v3.1 bundle for doc_151fad1677f63120.

This is a deterministic transcription aid for the visually inspected paper pages.
It does not claim chemical or semantic verification.
"""
from __future__ import annotations

import copy
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

RUN = ROOT / "runs/current"
DATASET = RUN / "dataset.json"
OUT = RUN / "submitted.json"
DOC = "doc_151fad1677f63120"
E = {p: f"{DOC}:page{p}:image" for p in range(1, 6)}


def f(raw=None, norm=None, ev=(), method="original_extraction", status=None, derivation=None):
    return st.fact(raw, norm, ev, method=method, status=status, derivation=derivation)


def q(number, unit, qualifier=None, raw_text=None):
    value = {"number": number, "unit": unit}
    if qualifier:
        value["qualifier"] = qualifier
    if raw_text:
        value["raw_text"] = raw_text
    return value


def mid(key):
    return f"{DOC}:mol:{key}"


bundle = st.empty(json.loads(DATASET.read_text())["documents"])
molecules = {}
stereo_diagnostics = {}


def molecule(key, label, name=None, ev=(E[1],), smiles=None, chirality="unresolved",
             geometry="unresolved", representation=None, assignment=None,
             axial=None, atom_centers=None, double_bonds=None, mixture=None,
             identity_note=None):
    ident = mid(key)
    m = st.new_molecule(ident, DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f(label, label, ev)
    m["occurrences"][0]["evidence_ids"] = list(ev)
    original = name or label
    m["names"]["original"] = f(original, original, ev)
    m["names"]["english"] = f(original, original, ev)
    m["cas_numbers"] = f()
    if smiles:
        m["representation_type"] = "smiles"
        m["smiles"] = f(
            f"candidate transcription for {label}", smiles, ev,
            method="structure_derivation",
            derivation="Candidate SMILES transcribed from the source drawing/name; graph parsing is not source identity verification.",
        )
        diag = st.stereo_inspect(smiles)
        stereo_diagnostics[ident] = diag
        m["stereochemistry"]["isomeric_smiles"] = f(
            None, diag["canonical_isomeric_smiles"], ev, method="structure_derivation",
            derivation="RDKit canonical isomeric SMILES for the unverified candidate graph.",
        )
        m["stereochemistry"]["identity_key"] = f(
            None, diag["stereo_graph_key"], ev, method="structure_derivation",
            derivation="Stereo graph key calculated from the unverified candidate SMILES.",
        )
        if representation is None:
            representation = "partial" if any(x["specified"] == "Unspecified" for x in diag["potential_stereo"]) else "encoded"
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = f()
        representation = representation or "not_encoded"
    m["stereochemistry"]["chirality_status"] = f(chirality, chirality, ev)
    m["stereochemistry"]["geometry_status"] = f(geometry, geometry, ev)
    m["stereochemistry"]["representation_status"] = f(representation, representation, ev)
    if assignment is not None:
        m["stereochemistry"]["assignment"] = f(assignment, assignment, ev)
    if axial is not None:
        m["stereochemistry"]["axial"] = f(axial, axial, ev)
    if atom_centers is not None:
        m["stereochemistry"]["atom_centers"] = f(atom_centers, atom_centers, ev)
    if double_bonds is not None:
        m["stereochemistry"]["double_bonds"] = f(double_bonds, double_bonds, ev)
    if mixture is not None:
        m["stereochemistry"]["mixture"] = f(mixture, mixture, ev)
    if identity_note:
        m["identity_note"] = f(identity_note, identity_note, ev)
    molecules[ident] = m
    return ident


def generic_molecule(key, label, description, ev=(E[1],), variable_labels=("R",)):
    ident = mid(key)
    m = st.new_molecule(ident, DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f(label, label, ev)
    m["occurrences"][0]["evidence_ids"] = list(ev)
    m["names"]["original"] = f(label, label, ev)
    m["names"]["english"] = f(description, description, ev)
    m["cas_numbers"] = f()
    m["representation_type"] = "generic_template"
    m["smiles"] = f()
    variables = []
    for variable in variable_labels:
        variables.append({
            "label": variable,
            "attachment_points": f(f"{variable} attachment shown in source drawing", {"label": variable, "position": "source-defined"}, ev),
            "definition": f(f"{variable} as drawn", "source-defined substituent", ev),
            "allowed_values": f(),
            "constraints": f(),
        })
    m["generic_template"] = {
        "fixed_structure": f(description, {"description": description}, ev),
        "variables": variables,
        "instance_mapping": [],
        "reason": f("R-group or abbreviated source drawing", "A concrete single-molecule SMILES cannot represent the source family safely.", ev),
    }
    for key2, value in (("chirality_status", "unresolved"), ("geometry_status", "unresolved"),
                        ("representation_status", "not_encoded")):
        m["stereochemistry"][key2] = f(value, value, ev)
    m["identity_note"] = f(description, description, ev)
    molecules[ident] = m
    return ident


# Generic source families.
M1 = generic_molecule("family_enallene_1", "1", "R-substituted enallenol family 1 shown in Schemes 1-3", (E[1], E[3], E[4]))
M2 = generic_molecule("family_product_2", "2", "alkynyl-substituted cyclopentenol family 2 with source-drawn cis relationship", (E[1], E[3], E[4]), ("R", "R'"))
M3 = generic_molecule("family_product_3", "3", "oxygen-heterocycle family 3", (E[1], E[4]))
M4 = generic_molecule("family_product_4", "4", "alkynylated open-chain family 4", (E[1], E[4]), ("R", "R'"))
M5 = generic_molecule("family_product_5", "5", "cyclopentenol family 5", (E[1], E[4]), ("R",))
MALK = generic_molecule("family_terminal_alkyne", "R'-C≡CH", "terminal alkyne family shown with variable R'", (E[1], E[3], E[4]), ("R'",))
MINT = generic_molecule("family_proposed_intermediates", "Int-1/Int-2/Int-3", "author-proposed palladium intermediates in Scheme 5", (E[4],), ("R", "R'", "X"))

# Common concrete/candidate molecules and reaction agents.
M1A = molecule(
    "1a", "1a", "enallenol 1a", (E[1], E[2], E[3]),
    smiles="CC(C)=C=C(C1CCCCC1)CC(O)CC=C", chirality="unspecified",
    geometry="unspecified", representation="partial",
    atom_centers=[{"feature": "carbinol carbon", "configuration": "unspecified"}],
    axial=[{"feature": "allene axis", "status": "not assigned; candidate graph has two methyl groups at one terminus"}],
    identity_note="Candidate connectivity is a visual transcription and remains unverified; the alcohol center is not assigned in the source pages.",
)
MPH = molecule("phenylacetylene", "PhC≡CH", "phenylacetylene", (E[1], E[2], E[4]), smiles="C#Cc1ccccc1", chirality="achiral", geometry="not_applicable")
MBQ = molecule("BQ", "BQ", "1,4-benzoquinone", (E[1], E[2], E[3], E[4]), smiles="O=C1C=CC(=O)C=C1", chirality="achiral", geometry="not_applicable")
MPDOAC = molecule("PdOAc2", "Pd(OAc)2", "palladium(II) acetate", (E[1], E[2], E[3], E[4]), smiles="CC(=O)[O-].CC(=O)[O-].[Pd+2]", chirality="not_applicable", geometry="not_applicable")
MPDTFA = molecule("PdTFA2", "Pd(TFA)2", "palladium(II) trifluoroacetate", (E[2],), smiles="O=C([O-])C(F)(F)F.O=C([O-])C(F)(F)F.[Pd+2]", chirality="not_applicable", geometry="not_applicable")
MPDPH3 = generic_molecule("PdPPh3Cl2", "Pd(PPh3)2Cl2", "bis(triphenylphosphine)palladium(II) dichloride", (E[2],), ())
MPDMCF = generic_molecule("Pd_AmP_MCF", "Pd-AmP-MCF", "palladium immobilized on amino-functionalized siliceous mesocellular foam", (E[2], E[3], E[4]), ())
MPDCNC = generic_molecule("Pd_AmP_CNC", "Pd-AmP-CNC", "palladium immobilized on amino-functionalized crystalline nanocellulose foam", (E[2],), ())
MAGOTF = molecule("AgOTf", "AgOTf", "silver trifluoromethanesulfonate", (E[2], E[4], E[5]), smiles="O=S(=O)([O-])C(F)(F)F.[Ag+]", chirality="not_applicable", geometry="not_applicable")
MET3N = molecule("Et3N", "Et3N", "triethylamine", (E[2], E[3], E[4]), smiles="CCN(CC)CC", chirality="achiral", geometry="not_applicable")
MPR3N = molecule("Pr3N", "Pr3N", "tri-n-propylamine", (E[2],), smiles="CCCN(CCC)CCC", chirality="achiral", geometry="not_applicable")
MPR2NH = molecule("Pr2NH", "Pr2NH", "di-n-propylamine", (E[2],), smiles="CCCNCCC", chirality="achiral", geometry="not_applicable")
MPRNH2 = molecule("PrNH2", "PrNH2", "n-propylamine", (E[2],), smiles="CCCN", chirality="achiral", geometry="not_applicable")
MTMEDA = molecule("TMEDA", "TMEDA", "N,N,N',N'-tetramethylethylenediamine", (E[2],), smiles="CN(C)CCN(C)C", chirality="achiral", geometry="not_applicable")
MSOL = {
    "DCE": molecule("DCE", "DCE", "1,2-dichloroethane", (E[1], E[2], E[3], E[4]), smiles="ClCCCl", chirality="achiral", geometry="not_applicable"),
    "THF": molecule("THF", "THF", "tetrahydrofuran", (E[2],), smiles="C1CCOC1", chirality="achiral", geometry="not_applicable"),
    "Toluene": molecule("toluene", "Toluene", "toluene", (E[2],), smiles="Cc1ccccc1", chirality="achiral", geometry="not_applicable"),
    "CHCl3": molecule("chloroform", "CHCl3", "chloroform", (E[2], E[4], E[5]), smiles="ClC(Cl)Cl", chirality="achiral", geometry="not_applicable"),
    "Acetone": molecule("acetone", "Acetone", "acetone", (E[2],), smiles="CC(C)=O", chirality="achiral", geometry="not_applicable"),
    "CH3CN": molecule("acetonitrile", "CH3CN", "acetonitrile", (E[2], E[4]), smiles="CC#N", chirality="achiral", geometry="not_applicable"),
}
MB2PIN2 = molecule("B2pin2", "B2pin2", "bis(pinacolato)diboron", (E[4],), smiles="CC1(C)OB(B2OC(C)(C)C(C)(C)O2)OC1(C)C", chirality="achiral", geometry="not_applicable")
MANISOLE = molecule("anisole", "anisole", "anisole (NMR internal standard)", (E[2],), smiles="COc1ccccc1", chirality="achiral", geometry="not_applicable")


def product_molecule(series, label, ev, stereo_kind="unspecified", assignment=None):
    key = f"product_{label.replace("'", "prime").replace(' ', '_')}"
    return molecule(key, label, f"source-drawn product {label}", ev, smiles=None,
                    chirality=stereo_kind, geometry="not_applicable", representation="not_encoded",
                    assignment=assignment,
                    identity_note="The compact source drawing was not converted to a complete machine SMILES; label and stereochemical statements are retained.")


# Recurrent products.
P2 = {}
for label in [f"2{c}" for c in "abcdefghijklmnopqr"]:
    P2[label] = product_molecule("2", label, (E[3],), "partially_specified", "cis relative configuration is shown/reported; absolute configuration is unassigned")
P3A = product_molecule("3", "3a", (E[2], E[4]), "unspecified")
P4A = product_molecule("4", "4a", (E[2], E[4]), "unspecified")
P5A = product_molecule("5", "5a", (E[1], E[2], E[4], E[5]), "unspecified")

# Source labels for Scheme 3 products not already registered.
SCHEME3_PRODUCTS = {"3a": P3A, "4a": P4A, "5a": P5A}
for series in ("3", "4", "5"):
    for letter in "bcd":
        label = series + letter
        SCHEME3_PRODUCTS[label] = product_molecule(series, label, (E[4], E[5]) if series == "5" else (E[4],), "unspecified")


def alkyne_molecule(key, label, smiles, ev=(E[3],)):
    return molecule(key, label, label, ev, smiles=smiles, chirality="achiral", geometry="not_applicable")


MALK_2B = alkyne_molecule("alkyne_2b", "4-ethynylanisole", "COc1ccc(C#C)cc1")
ALKYNES = {
    "2a": MPH,
    "2b": MALK_2B,
    "2c": alkyne_molecule("alkyne_2c", "2-ethynylanisole", "COc1ccccc1C#C"),
    "2d": alkyne_molecule("alkyne_2d", "3-ethynylanisole", "COc1cccc(C#C)c1"),
    "2e": alkyne_molecule("alkyne_2e", "1-ethynyl-4-nitrobenzene", "C#Cc1ccc([N+](=O)[O-])cc1"),
    "2f": alkyne_molecule("alkyne_2f", "1-bromo-4-ethynylbenzene", "C#Cc1ccc(Br)cc1"),
    "2g": alkyne_molecule("alkyne_2g", "1-chloro-4-ethynylbenzene", "C#Cc1ccc(Cl)cc1"),
    "2h": alkyne_molecule("alkyne_2h", "methyl 4-ethynylbenzoate", "COC(=O)c1ccc(C#C)cc1"),
    "2i": alkyne_molecule("alkyne_2i", "1-ethynyl-4-fluorobenzene", "C#Cc1ccc(F)cc1"),
    "2j": alkyne_molecule("alkyne_2j", "1-ethynyl-4-(trifluoromethyl)benzene", "C#Cc1ccc(C(F)(F)F)cc1"),
    "2k": generic_molecule("alkyne_2k", "ethynylthiophene isomer for 2k", "source-drawn ethynylthiophene for product 2k; attachment isomer not asserted", (E[3],), ()),
    "2l": generic_molecule("alkyne_2l", "ethynylthiophene isomer for 2l", "source-drawn ethynylthiophene for product 2l; attachment isomer not asserted", (E[3],), ()),
    "2m": MALK_2B,
    "2n": MPH,
    "2o": MPH,
    "2p": MPH,
    "2q": MALK_2B,
    "2r": alkyne_molecule("alkyne_2r", "4-ethynyl-N,N-dimethylaniline", "C#Cc1ccc(N(C)C)cc1"),
}

SCOPE_SUBSTRATES = {label: M1A for label in P2}
for label, desc in {
    "2m": "R=propyl enallenol 1", "2n": "R=butyl enallenol 1", "2o": "R=cyclopropyl enallenol 1",
    "2p": "R=cyclopentyl enallenol 1", "2q": "R=phenylethyl enallenol 1", "2r": "tertiary-amine-containing enallenol 1",
}.items():
    SCOPE_SUBSTRATES[label] = molecule(f"substrate_{label}", desc, desc, (E[3],), smiles=None,
                                       chirality="unspecified", geometry="unspecified", representation="not_encoded",
                                       identity_note="Specific compact substrate drawing not converted to SMILES.")


def classification(ev, codes=(), coverage="partial", limitation=None):
    c = st._skeleton("classification", ev)
    c["standard_id"] = st.STANDARD
    c["mapped_reaction_smiles"] = f()
    if coverage in ("unresolved", "not_applicable"):
        c["labels"] = f(None, [], ev, method="inference", derivation="No observed bond-formation labels asserted for this record.")
        c["coverage"] = f(None, coverage, ev, method="inference", derivation="Source outcome does not support a mapped product transformation.")
        c["bond_changes"] = []
    else:
        changes = []
        labels = []
        for index, code in enumerate(codes, 1):
            order = 2 if code.startswith("DB-") else 3 if code.startswith("TB-") else 1
            pair = code.removeprefix("DB-").removeprefix("TB-").split("-")
            if code == "C-B": pair = ["C", "B"]
            if code == "C-H": pair = ["C", "H"]
            change_type = "order_increased" if order > 1 else ("hydrogen_gain" if code == "C-H" else "formed")
            before = 1 if order > 1 else 0
            change = {"change_id": f"b{index}", "atoms": [{"element": pair[0]}, {"element": pair[1]}],
                      "before_order": before, "after_order": order, "change_type": change_type}
            if code == "C-H": change["hydrogen_count"] = 1
            changes.append(f(None, change, ev, method="inference", derivation="Class-level net bond change inferred from the source transformation drawing; atom identities are unresolved."))
            row = st.bond_label(pair[0], pair[1], order)
            labels.append({**row, "supporting_change_ids": [f"b{index}"]})
        c["bond_changes"] = changes
        c["labels"] = f(None, labels, ev, method="inference", derivation="Formation labels inferred at element/bond-order class level from the source drawing.")
        c["coverage"] = f(None, "partial", ev, method="inference", derivation="No full atom mapping; coproducts, exact atom identities, and complete hydrogen changes remain unresolved.")
    c["limitations"] = [limitation or "Class-level inference only; no complete atom mapping or complete hydrogen/coproduct accounting."]
    c["model_votes"] = []
    return c


reactions = []


def reaction(suffix, ev, record_type="specific_experiment", categories=(), name=None, transform=None,
             location=None, class_codes=(), class_coverage="partial", unresolved=()):
    rid = f"{DOC}:{suffix}"
    r = st.new_reaction(rid, DOC, ev)
    r["record_type"] = record_type
    r["experiment_categories"] = f(", ".join(categories) if categories else "not applicable", list(categories), ev)
    r["location"] = f(location or suffix, {"source_location": location or suffix}, ev)
    r["reaction_name"] = f(name, name, ev) if name else f()
    if transform:
        r["transformation_type"] = f(transform, transform, ev)
    r["classification"] = classification(ev, class_codes, class_coverage)
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None, "checked_fields": [],
                   "unresolved": list(unresolved) or ["Machine record requires expert source review; no SI was supplied."]}
    r["completeness"] = "partial"
    reactions.append(r)
    return r


def add_participant(r, molecule_id, role, ev, stage="reaction", product_role=None, charges=None, reference=None):
    pid = f"p{len(r['participants']) + 1}"
    p = st.new_participant(pid, molecule_id, role, ev, stage)
    if product_role:
        p["product_role"] = f(product_role, product_role, ev)
    for key, value in (charges or {}).items():
        raw, norm = value
        p["charges"][key] = f(raw, norm, ev)
    if any(p["charges"][k]["normalized_value"] is not None for k in ("equivalents", "mol_percent")):
        if reference is None:
            raise ValueError(f"reference required for {molecule_id}")
        p["equivalence_reference"] = f(f"relative to {reference}", reference, ev)
    r["participants"].append(p)
    return p


def add_step(r, ev, operation="reaction", temp=None, time=None, atmosphere=None, stage="reaction", participants=None):
    s = st.new_step(f"step_{len(r['steps']) + 1}", ev, stage)
    s["participant_ids"] = list(participants or [p["participant_id"] for p in r["participants"] if p["stage"] == stage])
    s["operation"] = f(operation, operation, ev)
    if temp is not None:
        s["temperature"] = f(temp[0], temp[1], ev)
    if time is not None:
        s["time"] = f(time[0], time[1], ev)
    if atmosphere is not None:
        s["atmosphere"] = f(atmosphere, atmosphere, ev)
    r["steps"].append(s)
    return s


def add_outcome(r, kind, target_mids, ev, raw, norm, basis=None, yield_type=None, scope=None, product_role=None):
    o = st.new_outcome(f"outcome_{len(r['outcomes']) + 1}", kind, target_mids, ev)
    o["value"] = f(raw, norm, ev)
    if basis is not None:
        o["measurement_basis"] = f(basis, basis, ev)
    if yield_type is not None:
        o["yield_type"] = f(yield_type, yield_type, ev)
    if scope is not None:
        o["scope"] = f(scope, scope, ev)
    if product_role is not None:
        o["product_role"] = f(product_role, product_role, ev)
    r["outcomes"].append(o)
    return o


def set_solvent(r, label, molecule_id, ev, volume=None):
    charges = {"volume": (f"{volume} mL", q(volume, "mL"))} if volume is not None else {}
    add_participant(r, molecule_id, "solvent", ev, charges=charges)
    r["solvent_conditions"]["composition"] = f(label, [label], ev)
    if volume is not None:
        r["solvent_conditions"]["total_volume"] = f(f"{volume} mL", q(volume, "mL"), ev)


# Scheme 1: cited prior family and current generic scheme.
r = reaction("scheme1_prior_oxidative_domino", (E[1],), "literature_example", (),
             "Pd-catalyzed oxidative domino process via a vinyl-palladium intermediate",
             "cited prior-work reaction family", "Scheme 1a", (), "unresolved",
             ("Generic cited process; no concrete substrate, quantities, product, or yield is reported here.",))
add_participant(r, M1, "reactant", (E[1],))
add_participant(r, MINT, "product", (E[1],), product_role="proposed intermediate family")

GENERIC_ID = f"{DOC}:scheme1_current_generic"
r = reaction("scheme1_current_generic", (E[1], E[3]), "generic_scheme", (),
             "oxidative domino route to cyclopentenol enynes", "chemodivergent domino reaction",
             "Scheme 1b and generic header of Scheme 2", ("C-C", "C-O", "DB-C-C"), "partial",
             ("R/R' variables and proposed intermediates are not concrete experiments.",))
for molecule_id, role in ((M1, "reactant"), (MALK, "reactant"), (M2, "product"), (M3, "product"), (M4, "product"), (M5, "product")):
    add_participant(r, molecule_id, role, (E[1], E[3]), product_role="possible product family" if role == "product" else None)


# Table 1 optimization rows.
table1 = [
    (1, MPDOAC, "DCE", [10, None, 12, 69]),
    (2, MPDTFA, "DCE", [None, 22, 10, 60]),
    (3, MPDPH3, "DCE", [None, 89, None, None]),
    (4, MPDOAC, "THF", [5, None, 20, 71]),
    (5, MPDOAC, "Toluene", [None, None, 20, 74]),
    (6, MPDOAC, "CHCl3", [None, None, None, 90]),
    (7, MPDOAC, "Acetone", [None, 10, 25, 54]),
    (8, MPDOAC, "CH3CN", [None, None, 85, 3]),
    (9, MPDMCF, "DCE", [44, None, 5, 35]),
    (10, MPDMCF, "DCE", [None, 94, None, None]),
    (11, MPDMCF, "DCE", [84, None, None, 4]),
    (12, MPDMCF, "DCE", [83, None, None, 4]),
    (13, MPDCNC, "DCE", [80, None, None, 6]),
]
table1_records = {}
products = [P2["2a"], P3A, P4A, P5A]
product_labels = ["2a", "3a", "4a", "5a"]
for entry, catalyst, solvent, values in table1:
    ev = (E[2], E[3]) if entry == 12 else ((E[1], E[2]) if entry == 1 else (E[2],))
    cats = ["optimization"] + (["substrate_scope"] if entry == 12 else [])
    codes = []
    if values[0] is not None or values[2] is not None or values[3] is not None: codes.append("C-C")
    if values[1] is not None: codes.append("C-O")
    if values[3] is not None: codes.append("DB-C-C")
    r = reaction(f"table1_entry{entry}", ev, "specific_experiment", cats,
                 "optimization of the oxidative domino reaction", "chemodivergent oxidative domino reaction",
                 f"Table 1 entry {entry}" + ("; Scheme 2 product 2a duplicate" if entry == 12 else ""),
                 tuple(dict.fromkeys(codes)), "partial")
    r["general_scheme_ids"] = [GENERIC_ID]
    add_participant(r, M1A, "reactant", ev, charges={"amount": ("0.2 mmol", q(0.2, "mmol")), "equivalents": ("1.0 equiv", q(1.0, "equiv"))}, reference=M1A)
    add_participant(r, MPH, "reactant", ev, charges={"amount": ("0.25 mmol", q(0.25, "mmol")), "equivalents": ("1.25 equiv", q(1.25, "equiv"))}, reference=M1A)
    loading = 1 if entry in (12, 13) else 5
    add_participant(r, catalyst, "catalyst", ev, charges={"mol_percent": (f"{loading} mol%", q(loading, "mol%"))}, reference=M1A)
    add_participant(r, MBQ, "reagent", ev, charges={"equivalents": ("1.1 equiv", q(1.1, "equiv"))}, reference=M1A)
    if entry == 10:
        add_participant(r, MAGOTF, "additive", ev, charges={"mol_percent": ("5 mol%", q(5, "mol%"))}, reference=M1A)
    if entry in (11, 12, 13):
        add_participant(r, MET3N, "additive", ev, charges={"equivalents": ("0.1 equiv", q(0.1, "equiv"))}, reference=M1A)
    set_solvent(r, solvent, MSOL[solvent], ev, 1.0)
    for pmid in products:
        add_participant(r, pmid, "product", ev, product_role="reported product channel")
    add_participant(r, MANISOLE, "reagent", ev, stage="characterization")
    add_step(r, ev, "stir reaction mixture", ("room temperature", "r.t."),
             ("8 h" if entry in (12, 13) else "12 h", q(8 if entry in (12, 13) else 12, "h")))
    for label, pmid, value in zip(product_labels, products, values):
        if value is None:
            add_outcome(r, "other", [pmid], ev, f"Table 1 {label} cell: –", "table dash (–)", scope={"product": label})
        else:
            add_outcome(r, "yield", [pmid], ev, f"{value}%", q(value, "%"), basis="NMR spectroscopy using anisole as internal standard", yield_type="NMR yield", scope={"product": label})
    if entry == 1:
        add_outcome(r, "selectivity", [P2["2a"]], ev, ">20:1 d.r.", ">20:1", scope={"metric": "d.r.", "major": "cis-diastereomer", "order": "cis:other diastereomer(s)"})
    if entry == 12:
        add_outcome(r, "selectivity", [P2["2a"]], (E[3],), ">20:1 d.r.", ">20:1", scope={"metric": "d.r.", "major": "source-drawn cis product", "order": "cis:other diastereomer(s)"})
        r["notes"].append(f("The 83% 2a result in Scheme 2 matches Table 1 entry 12 and is represented once.", "cross-page duplicate merged", ev))
    table1_records[entry] = r


# Table 2: each catalyst block is a separate experiment.
table2_rows = [
    ("none", None, (44, 35), (5, 35), True),
    ("Pr3N", MPR3N, (82, 4), (46, 0), True),
    ("Pr2NH", MPR2NH, (47, 38), (28, 14), True),
    ("PrNH2", MPRNH2, (45, 40), (25, 15), True),
    ("TMEDA", MTMEDA, (0, 0), (0, 0), False),
]
for amine_label, amine, hetero, homo, partial_homo in table2_rows:
    for catalyst_label, catalyst, vals in (("Pd-AmP-MCF", MPDMCF, hetero), ("Pd(OAc)2", MPDOAC, homo)):
        slug = amine_label.lower().replace("2", "2").replace("3", "3") + "_" + ("mcf" if catalyst == MPDMCF else "pdoac2")
        r = reaction(f"table2_{slug}", (E[2],), "specific_experiment", ("optimization", "control"),
                     "amine effect under homogeneous and heterogeneous conditions", "amine-effect control",
                     f"Table 2 row {amine_label}, catalyst {catalyst_label}", ("C-C", "DB-C-C") if any(vals) else (),
                     "partial" if any(vals) else "not_applicable")
        r["general_scheme_ids"] = [GENERIC_ID]
        add_participant(r, M1A, "reactant", (E[2],), charges={"amount": ("0.2 mmol", q(0.2, "mmol")), "equivalents": ("1.0 equiv", q(1.0, "equiv"))}, reference=M1A)
        add_participant(r, MPH, "reactant", (E[2],), charges={"amount": ("0.25 mmol", q(0.25, "mmol")), "equivalents": ("1.25 equiv", q(1.25, "equiv"))}, reference=M1A)
        add_participant(r, catalyst, "catalyst", (E[2],), charges={"mol_percent": ("1 mol%", q(1, "mol%"))}, reference=M1A)
        if amine:
            add_participant(r, amine, "additive", (E[2],), charges={"equivalents": ("0.1 equiv", q(0.1, "equiv"))}, reference=M1A)
        add_participant(r, MBQ, "reagent", (E[2],), charges={"equivalents": ("1.1 equiv", q(1.1, "equiv"))}, reference=M1A)
        set_solvent(r, "DCE", MSOL["DCE"], (E[2],), 1.0)
        add_participant(r, P2["2a"], "product", (E[2],), product_role="product 2a")
        add_participant(r, P5A, "product", (E[2],), product_role="product 5a")
        add_step(r, (E[2],), "run amine-effect reaction", None, ("8 h", q(8, "h")))
        add_outcome(r, "yield", [P2["2a"]], (E[2],), f"{vals[0]}%", q(vals[0], "%"), scope={"product": "2a"})
        add_outcome(r, "yield", [P5A], (E[2],), f"{vals[1]}%", q(vals[1], "%"), scope={"product": "5a"})
        if amine_label == "TMEDA":
            add_outcome(r, "recovery", [M1A], (E[2],), ">90% recovery of 1a", q(90, "%", "gt", ">90%"), scope={"material": "1a"})
        elif catalyst == MPDOAC and partial_homo:
            add_outcome(r, "recovery", [M1A], (E[2],), "1a was partially recovered", "partially recovered", scope={"material": "1a"})


# Scheme 2 substrate scope, excluding 2a because it is merged with Table 1 entry 12.
scope_yields = {"2b": 85, "2c": 76, "2d": 80, "2e": 79, "2f": 70, "2g": 67, "2h": 76,
                "2i": 80, "2j": 83, "2k": 73, "2l": 75, "2m": 84, "2n": 81, "2o": 84,
                "2p": 85, "2q": 80, "2r": 84}
for label, yield_value in scope_yields.items():
    substrate = SCOPE_SUBSTRATES[label]
    alkyne = ALKYNES[label]
    r = reaction(f"scheme2_{label}", (E[3],), "specific_experiment", ("substrate_scope",),
                 "synthesis of cyclopentenol enyne " + label, "oxidative carbocyclization-alkynylation",
                 f"Scheme 2 product {label}", ("C-C",), "partial",
                 ("Yield and d.r. measurement methods are not reported on the supplied pages; complete product SMILES is not asserted.",))
    r["general_scheme_ids"] = [GENERIC_ID]
    add_participant(r, substrate, "reactant", (E[3],), charges={"equivalents": ("1.0 equiv (reference)", q(1.0, "equiv"))}, reference=substrate)
    add_participant(r, alkyne, "reactant", (E[3],), charges={"equivalents": ("1.25 equiv", q(1.25, "equiv"))}, reference=substrate)
    add_participant(r, MPDMCF, "catalyst", (E[3],), charges={"mol_percent": ("1 mol%", q(1, "mol%"))}, reference=substrate)
    add_participant(r, MBQ, "reagent", (E[3],), charges={"equivalents": ("1.1 equiv", q(1.1, "equiv"))}, reference=substrate)
    if label != "2r":
        add_participant(r, MET3N, "additive", (E[3],), charges={"equivalents": ("0.1 equiv", q(0.1, "equiv"))}, reference=substrate)
    set_solvent(r, "DCE", MSOL["DCE"], (E[3],))
    add_participant(r, P2[label], "product", (E[3],), product_role="desired cis product")
    add_step(r, (E[3],), "heterogeneous palladium-catalyzed domino reaction", ("room temperature", "r.t."), ("8 h", q(8, "h")))
    add_outcome(r, "yield", [P2[label]], (E[3],), f"{yield_value}%", q(yield_value, "%"), scope={"product": label})
    add_outcome(r, "selectivity", [P2[label]], (E[3],), ">20:1 d.r.", ">20:1", scope={"metric": "d.r.", "major": "source-drawn cis product", "order": "cis:other diastereomer(s)"})
    if label == "2r":
        r["notes"].append(f("Scheme 2 footnote a: Without Et3N.", "Et3N omitted", (E[3],)))


# Scheme 3 chemodivergent products. The compact R/R' identities remain source-labelled.
scheme3_yields = {
    "3": {"AgOTf_5mol_DCE": [94, 90, 95, 88]},
    "4": {"PdOAc2_5mol": [86, 77, 70, 80], "PdOAc2_1mol": [48, 30, 33, 40], "PdAmPMCF_1mol": [88, 80, 67, 84]},
    "5": {"PdOAc2_5mol": [90, 87, 84, 89], "PdOAc2_1mol": [44, 40, 36, 45], "PdAmPMCF_1mol": [59, 45, 44, 50]},
}
ratio_5_to_2 = ["59:40", "45:40", "44:42", "50:38"]
byproducts2 = [P2["2a"], P2["2p"], P2["2o"], P2["2n"]]
for series, variants in scheme3_yields.items():
    for variant, yields in variants.items():
        for i, letter in enumerate("abcd"):
            label = series + letter
            ev = (E[4], E[5]) if series == "5" and variant == "PdAmPMCF_1mol" else (E[4],)
            suffix = f"scheme3_{label}_{variant}"
            codes = ("C-O",) if series == "3" else (("C-C",) if series == "4" else ("C-C", "DB-C-C"))
            r = reaction(suffix, ev, "specific_experiment", ("substrate_scope",),
                         "chemodivergent synthesis of product " + label, "chemodivergent domino reaction",
                         f"Scheme 3 product {label}, condition {variant}", codes, "partial",
                         ("Specific compact R/R' reactant and product graphs are not fully encoded; yield method is not reported.",))
            r["general_scheme_ids"] = [GENERIC_ID]
            add_participant(r, M1, "reactant", ev, charges={"equivalents": ("1.0 equiv (reference)", q(1.0, "equiv"))}, reference=M1)
            add_participant(r, MALK, "reactant", ev, charges={"equivalents": ("1.25 equiv", q(1.25, "equiv"))}, reference=M1)
            add_participant(r, MBQ, "reagent", ev, charges={"equivalents": ("1.1 equiv", q(1.1, "equiv"))}, reference=M1)
            if variant == "AgOTf_5mol_DCE":
                add_participant(r, MAGOTF, "catalyst", ev, charges={"mol_percent": ("5 mol%", q(5, "mol%"))}, reference=M1)
                set_solvent(r, "DCE", MSOL["DCE"], ev)
            elif variant.startswith("PdOAc2"):
                loading = 5 if "5mol" in variant else 1
                add_participant(r, MPDOAC, "catalyst", ev, charges={"mol_percent": (f"{loading} mol%", q(loading, "mol%"))}, reference=M1)
                set_solvent(r, "CH3CN" if series == "4" else "CHCl3", MSOL["CH3CN"] if series == "4" else MSOL["CHCl3"], ev)
            else:
                add_participant(r, MPDMCF, "catalyst", ev, charges={"mol_percent": ("1 mol%", q(1, "mol%"))}, reference=M1)
                set_solvent(r, "CH3CN" if series == "4" else "CHCl3", MSOL["CH3CN"] if series == "4" else MSOL["CHCl3"], ev)
            product_id = SCHEME3_PRODUCTS[label]
            add_participant(r, product_id, "product", ev, product_role="reported product")
            add_step(r, ev, "chemodivergent reaction", ("room temperature", "r.t."), ("8 h", q(8, "h")))
            add_outcome(r, "yield", [product_id], ev, f"{yields[i]}%", q(yields[i], "%"), scope={"product": label})
            if "PdOAc2_1mol" == variant:
                add_outcome(r, "recovery", [M1], ev, "starting materials were partially recovered", "partially recovered", scope={"material": "starting materials"})
            if series == "5" and variant == "PdAmPMCF_1mol":
                byproduct = byproducts2[i]
                add_participant(r, byproduct, "product", ev, product_role="product 2 byproduct")
                add_outcome(r, "conversion", [M1], (E[4],), ">90% conversion", q(90, "%", "gt", ">90%"), scope={"starting material": "enallene 1"})
                add_outcome(r, "selectivity", [product_id, byproduct], ev, ratio_5_to_2[i], ratio_5_to_2[i], scope={"metric": "product ratio", "order": f"{label}:corresponding product 2"})


# Scheme 4a asymmetric syntheses.
MR1C = molecule("R_1c", "(R)-1c", "enantiomerically pure enallene (R)-1c", (E[4],), smiles=None,
                 chirality="specified", geometry="unspecified", representation="not_encoded",
                 atom_centers=[{"source_label": "R", "site": "source-drawn alcohol-bearing center"}],
                 assignment="(R)-1c as labelled; exact atom mapping is not encoded")
ASYM_PRODUCTS = {
    "2n": molecule("asym_RS_2n", "(R,S)-2n", "optically pure (R,S)-2n", (E[4],), smiles=None,
                   chirality="specified", geometry="not_applicable", representation="not_encoded",
                   atom_centers=[{"source_label": "R,S", "status": "specified by source label"}],
                   assignment="(R,S)-2n; source-drawn relative cis product"),
    "4e": molecule("asym_R_4e", "(R)-4e", "optically pure (R)-4e", (E[4],), smiles=None,
                   chirality="specified", geometry="not_applicable", representation="not_encoded",
                   atom_centers=[{"source_label": "R", "status": "specified by source label"}], assignment="(R)-4e"),
    "5e": molecule("asym_R_5e", "(R)-5e", "optically pure (R)-5e", (E[4],), smiles=None,
                   chirality="specified", geometry="not_applicable", representation="not_encoded",
                   atom_centers=[{"source_label": "R", "status": "specified by source label"}], assignment="(R)-5e"),
}
asym_specs = [
    ("2n", 79, MPDMCF, "DCE", True, ("C-C",), ">20:1 d.r."),
    ("4e", 82, MPDOAC, "CH3CN", False, ("C-C",), None),
    ("5e", 87, MPDOAC, "CHCl3", False, ("C-C", "DB-C-C"), None),
]
for label, yld, catalyst, solvent, use_amine, codes, dr in asym_specs:
    r = reaction(f"scheme4a_asymmetric_{label}", (E[4],), "specific_experiment", ("substrate_scope",),
                 "asymmetric synthesis of " + label, "enantiospecific chemodivergent domino reaction",
                 f"Scheme 4a product {label}", codes, "partial",
                 ("Catalyst loading, reagent quantities, time, temperature, and yield/ee measurement methods are not reported in Scheme 4a.",))
    r["general_scheme_ids"] = [GENERIC_ID]
    add_participant(r, MR1C, "reactant", (E[4],))
    add_participant(r, MPH, "reactant", (E[4],))
    add_participant(r, MBQ, "reagent", (E[4],))
    add_participant(r, catalyst, "catalyst", (E[4],))
    if use_amine:
        add_participant(r, MET3N, "additive", (E[4],))
    set_solvent(r, solvent, MSOL[solvent], (E[4],))
    add_participant(r, ASYM_PRODUCTS[label], "product", (E[4],), product_role="optically pure product")
    add_step(r, (E[4],), "asymmetric chemodivergent transformation")
    add_outcome(r, "yield", [ASYM_PRODUCTS[label]], (E[4],), f"{yld}%", q(yld, "%"), scope={"product": label})
    if dr:
        add_outcome(r, "selectivity", [ASYM_PRODUCTS[label]], (E[4],), dr, ">20:1", scope={"metric": "d.r.", "order": "source-drawn major:other diastereomer(s)"})
    add_outcome(r, "selectivity", [ASYM_PRODUCTS[label]], (E[4],), ">99% ee", q(99, "%", "gt", ">99%"), scope={"metric": "ee", "configuration": label.split(")")[0] + ")"})


# Scheme 4b borylation and footnote-38 unsuccessful aromatic-R case.
M1I = molecule("1i", "1i", "aromatic R-substituted enallenol 1i (R=Ph)", (E[4], E[5]), smiles=None,
                chirality="unspecified", geometry="unspecified", representation="not_encoded",
                identity_note="Source-drawn specific enallenol with R=Ph; compact graph not encoded.")
P6 = molecule("product_6", "6", "cyclopentenol boron compound 6", (E[4],), smiles=None,
              chirality="partially_specified", geometry="not_applicable", representation="not_encoded",
              assignment="source-drawn relative configuration; >20:1 d.r. reported at experiment level")
r = reaction("scheme4b_borylation_6", (E[4],), "specific_experiment", ("substrate_scope", "derivatization"),
             "oxidative carbocyclization-borylation", "oxidative carbocyclization-borylation",
             "Scheme 4b", ("C-C", "C-B"), "partial",
             ("Yield and d.r. measurement methods are not reported; complete structures and atom mapping remain unresolved.",))
r["general_scheme_ids"] = [GENERIC_ID]
add_participant(r, M1I, "reactant", (E[4],), charges={"equivalents": ("1.0 equiv (reference)", q(1, "equiv"))}, reference=M1I)
add_participant(r, MB2PIN2, "reactant", (E[4],), charges={"equivalents": ("1.25 equiv", q(1.25, "equiv"))}, reference=M1I)
add_participant(r, MPDMCF, "catalyst", (E[4],), charges={"mol_percent": ("1 mol%", q(1, "mol%"))}, reference=M1I)
add_participant(r, MBQ, "reagent", (E[4],), charges={"equivalents": ("1.1 equiv", q(1.1, "equiv"))}, reference=M1I)
add_participant(r, MET3N, "additive", (E[4],), charges={"equivalents": ("0.1 equiv", q(0.1, "equiv"))}, reference=M1I)
set_solvent(r, "DCE", MSOL["DCE"], (E[4],))
add_participant(r, P6, "product", (E[4],), product_role="desired boron compound")
add_step(r, (E[4],), "oxidative carbocyclization-borylation", ("room temperature", "r.t."), ("8 h", q(8, "h")))
add_outcome(r, "yield", [P6], (E[4],), "89%", q(89, "%"), scope={"product": "6"})
add_outcome(r, "selectivity", [P6], (E[4],), ">20:1 d.r.", ">20:1", scope={"metric": "d.r.", "order": "source-drawn major:other diastereomer(s)"})

MMESSY = generic_molecule("messy_mixture_1i", "messy mixture", "unidentified mixture from aromatic R-substituted enallenol 1i and a terminal alkyne", (E[5],), ())
r = reaction("footnote38_aromatic_R_1i_failure", (E[5],), "specific_experiment", ("substrate_scope", "control"),
             "unsuccessful reaction of aromatic R-substituted enallenol 1i", "failed scope reaction",
             "Footnote 38", (), "unresolved",
             ("The terminal alkyne identity and exact conditions are not stated in footnote 38; unavailable SI was not inferred. Product identities and yields are unknown.",))
r["general_scheme_ids"] = [GENERIC_ID]
add_participant(r, M1I, "reactant", (E[5],))
add_participant(r, MALK, "reactant", (E[5],))
add_participant(r, MMESSY, "product", (E[5],), product_role="unidentified mixture")
add_step(r, (E[5],), "reaction under the reported reaction conditions; exact conditions not restated")
add_outcome(r, "other", [MMESSY], (E[5],), "messy mixture as determined by TLC and 1H NMR spectroscopy", "messy mixture", basis="TLC and 1H NMR spectroscopy", scope={"observation": "unidentified mixture"})


# Scheme 5 is an author-proposed mechanism, not a set of executed reactions.
r = reaction("scheme5_proposed_mechanism", (E[4],), "proposed_mechanism", (),
             "proposed mechanism for the Pd-catalyzed domino reactions", "author-proposed catalytic cycle",
             "Scheme 5", (), "unresolved",
             ("Int-1, Int-2, and Int-3 are proposed intermediates; no elementary step is represented as a concrete experiment.",))
for molecule_id, role in ((M1, "reactant"), (MALK, "reactant"), (MINT, "product"), (M2, "product"), (M4, "product"), (M5, "product")):
    add_participant(r, molecule_id, role, (E[4],), product_role="proposed species/product family" if role == "product" else None)


bundle["molecules"] = list(molecules.values())
bundle["reactions"] = reactions
bundle["source_context"] = [
    {"context_id": f"{DOC}:context:metadata", "document_id": DOC, "kind": "bibliographic_metadata",
     "evidence_ids": [E[1]], "title": "Amino-Supported Palladium Catalyst for Chemo- and Stereoselective Domino Reactions",
     "doi": "10.1002/anie.202011708"},
    {"context_id": f"{DOC}:context:xps", "document_id": DOC, "kind": "catalyst_characterization",
     "evidence_ids": [E[2], E[3], E[5]],
     "observation": "Figure 1 and footnote 37 report Pd 3d XPS comparisons before/after homogeneous and heterogeneous reactions; PdII/Pd0 for Pd-AmP-MCF changes from about 90/10 to about 85/15 after run 1 and is maintained through run 7.",
     "decision": "Retained as catalyst/spectroscopic context, not a separate synthesis record."},
    {"context_id": f"{DOC}:context:heterogeneity", "document_id": DOC, "kind": "catalyst_performance_context",
     "evidence_ids": [E[3]],
     "observation": "ICP-OES (<0.1 ppm Pd), hot filtration, recycling through seven cycles, and kinetic plots support the heterogeneous-catalyst discussion.",
     "decision": "The supplied pages do not give sufficient discrete protocols/row values to create reliable additional synthesis records; observations remain source context."},
    {"context_id": f"{DOC}:context:no_mass_spectrum", "document_id": DOC, "kind": "negative_inventory",
     "evidence_ids": [E[1], E[2], E[3], E[4], E[5]],
     "observation": "No standalone mass-spectrum peak table, molecular-weight table, or fragment-relative-abundance block occurs on pages 1-5."},
    {"context_id": f"{DOC}:context:si", "document_id": DOC, "kind": "supplement_status",
     "evidence_ids": [E[1], E[5]],
     "observation": "The paper points to Supporting Information, but no SI file was supplied in this workspace; association remains unestablished."},
]

OUT.write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")
(RUN / "raw_model/stereo_checks.json").write_text(json.dumps(stereo_diagnostics, ensure_ascii=False, indent=2) + "\n")
reaction_ids = [r["reaction_id"] for r in reactions]
def ids(prefixes):
    return [rid for rid in reaction_ids if any(rid.split(":", 1)[1].startswith(p) for p in prefixes)]

coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {"page": 1, "state": "attempted", "regions": [
            {"region": "title, abstract, and introductory prose", "finding": "Paper identity and discussion of the current chemodivergent oxidative domino project; no standalone experimental procedure.", "evidence_ids": [E[1]]},
            {"region": "right-column initial-attempt paragraph", "finding": "Table 1 entry 1 is repeated in prose (2a 10%, cis >20:1 d.r.; 5a 69%) and merged with that table row.", "evidence_ids": [E[1]]},
            {"region": "Scheme 1a-c", "finding": "One cited prior-work generic family, one current generic chemodivergent scheme, proposed intermediates, and bioactive context; no extra concrete experiment invented.", "evidence_ids": [E[1]]},
        ], "record_ids": ids(["scheme1_", "table1_entry1"]),
         "unresolved": ["Scheme 1 R/R' structures and intermediates are generic; exact atom mappings are not asserted."]},
        {"page": 2, "state": "attempted", "regions": [
            {"region": "Table 1 including structures, 13 rows, and footnotes a-e", "finding": "All 13 optimization rows recorded independently; dashes preserved qualitatively, additives/loadings and NMR/anisole method structured.", "evidence_ids": [E[2]]},
            {"region": "Table 1 discussion", "finding": "Catalyst/solvent selectivity discussion and recovery/Pd-black observations mapped to existing rows or source context.", "evidence_ids": [E[2]]},
            {"region": "Table 2 including five amine rows and two catalyst blocks", "finding": "Each row contains two distinct catalyst experiments; represented as 10 records with 2a/5a outcomes and recovery statements.", "evidence_ids": [E[2]]},
        ], "record_ids": ids(["table1_", "table2_"]),
         "unresolved": ["Table 2 does not report yield measurement basis or reaction temperature; homogeneous partial recoveries have no numerical percentage except TMEDA >90%."]},
        {"page": 3, "state": "attempted", "regions": [
            {"region": "Figure 1 XPS spectra and catalyst illustration", "finding": "Catalyst characterization retained in source_context, not converted to a synthesis record.", "evidence_ids": [E[3]]},
            {"region": "Scheme 2 generic header and products 2a-r", "finding": "All 18 product outcomes inventoried; 2a merged with identical Table 1 entry 12, and 2b-r recorded separately with yields and >20:1 d.r.; 2r omits Et3N per footnote.", "evidence_ids": [E[3]]},
            {"region": "Figure 2 recycling/kinetic plot and adjacent ICP/hot-filtration prose", "finding": "Performance and heterogeneity observations retained in source_context because discrete protocols and exact plotted values are not supplied on these pages.", "evidence_ids": [E[3]]},
        ], "record_ids": [GENERIC_ID, f"{DOC}:table1_entry12"] + ids(["scheme2_"]),
         "unresolved": ["Scheme 2 yield and d.r. measurement methods are unreported; several compact structures remain unencoded rather than guessed."]},
        {"page": 4, "state": "attempted", "regions": [
            {"region": "Scheme 3 chemodivergent syntheses", "finding": "Four 3-series AgOTf experiments and three catalyst conditions for each 4a-d and 5a-d were separated, totaling 28 records; partial recovery/conversion footnotes retained.", "evidence_ids": [E[4]]},
            {"region": "Scheme 4a asymmetric syntheses", "finding": "Three independent products 2n, 4e, and 5e recorded with source configurations, yields, ee, and 2n d.r.", "evidence_ids": [E[4]]},
            {"region": "Scheme 4b borylation", "finding": "One oxidative carbocyclization-borylation experiment recorded with all displayed loadings, 89% yield, and >20:1 d.r.", "evidence_ids": [E[4]]},
            {"region": "Scheme 5 and mechanistic prose", "finding": "Author-proposed catalytic cycle represented once as proposed_mechanism; Int-1 to Int-3 are not executed experiments.", "evidence_ids": [E[4]]},
        ], "record_ids": ids(["scheme3_", "scheme4", "scheme5_"]),
         "unresolved": ["Scheme 3/4 yield/selectivity measurement methods and many exact compact structures are not reported/encoded; Scheme 4a omits loadings, time, and temperature."]},
        {"page": 5, "state": "attempted", "regions": [
            {"region": "acknowledgements, conflict statement, keywords, and references", "finding": "No current-paper synthesis arrow or experimental procedure; bibliography was not imported as experiments.", "evidence_ids": [E[5]]},
            {"region": "footnotes 29-42", "finding": "Footnote 38 supplies one failed aromatic-R substrate experiment; footnote 41 supplies product 5:2 ratios for four Scheme 3 records; XPS/SI notices remain source context.", "evidence_ids": [E[5]]},
            {"region": "full-page negative analytical inventory", "finding": "No standalone mass-spectrum peak table, molecular-weight table, or fragment-relative-abundance block.", "evidence_ids": [E[5]]},
        ], "record_ids": [f"{DOC}:footnote38_aromatic_R_1i_failure"] + ids(["scheme3_5a_PdAmPMCF", "scheme3_5b_PdAmPMCF", "scheme3_5c_PdAmPMCF", "scheme3_5d_PdAmPMCF"]),
         "unresolved": ["Footnote 38 does not identify the terminal alkyne or exact conditions; cited SI/CCDC-style external details were not supplied and were not inferred."]},
    ],
}
(RUN / "coverage/pages_0001_0005.json").write_text(json.dumps(coverage, ensure_ascii=False, indent=2) + "\n")

page_extraction = {
    "document_id": DOC,
    "source_pdf": "assets/doc_151fad1677f63120/source.pdf",
    "pages_attempted": [1, 2, 3, 4, 5],
    "image_evidence_ids": [E[p] for p in range(1, 6)],
    "visual_crops": [
        "raw_model/crops/page1_scheme1.jpg", "raw_model/crops/page2_table1.jpg", "raw_model/crops/page2_table2.jpg",
        "raw_model/crops/page3_scheme2.jpg", "raw_model/crops/page4_scheme3.jpg", "raw_model/crops/page4_scheme4.jpg",
        "raw_model/crops/page4_scheme5.jpg",
    ],
    "inventory": {
        "page1": {"table_rows": 0, "concrete_arrows": 0, "generic_or_prior_schemes": 2, "decisions": ["Table 1 entry 1 prose duplicate merged"]},
        "page2": {"table1_rows": 13, "table2_physical_rows": 5, "table2_distinct_catalyst_experiments": 10},
        "page3": {"scheme2_drawn_products": 18, "new_scheme2_records": 17, "duplicates": ["2a = Table 1 entry 12"]},
        "page4": {"scheme3_specific_records": 28, "scheme4a_specific_records": 3, "scheme4b_specific_records": 1, "proposed_mechanisms": 1},
        "page5": {"new_specific_records": 1, "enriched_existing_records": 4, "synthesis_arrows": 0},
    },
    "deduplication": [
        {"source_locations": ["Table 1 entry 12", "Scheme 2 product 2a"], "record_id": f"{DOC}:table1_entry12", "reason": "same displayed reagents, loadings, solvent, temperature, time, and 83% 2a outcome"},
        {"source_locations": ["page-1 initial-attempt prose", "Table 1 entry 1"], "record_id": f"{DOC}:table1_entry1", "reason": "same Pd(OAc)2/DCE initial experiment and 2a/5a yields"},
    ],
    "semantic_accuracy_verified": False,
    "chemical_accuracy_verified": False,
}
(RUN / "raw_model/page_extraction.json").write_text(json.dumps(page_extraction, ensure_ascii=False, indent=2) + "\n")
print(json.dumps({"molecules": len(bundle["molecules"]), "reactions": len(bundle["reactions"]),
                  "source_context": len(bundle["source_context"]), "out": str(OUT)}, ensure_ascii=False, indent=2))
