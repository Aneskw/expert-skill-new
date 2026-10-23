#!/usr/bin/env python3
"""Build the evidence-linked v3.1 submission for Dong and Breit, Org. Lett. 2021.

This is a paper-specific transcription artifact. It does not modify extraction tools,
schemas, taxonomies, or quality rules. All generated records remain needs_review.
"""
from __future__ import annotations

import copy
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

RUN = ROOT / "runs/current"
DOC = "doc_833ef6f3ab8863df"
E = {
    p: [f"{DOC}:page{p}:image", f"{DOC}:page{p}:text"] for p in range(1, 6)
}
SENTINEL = object()
diagnostics = {}


def f(raw=None, normalized=SENTINEL, evidence=(), method="original_extraction",
      status=None, derivation=None, alternatives=()):
    if normalized is SENTINEL:
        normalized = raw
    return st.fact(raw, normalized, evidence, method=method, status=status,
                   derivation=derivation, alternatives=alternatives)


def inferred(raw, normalized, evidence, derivation):
    return f(raw, normalized, evidence, method="inference", derivation=derivation)


def quantity(number, unit, qualifier=None, raw_text=None):
    value = {"number": number, "unit": unit}
    if qualifier:
        value["qualifier"] = qualifier
    if raw_text:
        value["raw_text"] = raw_text
    return value


def mid(label):
    safe = label.replace("′", "prime").replace("'", "prime").replace("-", "_").replace(" ", "_")
    return f"{DOC}:mol_{safe}"


def set_stereo_status(m, chirality, geometry, representation, ev, reason):
    m["stereochemistry"]["chirality_status"] = inferred(reason, chirality, ev, reason)
    m["stereochemistry"]["geometry_status"] = inferred(reason, geometry, ev, reason)
    m["stereochemistry"]["representation_status"] = inferred(reason, representation, ev, reason)


def molecule(label, ev, name=None, *, chirality="unresolved", geometry="unresolved",
             representation="not_encoded", note=None, smiles=None, axial=None,
             mixture=None, abbreviation=None):
    m = st.new_molecule(mid(label), DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f(label, label, ev)
    m["names"]["original"] = f(name or label, name or label, ev)
    m["names"]["english"] = f(name or label, name or label, ev)
    if abbreviation:
        m["names"]["abbreviations"] = f(abbreviation, [abbreviation], ev)
    reason = note or "Source label/drawing inspected; a complete stereochemical graph was not safely encoded."
    m["identity_note"] = inferred({"label": label, "source": "paper drawing or name"}, reason, ev, reason)
    if axial is not None:
        m["stereochemistry"]["axial"] = inferred("Source allene drawing", axial, ev, reason)
    if mixture is not None:
        m["stereochemistry"]["mixture"] = inferred("Source reports a mixture/result distribution", mixture, ev, reason)
    if geometry == "mixture":
        m["stereochemistry"]["double_bonds"] = inferred(
            "E/Z ratio is reported for this product in experiment outcomes",
            [{"bond": "allylic product C=C", "assignment": "E/Z mixture; ratio is experiment-specific"}],
            ev, "The product label denotes the E/Z material; exact ratios are stored per experiment."
        )
    if chirality == "unspecified":
        m["stereochemistry"]["assignment"] = inferred(
            "No wedge/dash or absolute assignment is printed for the depicted potentially stereogenic product",
            {"assignment": "unspecified"}, ev,
            "Potential stereogenicity is retained without inventing racemic or absolute configuration."
        )
    if smiles:
        diag = st.stereo_inspect(smiles)
        diagnostics[label] = diag
        m["representation_type"] = "smiles"
        m["smiles"] = f(
            {"source": name or label, "candidate_smiles": smiles}, smiles, ev,
            method="structure_derivation",
            derivation="Name/source drawing converted to a machine candidate; RDKit parsing is not source verification."
        )
        m["stereochemistry"]["isomeric_smiles"] = f(
            {"candidate_smiles": smiles}, diag["canonical_isomeric_smiles"], ev,
            method="structure_derivation", derivation="Generated with the standard stereo diagnostic."
        )
        m["stereochemistry"]["identity_key"] = f(
            {"canonical_isomeric_smiles": diag["canonical_isomeric_smiles"]}, diag["stereo_graph_key"], ev,
            method="structure_derivation", derivation="Standard stereo graph key; not a sample identity."
        )
    set_stereo_status(m, chirality, geometry, representation, ev, reason)
    return m


def generic_molecule(label, ev, fixed, variables, reason):
    m = molecule(label, ev, name=label, chirality="unresolved", geometry="unresolved",
                 representation="not_encoded", note=reason)
    m["representation_type"] = "generic_template"
    m["generic_template"] = {
        "fixed_structure": f(fixed, fixed, ev),
        "variables": [
            {
                "label": var,
                "attachment_points": inferred(f"Attachment point for {var} is shown in the scheme", [var], ev,
                                               "Attachment is retained descriptively because the drawing is generic."),
                "definition": f(desc, desc, ev),
                "allowed_values": f(),
                "constraints": f(),
            }
            for var, desc in variables
        ],
        "instance_mapping": [],
        "reason": f(reason, reason, ev),
    }
    return m


def reaction(suffix, ev, record_type, categories, location, reaction_name, transformation,
             general=()):
    r = st.new_reaction(f"{DOC}:{suffix}", DOC, ev)
    r["record_type"] = record_type
    r["experiment_categories"] = f(categories, categories, ev)
    r["location"] = f(location, location, ev)
    r["general_scheme_ids"] = list(general)
    r["reaction_name"] = f(reaction_name, reaction_name, ev)
    r["transformation_type"] = inferred(
        reaction_name, transformation, ev,
        "Transformation description derived from the depicted substrate/product connectivity and article wording."
    )
    r["review"]["unresolved"] = [
        "No explicitly associated SI was supplied; detailed workup, purification, and characterization remain unavailable.",
        "Machine transcription and category-level bond changes require expert source review."
    ]
    return r


def participant(pid, molecule_label, role, ev, *, amount=None, volume=None,
                equivalents=None, mol_percent=None, ref=None, product_role=None):
    p = st.new_participant(pid, mid(molecule_label), role, ev)
    if product_role:
        p["product_role"] = f(product_role, product_role, ev)
    for key, item in (("amount", amount), ("volume", volume), ("equivalents", equivalents),
                      ("mol_percent", mol_percent)):
        if item is not None:
            raw, normalized = item
            p["charges"][key] = f(raw, normalized, ev)
    if ref is not None:
        p["equivalence_reference"] = inferred(
            f"Loading/equivalents referenced to {ref}", mid(ref), ev,
            "Reference substrate assigned from the scheme/table loading convention."
        )
    return p


def add_conditions(r, ev, solvent, *, temperature="rt", time_h=None, atmosphere=None,
                   volume_ml=None, concentration_m=None, inherited=False):
    psol = participant(f"p_{solvent}", solvent, "solvent", ev,
                       volume=(f"{volume_ml} mL", quantity(volume_ml, "mL")) if volume_ml is not None else None)
    r["participants"].append(psol)
    method = "inherited_general_condition" if inherited else "original_extraction"
    r["solvent_conditions"]["composition"] = f(solvent, [mid(solvent)], ev, method=method)
    if volume_ml is not None:
        r["solvent_conditions"]["total_volume"] = f(f"{volume_ml} mL", quantity(volume_ml, "mL"), ev, method=method)
    if concentration_m is not None:
        r["solvent_conditions"]["concentration"] = f(
            f"{solvent} ({concentration_m} M)", quantity(concentration_m, "M"), ev, method=method
        )
    step = st.new_step("reaction", ev)
    step["participant_ids"] = [p["participant_id"] for p in r["participants"]]
    step["operation"] = f("Reaction under the conditions printed in the scheme/table",
                          "combine/hold under stated conditions", ev, method=method)
    if temperature == "rt":
        step["temperature"] = f("room temperature (rt)", "rt", ev, method=method)
    elif temperature is not None:
        step["temperature"] = f(str(temperature), temperature, ev, method=method)
    if time_h is not None:
        step["time"] = f(f"{time_h} h", quantity(time_h, "h"), ev, method=method)
    if atmosphere:
        step["atmosphere"] = f(atmosphere, atmosphere, ev, method=method)
    r["steps"] = [step]


def classification(r, ev, observed=True):
    c = r["classification"]
    if not observed:
        c["labels"] = inferred("0% yield/no observed product", [], ev,
                               "No formation label is assigned to an unsuccessful experiment.")
        c["coverage"] = inferred("No observed transformation", "not_applicable", ev,
                                 "Negative/zero-yield experiment; expected product is not treated as observed chemistry.")
        c["limitations"] = ["No bond formation inferred from an intended but unobserved transformation."]
        return
    changes = [
        {"change_id": "c1", "atoms": [{"element": "C", "source_position": "benzimidazole C2"},
                                           {"element": "C", "source_position": "allene-derived allyl carbon"}],
         "before_order": 0, "after_order": 1, "change_type": "formed"},
        {"change_id": "c2", "atoms": [{"element": "C", "source_position": "allene-derived carbon"},
                                           {"element": "H", "source_position": "net transferred hydrogen; exact source unresolved"}],
         "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain", "hydrogen_count": 1},
    ]
    c["bond_changes"] = [
        inferred("Category-level change inferred from the depicted net transformation", ch, ev,
                 "No atom-complete mapping is asserted; connectivity and implicit-H change are partial.")
        for ch in changes
    ]
    labels = [
        {"category_number": 2, "code": "C-C", "supporting_change_ids": ["c1"]},
        {"category_number": 1, "code": "C-H", "supporting_change_ids": ["c2"]},
    ]
    c["labels"] = inferred("Depicted C2 allylation/hydroallylation", labels, ev,
                           "Labels follow the 41-category taxonomy from the two retained partial bond changes.")
    c["coverage"] = inferred("Category-level only", "partial", ev,
                             "Dense structures were not atom-mapped and the proton source is not assigned.")
    c["limitations"] = [
        "Category-level net changes only; no atom-complete mapped reaction is claimed.",
        "Implicit-H gain is retained, but its exact atom provenance and all bond-order changes remain unresolved.",
        "Unchanged product bonds and reduced allene bond order are not mislabeled as newly formed sigma bonds."
    ]


def add_yield(r, oid, product, raw, value, ev, *, isolated=True, role="main_product", scope=None):
    kind = "trace" if value == "trace" else "yield"
    o = st.new_outcome(oid, kind, [mid(product)], ev)
    o["product_role"] = f(role, role, ev)
    o["value"] = f(raw, value, ev)
    if isolated:
        o["yield_type"] = f("isolated yield", "isolated", ev)
        o["measurement_basis"] = f("isolated yield", "isolated", ev)
    if scope:
        o["scope"] = f(scope, scope, ev)
    r["outcomes"].append(o)


def add_ambiguous_yield(r, oid, product, raw, ev):
    o = st.new_outcome(oid, "yield", [mid(product)], ev)
    o["product_role"] = f("product", "main_product", ev)
    o["value"] = f(raw, None, ev, status="unreadable")
    o["scope"] = inferred(raw, {"note": "printed beneath product without an explicit unit or yield label"}, ev,
                          "The numeral is retained but not normalized to percent.")
    r["outcomes"].append(o)


def add_selectivity(r, oid, product, ratio, ev, method_text):
    rel = "gt" if ratio.startswith(">") else "exact"
    clean = ratio.lstrip(">")
    a, b = [int(x) for x in clean.split(":")]
    o = st.new_outcome(oid, "selectivity", [mid(product)], ev)
    o["product_role"] = f("E/Z product mixture", "main_product_mixture", ev)
    o["value"] = f(f"{ratio} (E/Z)", {"metric": "E/Z", "order": ["E", "Z"],
                                           "ratio": [a, b], "qualifier": rel}, ev)
    o["measurement_basis"] = f(method_text, method_text, ev)
    o["scope"] = f("E/Z", {"metric": "E/Z", "isomer_order": ["E", "Z"],
                                  "absolute_configuration": "not assigned"}, ev)
    r["outcomes"].append(o)


def add_isotope_outcome(r, oid, product, percent, ev):
    o = st.new_outcome(oid, "other", [mid(product)], ev)
    o["product_role"] = f("deuterium-labeled product", "main_product", ev)
    o["value"] = f(f"{percent}% D", quantity(percent, "%"), ev)
    o["scope"] = f("deuterium incorporation at the alkene position shown in Scheme 4",
                    {"metric": "D incorporation", "position": "source-depicted alkene position"}, ev)
    r["outcomes"].append(o)


bundle = json.loads((RUN / "submitted.json").read_text())
molecules = []

# Source-named reagents and solvents. Only low-risk candidates are encoded; parsing is not source verification.
concrete = [
    ("THF", "tetrahydrofuran", "C1CCOC1", E[1]),
    ("DCE", "1,2-dichloroethane", "ClCCCl", E[2]),
    ("CH3CN", "acetonitrile", "CC#N", E[2]),
    ("toluene", "toluene", "Cc1ccccc1", E[1]),
    ("1,4-dioxane", "1,4-dioxane", "O1CCOCC1", E[2]),
    ("MTBE", "methyl tert-butyl ether", "COC(C)(C)C", E[1]),
    ("(MeO)2MeSiH", "dimethoxymethylsilane", "CO[SiH](C)OC", E[1]),
    ("(EtO)2MeSiH", "diethoxymethylsilane", "CCO[SiH](C)OCC", E[2]),
    ("3a", "1-methylbenzimidazole", "Cn1c2ccccc2nc1", E[2] + E[3]),
]
for label, name, smi, ev in concrete:
    molecules.append(molecule(label, ev, name=name, smiles=smi, chirality="achiral",
                              geometry="not_applicable", representation="encoded",
                              note="Candidate graph has no tetrahedral stereocenter or E/Z element; this is a derived graph diagnostic."))
for label, name, smi, ev in [
    ("4a", "cyclohexylallene", "C=C=CC1CCCCC1", E[2] + E[3]),
    ("4b", "cyclopentylallene", "C=C=CC1CCCC1", E[3]),
    ("4c", "phenylallene", "C=C=Cc1ccccc1", E[3]),
]:
    molecules.append(molecule(label, ev, name=name, smiles=smi, chirality="achiral",
                              geometry="not_applicable", representation="partial",
                              axial=[{"axis": "C=C=C", "assignment": "not axially stereogenic because one terminus is CH2"}],
                              note="Terminal monosubstituted allene candidate; source has no axial assignment and toolkit axial support is limited."))

for label, name, abbr, ev, chirality in [
    ("Cu(OAc)2", "copper(II) acetate", None, E[1] + E[2] + E[3], "unresolved"),
    ("CuOAc", "copper(I) acetate", None, E[1], "unresolved"),
    ("Xantphos", "Xantphos", "Xantphos", E[1] + E[2] + E[3], "unresolved"),
    ("DCyPE", "1,2-bis(dicyclohexylphosphino)ethane", "DCyPE", E[1] + E[2], "unresolved"),
    ("(S,S)-Ph-BPE", "(S,S)-Ph-BPE", "(S,S)-Ph-BPE", E[1], "specified"),
]:
    m = molecule(label, ev, name=name, abbreviation=abbr, chirality=chirality,
                 geometry="not_applicable", representation="not_encoded",
                 note="Named/drawn catalyst or ligand retained without an unchecked complete SMILES.")
    if chirality == "specified":
        m["stereochemistry"]["assignment"] = f("(S,S)-Ph-BPE", {"source_assignment": "(S,S)"}, ev)
    molecules.append(m)

molecules.append(generic_molecule(
    "PMHS", E[1], "poly(methylhydrosiloxane) repeating siloxane/hydrosilane material",
    [("n", "polymer repeat count and distribution are not reported")],
    "Named polymeric hydrosilane; composition/molecular-weight distribution is not specified in the main text."
))

# Generic reaction families from Scheme 1.
generic_specs = [
    ("lit1a_substrate", "N-OPiv benzimidazole plus 1,3-diene family", [("R1", "benzimidazole substituent"), ("Ar", "aryl group")]),
    ("lit1a_product", "branched allylated benzimidazole family", [("R1", "benzimidazole substituent"), ("Ar", "aryl group")]),
    ("lit1b_substrate", "N-substituted benzimidazole plus alkene family", [("R1", "benzimidazole substituent"), ("R2", "N substituent"), ("Ar", "aryl group")]),
    ("lit1b_product", "branched alkylated benzimidazole family", [("R1", "benzimidazole substituent"), ("R2", "N substituent"), ("Ar", "aryl group")]),
    ("lit1b_tethered_substrate", "allene/alkene-tethered N-substituted benzimidazole family", [("R1", "benzimidazole substituent"), ("Ar", "aryl group")]),
    ("lit1b_cyclic_product", "cyclized benzimidazole family", [("R1", "benzimidazole substituent"), ("Ar", "aryl group")]),
    ("generic_1", "N-tethered terminal-allene benzimidazole family", [("R1", "benzimidazole substituent"), ("R2", "allene-terminal substituent")]),
    ("generic_2", "cyclized allylated benzimidazole family", [("R1", "benzimidazole substituent"), ("R2", "new stereocenter substituent")]),
    ("generic_3", "N-substituted benzimidazole family", [("R1", "benzimidazole substituent")]),
    ("generic_4", "terminal allene family", [("R2", "allene substituent")]),
    ("generic_5", "C2-allylated benzimidazole E/Z family", [("R1", "benzimidazole substituent"), ("R2", "allyl substituent")]),
]
for label, fixed, variables in generic_specs:
    molecules.append(generic_molecule(label, E[1], fixed, variables,
                                      "Scheme 1 uses variable substituents; a wildcard SMILES is intentionally not used."))

# Specific cyclization inputs/products. Product stereocenters are not assigned in the source.
for label in ["1a", "1b", "1c"]:
    molecules.append(molecule(label, E[2] + (E[3] if label == "1c" else []),
                              chirality="achiral", geometry="not_applicable", representation="not_encoded",
                              axial=[{"axis": "terminal allene", "assignment": "not axially stereogenic; CH2 terminus"}],
                              note="Source-depicted terminal allene; complete graph retained as unresolved."))
for letter in "defghijklmnopqrst":
    label = "1" + letter
    molecules.append(molecule(label, E[2], chirality="achiral", geometry="not_applicable",
                              representation="not_encoded",
                              axial=[{"axis": "terminal allene", "assignment": "not axially stereogenic; CH2 terminus"}],
                              note="Scheme 2 substrate inferred by one-to-one correspondence with the labeled product; full graph not encoded."))
for label in ["2a", "2b", "2b′"] + ["2" + x for x in "cdefghijklmnopqrst"]:
    molecules.append(molecule(label, E[2] + (E[3] if label in ("2c", "2o") else []),
                              chirality="unspecified", geometry="not_applicable", representation="not_encoded",
                              note="Cyclized product drawing contains a potentially stereogenic ring carbon but no wedge/dash or absolute assignment."))

# Intermolecular scope substrate and product labels.
for letter in "bcdefghijklmnop":
    molecules.append(molecule("3" + letter, E[3], chirality="achiral", geometry="not_applicable",
                              representation="not_encoded",
                              note="Source-depicted substituted benzimidazole; no stereocenter or geometry assignment is shown."))
for label in ["4d", "4e"]:
    molecules.append(molecule(label, E[3], chirality="achiral", geometry="not_applicable",
                              representation="not_encoded",
                              axial=[{"axis": "terminal allene", "assignment": "not axially stereogenic; CH2 terminus"}],
                              note="Terminal allene structure is source-depicted but not fully encoded."))
molecules.append(molecule("4f", E[3], chirality="unspecified", geometry="not_applicable",
                          representation="not_encoded",
                          axial=[{"axis": "internal allene", "assignment": "unspecified"}],
                          note="Internal allene may carry axial stereochemistry; the source gives no axial assignment."))
product_labels = ["5" + x + "a" for x in "abcdefghijklmnop"] + ["5a" + x for x in "bcdef"]
success_ratios = {"5aa": "15:1", "5ba": ">20:1", "5ca": ">20:1", "5da": ">20:1",
                  "5ea": ">20:1", "5fa": ">20:1", "5ga": ">20:1", "5ha": "9:1",
                  "5ia": "14:1", "5ja": "13:1", "5ka": "15:1", "5ab": ">20:1",
                  "5ac": ">20:1", "5ad": ">20:1", "5ae": ">20:1"}
for label in product_labels:
    molecules.append(molecule(label, E[2] + E[3] if label == "5aa" else E[3],
                              chirality="unresolved", geometry="mixture" if label in success_ratios else "unspecified",
                              representation="not_encoded",
                              mixture={"kind": "E/Z mixture", "ratio": "experiment-specific"} if label in success_ratios else None,
                              note="Allylated product retained as the source-labeled E/Z material; individual isomer graphs are not assigned."))

# Deuterated species and proposed intermediates.
for label, note, chiral, geom in [
    ("1c-d", "C2-deuterated 1c; 83% D in Scheme 4.", "achiral", "not_applicable"),
    ("2c-d", "Deuterated cyclization product; incorporation varies by experiment.", "unspecified", "not_applicable"),
    ("2o-d", "Deuterated cyclization product from 1o competition.", "unspecified", "not_applicable"),
    ("3a-d", "C2-deuterated 3a; 90% D in Scheme 4.", "achiral", "not_applicable"),
    ("5aa-d", "Deuterated intermolecular product; E/Z and D distribution not encoded.", "unresolved", "unspecified"),
    ("5fa-d", "Deuterated intermolecular product; E/Z and D distribution not encoded.", "unresolved", "unspecified"),
    ("LCuH_A", "Proposed LCuH intermediate A, not an isolated compound.", "unresolved", "unresolved"),
    ("allylCu_B", "Proposed allylcopper intermediate B, not an isolated compound.", "unresolved", "unresolved"),
    ("copper_nitride_C", "Proposed copper nitride intermediate C, not an isolated compound.", "unspecified", "unresolved"),
    ("silylated_D", "Proposed silylated intermediate D, not an isolated compound.", "unspecified", "unresolved"),
]:
    molecules.append(molecule(label, E[3] if not label.endswith(("_A", "_B", "_C", "_D")) else E[4],
                              chirality=chiral, geometry=geom, representation="not_encoded", note=note))

bundle["molecules"] = molecules
reactions = []

# Page 1: cited/prior work and the two generic current transformations.
lit_records = [
    ("scheme1a_buchwald", "lit1a_substrate", "lit1a_product", "MTBE",
     [("Cu(OAc)2", "catalyst", 5), ("(S,S)-Ph-BPE", "ligand", 6)],
     [("(MeO)2MeSiH", "reagent", 4)], "Buchwald's work"),
    ("scheme1b_xu_acyclic", "lit1b_substrate", "lit1b_product", "toluene",
     [("CuOAc", "catalyst", 10), ("DCyPE", "ligand", 12)],
     [("PMHS", "reagent", 5)], "Xu's work, acyclic example"),
    ("scheme1b_xu_cyclization", "lit1b_tethered_substrate", "lit1b_cyclic_product", "toluene",
     [("CuOAc", "catalyst", 10), ("DCyPE", "ligand", 12)],
     [("(MeO)2MeSiH", "reagent", 5)], "Xu's work, tethered example"),
]
for suffix, sub, prod, solvent, cats, reagents, loc in lit_records:
    r = reaction(suffix, E[1], "literature_example", ["other"], {"page": 1, "scheme": "1", "region": loc},
                 "CuH-catalyzed benzimidazole functionalization", "C2 hydrofunctionalization")
    r["participants"] = [participant("p_sub", sub, "reactant", E[1]), participant("p_prod", prod, "product", E[1], product_role="main_product")]
    for label, role, loading in cats:
        r["participants"].append(participant("p_" + label, label, role, E[1], mol_percent=(f"{loading} mol %", quantity(loading, "mol%")), ref=sub))
    for label, role, equiv in reagents:
        r["participants"].append(participant("p_" + label, label, role, E[1], equivalents=(f"{equiv} equiv", quantity(equiv, "equiv")), ref=sub))
    add_conditions(r, E[1], solvent)
    classification(r, E[1], True)
    r["review"]["unresolved"].append("Cited literature example: yield, temperature, time, and detailed procedure are not reported in this source.")
    reactions.append(r)

generic_ids = {}
for suffix, sub, prod, solvent, kind in [
    ("scheme1c_cyclization_generic", "generic_1", "generic_2", "THF", "intramolecular"),
    ("scheme1c_intermolecular_generic", "generic_3", "generic_5", "toluene", "intermolecular"),
]:
    r = reaction(suffix, E[1], "generic_scheme", ["other"], {"page": 1, "scheme": "1c", "mode": kind},
                 "CuH-catalyzed C-H allylation of benzimidazoles with allenes", f"{kind} C2 allylation")
    r["participants"] = [participant("p_sub", sub, "reactant", E[1]), participant("p_prod", prod, "product", E[1], product_role="main_product")]
    if kind == "intermolecular":
        r["participants"].insert(1, participant("p_allene", "generic_4", "reactant", E[1]))
    for label, role, loading in [("Cu(OAc)2", "catalyst", 5), ("Xantphos", "ligand", 6), ("(MeO)2MeSiH", "reagent", 50)]:
        r["participants"].append(participant("p_" + label, label, role, E[1], mol_percent=(f"{loading} mol %", quantity(loading, "mol%")), ref=sub))
    add_conditions(r, E[1], solvent)
    classification(r, E[1], True)
    r["notes"].append(f("The page-1 graphical abstract repeats this current-work transformation.",
                        "duplicate depiction; not a second experiment", E[1]))
    reactions.append(r)
    generic_ids[kind] = r["reaction_id"]

# Table 1 and Scheme 2.
table1 = [
    (1, "1a", 200, "DCyPE", "2a", 0, None),
    (2, "1a", 200, "Xantphos", "2a", 0, None),
    (3, "1b", 200, "DCyPE", "2b", 0, None),
    (4, "1b", 200, "Xantphos", "2b", 54, ("2b′", 16)),
    (5, "1c", 200, "Xantphos", "2c", ">99", None),
    (6, "1c", 100, "Xantphos", "2c", ">99", None),
    (7, "1c", 50, "Xantphos", "2c", ">99", None),
    (8, "1c", 40, "Xantphos", "2c", 74, None),
    (9, "1c", 30, "Xantphos", "2c", "trace", None),
    (10, "1c", 50, "Xantphos", "2c", 61, None),
]
table1_ids = {}
for entry, sub, x, ligand, prod, yld, byproduct in table1:
    suffix = f"table1_entry_{entry:02d}" + ("_scheme2_2c" if entry == 7 else "")
    r = reaction(suffix, E[2], "specific_experiment", ["optimization"],
                 {"page": 2, "table": 1, "entry": entry, "duplicate_scope_label": "2c" if entry == 7 else None},
                 "CuH-catalyzed cyclization of benzimidazoles", "intramolecular C2 allylation",
                 [generic_ids["intramolecular"]])
    r["participants"] = [
        participant("p_sub", sub, "reactant", E[2], amount=("0.4 mmol", quantity(0.4, "mmol"))),
        participant("p_prod", prod, "product", E[2], product_role="main_product"),
        participant("p_cu", "Cu(OAc)2", "catalyst", E[2], mol_percent=("5 mol %", quantity(5, "mol%")), ref=sub),
        participant("p_lig", ligand, "ligand", E[2], mol_percent=("6 mol %", quantity(6, "mol%")), ref=sub),
        participant("p_silane", "(EtO)2MeSiH" if entry == 10 else "(MeO)2MeSiH", "reagent", E[2],
                    mol_percent=(f"{x} mol %", quantity(x, "mol%")), ref=sub),
    ]
    if byproduct:
        r["participants"].append(participant("p_byproduct", byproduct[0], "product", E[2], product_role="byproduct"))
    add_conditions(r, E[2], "THF", time_h=20, atmosphere="Ar", volume_ml=2, concentration_m=0.5)
    r["reaction_scale"] = f("compound 1 (0.4 mmol)", {"reference_molecule_id": mid(sub), "amount": quantity(0.4, "mmol")}, E[2])
    raw = "trace" if yld == "trace" else f"{yld}%"
    val = "trace" if yld == "trace" else quantity(yld if isinstance(yld, int) else 99, "%", qualifier="gt" if yld == ">99" else None, raw_text=raw)
    add_yield(r, "yield_main", prod, raw, val, E[2], isolated=True, scope=f"Table 1 entry {entry}")
    if byproduct:
        add_yield(r, "yield_byproduct", byproduct[0], f"{byproduct[1]}%", quantity(byproduct[1], "%"), E[2], isolated=True,
                  role="byproduct", scope="Table 1 footnote c")
    observed = yld not in (0,)
    classification(r, E[2], observed)
    r["notes"].append(inferred(
        ["THF (0.5 M)", "compound 1 (0.4 mmol), THF (2 mL)"],
        "The printed concentration basis is not stated; the arrow and footnote are retained without arithmetic correction.",
        E[2], "0.5 M cannot be safely reconciled to a unique basis from the main-text quantities alone."
    ))
    r["review"]["unresolved"].append("THF is printed as 0.5 M while the footnote gives 0.4 mmol substrate and 2 mL THF; concentration basis is unreported.")
    if entry == 7:
        r["notes"].append(f("Scheme 2 product 2c repeats Table 1 entry 7 (>99%).",
                            "same experiment represented once with duplicate location noted", E[2]))
    reactions.append(r)
    table1_ids[entry] = r["reaction_id"]

scheme2 = {"2c": ">99", "2d": 91, "2e": 96, "2f": 92, "2g": 86, "2h": 93,
           "2i": 81, "2j": 90, "2k": 74, "2l": 79, "2m": 91, "2n": 77,
           "2o": 84, "2p": 69, "2q": 92, "2r": 0, "2s": 0, "2t": 0}
scheme2_ids = {"2c": table1_ids[7]}
for prod, yld in scheme2.items():
    if prod == "2c":
        continue
    sub = "1" + prod[1:]
    r = reaction(f"scheme2_{prod}", E[2], "specific_experiment", ["substrate_scope"],
                 {"page": 2, "scheme": 2, "product": prod},
                 "CuH-catalyzed allylation/cyclization of benzimidazoles", "intramolecular C2 allylation",
                 [generic_ids["intramolecular"]])
    r["participants"] = [
        participant("p_sub", sub, "reactant", E[2], amount=("0.4 mmol", quantity(0.4, "mmol"))),
        participant("p_prod", prod, "product", E[2], product_role="main_product"),
        participant("p_cu", "Cu(OAc)2", "catalyst", E[2], mol_percent=("5 mol %", quantity(5, "mol%")), ref=sub),
        participant("p_lig", "Xantphos", "ligand", E[2], mol_percent=("6 mol %", quantity(6, "mol%")), ref=sub),
        participant("p_silane", "(MeO)2MeSiH", "reagent", E[2], mol_percent=("50 mol %", quantity(50, "mol%")), ref=sub),
    ]
    add_conditions(r, E[2], "THF", time_h=20, atmosphere="Ar", volume_ml=2, concentration_m=0.5, inherited=True)
    r["reaction_scale"] = f("conditions in entry 7 of Table 1: compound 1 (0.4 mmol)",
                            {"reference_molecule_id": mid(sub), "amount": quantity(0.4, "mmol")}, E[2],
                            method="inherited_general_condition")
    raw = f"{yld}%"
    add_yield(r, "yield", prod, raw, quantity(yld, "%"), E[2], isolated=True, scope=f"Scheme 2 {prod}")
    classification(r, E[2], yld != 0)
    r["review"]["unresolved"].append("Full substrate/product SMILES and any unreported stereoisomer distribution require SI/source-structure review.")
    reactions.append(r)
    scheme2_ids[prod] = r["reaction_id"]

# Table 2 and Scheme 3.
table2 = [(1, "THF", 73, "4:1"), (2, "DCE", 31, None), (3, "CH3CN", "trace", None),
          (4, "toluene", 93, "15:1"), (5, "1,4-dioxane", 0, None), (6, "MTBE", 0, None)]
table2_ids = {}
for entry, solvent, yld, ratio in table2:
    suffix = f"table2_entry_{entry:02d}" + ("_scheme3_5aa" if entry == 4 else "")
    ev = E[2] + (E[3] if entry == 4 else [])
    r = reaction(suffix, ev, "specific_experiment", ["optimization"],
                 {"page": 2, "table": 2, "entry": entry, "duplicate_scope_label": "5aa" if entry == 4 else None},
                 "CuH-catalyzed allylation of benzimidazoles", "intermolecular C2 allylation",
                 [generic_ids["intermolecular"]])
    r["participants"] = [
        participant("p_3a", "3a", "reactant", ev, amount=("0.4 mmol", quantity(0.4, "mmol"))),
        participant("p_4a", "4a", "reactant", ev, amount=("0.6 mmol", quantity(0.6, "mmol"))),
        participant("p_prod", "5aa", "product", ev, product_role="main_product"),
        participant("p_cu", "Cu(OAc)2", "catalyst", ev, mol_percent=("5 mol %", quantity(5, "mol%")), ref="3a"),
        participant("p_lig", "Xantphos", "ligand", ev, mol_percent=("6 mol %", quantity(6, "mol%")), ref="3a"),
        participant("p_silane", "(MeO)2MeSiH", "reagent", ev, mol_percent=("50 mol %", quantity(50, "mol%")), ref="3a"),
    ]
    add_conditions(r, ev, solvent, time_h=20, atmosphere="Ar", volume_ml=2, concentration_m=0.5)
    r["reaction_scale"] = f("compound 3a (0.4 mmol)", {"reference_molecule_id": mid("3a"), "amount": quantity(0.4, "mmol")}, ev)
    raw = "trace" if yld == "trace" else f"{yld}%"
    val = "trace" if yld == "trace" else quantity(yld, "%")
    add_yield(r, "yield", "5aa", raw, val, ev, isolated=True, scope=f"Table 2 entry {entry}")
    if ratio:
        add_selectivity(r, "ez", "5aa", ratio, ev, "E/Z isomers isolated and determined by 1H NMR spectroscopy")
    classification(r, ev, yld != 0)
    if entry == 4:
        r["notes"].append(f("Scheme 3 product 5aa repeats Table 2 entry 4 (93%, 15:1 E/Z).",
                            "same experiment represented once with both page-2 and page-3 evidence", ev))
    reactions.append(r)
    table2_ids[entry] = r["reaction_id"]

scheme3 = {
    "5aa": (93, "15:1"), "5ba": (81, ">20:1"), "5ca": (84, ">20:1"),
    "5da": (76, ">20:1"), "5ea": (69, ">20:1"), "5fa": (68, ">20:1"),
    "5ga": (89, ">20:1"), "5ha": (76, "9:1"), "5ia": (75, "14:1"),
    "5ja": (61, "13:1"), "5ka": (64, "15:1"), "5la": ("trace", None),
    "5ma": (0, None), "5na": (0, None), "5oa": (0, None), "5pa": (0, None),
    "5ab": (81, ">20:1"), "5ac": (63, ">20:1"), "5ad": (68, ">20:1"),
    "5ae": (53, ">20:1"), "5af": (0, None),
}
scheme3_ids = {"5aa": table2_ids[4]}
for prod, (yld, ratio) in scheme3.items():
    if prod == "5aa":
        continue
    if prod[1] == "a":
        sub1, sub2 = "3a", "4" + prod[2]
    else:
        sub1, sub2 = "3" + prod[1], "4a"
    r = reaction(f"scheme3_{prod}", E[3], "specific_experiment", ["substrate_scope"],
                 {"page": 3, "scheme": 3, "product": prod},
                 "CuH-catalyzed allylation of benzimidazoles", "intermolecular C2 allylation",
                 [generic_ids["intermolecular"]])
    r["participants"] = [
        participant("p_bzim", sub1, "reactant", E[3], amount=("0.4 mmol", quantity(0.4, "mmol"))),
        participant("p_allene", sub2, "reactant", E[3], amount=("0.6 mmol", quantity(0.6, "mmol"))),
        participant("p_prod", prod, "product", E[3], product_role="main_product"),
        participant("p_cu", "Cu(OAc)2", "catalyst", E[3], mol_percent=("5 mol %", quantity(5, "mol%")), ref=sub1),
        participant("p_lig", "Xantphos", "ligand", E[3], mol_percent=("6 mol %", quantity(6, "mol%")), ref=sub1),
        participant("p_silane", "(MeO)2MeSiH", "reagent", E[3], mol_percent=("50 mol %", quantity(50, "mol%")), ref=sub1),
    ]
    add_conditions(r, E[3], "toluene", time_h=20, atmosphere="Ar", volume_ml=2, concentration_m=0.5, inherited=True)
    r["reaction_scale"] = f("conditions in entry 4 of Table 2: compound 3 (0.4 mmol)",
                            {"reference_molecule_id": mid(sub1), "amount": quantity(0.4, "mmol")}, E[3],
                            method="inherited_general_condition")
    raw = "trace" if yld == "trace" else f"{yld}%"
    val = "trace" if yld == "trace" else quantity(yld, "%")
    add_yield(r, "yield", prod, raw, val, E[3], isolated=True, scope=f"Scheme 3 {prod}")
    if ratio:
        add_selectivity(r, "ez", prod, ratio, E[3], "E/Z isomers determined by crude 1H NMR spectroscopy")
    classification(r, E[3], yld != 0)
    reactions.append(r)
    scheme3_ids[prod] = r["reaction_id"]

# Scheme 4: four real mechanistic experiments.
def mech_base(suffix, location, ev, inter=False):
    return reaction(suffix, ev, "specific_experiment", ["mechanistic_study"], location,
                    "deuterium-labeling/competition experiment", "intermolecular C2 allylation" if inter else "intramolecular C2 allylation",
                    [generic_ids["intermolecular" if inter else "intramolecular"]])

r = mech_base("scheme4_a_deuterium_cyclization", {"page": 3, "scheme": 4, "equation": "a"}, E[3])
r["participants"] = [
    participant("p_sub", "1c-d", "reactant", E[3], amount=("0.4 mmol", quantity(0.4, "mmol"))),
    participant("p_prod", "2c-d", "product", E[3], product_role="main_product"),
    participant("p_cu", "Cu(OAc)2", "catalyst", E[3], mol_percent=("5 mol %", quantity(5, "mol%")), ref="1c-d"),
    participant("p_lig", "Xantphos", "ligand", E[3], mol_percent=("6 mol %", quantity(6, "mol%")), ref="1c-d"),
    participant("p_silane", "(MeO)2MeSiH", "reagent", E[3], mol_percent=("50 mol %", quantity(50, "mol%")), ref="1c-d"),
]
add_conditions(r, E[3], "THF", time_h=20, atmosphere="Ar", concentration_m=0.5)
add_yield(r, "yield", "2c-d", "87% yield", quantity(87, "%"), E[3], isolated=False, scope="Scheme 4a")
add_isotope_outcome(r, "d_incorporation", "2c-d", 68, E[3])
r["notes"].append(f("1c-d: 83% D", {"starting_material_D_incorporation": quantity(83, "%")}, E[3]))
classification(r, E[3], True)
reactions.append(r)

r = mech_base("scheme4_b_deuterium_competition", {"page": 3, "scheme": 4, "equation": "b"}, E[3])
r["participants"] = [
    participant("p_1cd", "1c-d", "reactant", E[3], amount=("0.2 mmol", quantity(0.2, "mmol"))),
    participant("p_1o", "1o", "reactant", E[3], amount=("0.2 mmol", quantity(0.2, "mmol"))),
    participant("p_2cd", "2c-d", "product", E[3], product_role="main_product"),
    participant("p_2od", "2o-d", "product", E[3], product_role="main_product"),
    participant("p_cu", "Cu(OAc)2", "catalyst", E[3], mol_percent=("5 mol %", quantity(5, "mol%")), ref="1c-d"),
    participant("p_lig", "Xantphos", "ligand", E[3], mol_percent=("6 mol %", quantity(6, "mol%")), ref="1c-d"),
    participant("p_silane", "(MeO)2MeSiH", "reagent", E[3], mol_percent=("50 mol %", quantity(50, "mol%")), ref="1c-d"),
]
add_conditions(r, E[3], "THF", time_h=20, atmosphere="Ar", concentration_m=0.5)
add_yield(r, "yield_2cd", "2c-d", "84% yield", quantity(84, "%"), E[3], isolated=False, scope="Scheme 4b")
add_yield(r, "yield_2od", "2o-d", "58% yield", quantity(58, "%"), E[3], isolated=False, scope="Scheme 4b")
add_isotope_outcome(r, "d_2cd", "2c-d", 30, E[3]); add_isotope_outcome(r, "d_2od", "2o-d", 28, E[3])
r["notes"].append(f("1c-d: 83% D", {"starting_material_D_incorporation": quantity(83, "%")}, E[3]))
classification(r, E[3], True)
reactions.append(r)

r = mech_base("scheme4_c_deuterium_allylation", {"page": 3, "scheme": 4, "equation": "c"}, E[3], inter=True)
r["participants"] = [
    participant("p_3ad", "3a-d", "reactant", E[3], amount=("0.4 mmol", quantity(0.4, "mmol"))),
    participant("p_4a", "4a", "reactant", E[3], amount=("0.6 mmol", quantity(0.6, "mmol"))),
    participant("p_prod", "5aa-d", "product", E[3], product_role="main_product"),
    participant("p_cu", "Cu(OAc)2", "catalyst", E[3], mol_percent=("5 mol %", quantity(5, "mol%")), ref="3a-d"),
    participant("p_lig", "Xantphos", "ligand", E[3], mol_percent=("6 mol %", quantity(6, "mol%")), ref="3a-d"),
    participant("p_silane", "(MeO)2MeSiH", "reagent", E[3], mol_percent=("50 mol %", quantity(50, "mol%")), ref="3a-d"),
]
add_conditions(r, E[3], "toluene", time_h=20, atmosphere="Ar", concentration_m=0.5)
add_ambiguous_yield(r, "printed_83", "5aa-d", "83 (printed beneath 5aa-d; no unit or yield label)", E[3])
add_isotope_outcome(r, "d_incorporation", "5aa-d", 64, E[3])
r["notes"].append(f("3a-d: 90% D", {"starting_material_D_incorporation": quantity(90, "%")}, E[3]))
r["review"]["unresolved"].append("The numeral 83 under 5aa-d has no printed unit/yield label in the supplied main-text figure and is not normalized to percent.")
classification(r, E[3], True)
reactions.append(r)

r = mech_base("scheme4_d_deuterium_competition", {"page": 3, "scheme": 4, "equation": "d"}, E[3], inter=True)
r["participants"] = [
    participant("p_3ad", "3a-d", "reactant", E[3], amount=("0.2 mmol", quantity(0.2, "mmol"))),
    participant("p_3f", "3f", "reactant", E[3], amount=("0.2 mmol", quantity(0.2, "mmol"))),
    participant("p_4a", "4a", "reactant", E[3], amount=("0.6 mmol", quantity(0.6, "mmol"))),
    participant("p_5aad", "5aa-d", "product", E[3], product_role="main_product"),
    participant("p_5fad", "5fa-d", "product", E[3], product_role="main_product"),
    participant("p_cu", "Cu(OAc)2", "catalyst", E[3], mol_percent=("5 mol %", quantity(5, "mol%")), ref="3a-d"),
    participant("p_lig", "Xantphos", "ligand", E[3], mol_percent=("6 mol %", quantity(6, "mol%")), ref="3a-d"),
    participant("p_silane", "(MeO)2MeSiH", "reagent", E[3], mol_percent=("50 mol %", quantity(50, "mol%")), ref="3a-d"),
]
add_conditions(r, E[3], "toluene", time_h=20, atmosphere="Ar", concentration_m=0.5)
add_ambiguous_yield(r, "printed_72", "5aa-d", "72 (printed beneath 5aa-d; no unit or yield label)", E[3])
add_ambiguous_yield(r, "printed_47", "5fa-d", "47 (printed beneath 5fa-d; no unit or yield label)", E[3])
add_isotope_outcome(r, "d_5aad", "5aa-d", 26, E[3]); add_isotope_outcome(r, "d_5fad", "5fa-d", 28, E[3])
r["notes"].append(f("3a-d: 90% D", {"starting_material_D_incorporation": quantity(90, "%")}, E[3]))
r["review"]["unresolved"].append("The numerals 72 and 47 under products lack printed units/yield labels and are not normalized to percent.")
classification(r, E[3], True)
reactions.append(r)

# Page 4: proposed catalytic cycle and prose-only E/Z transition-state proposal.
r = reaction("scheme5_proposed_cyclization_mechanism", E[4], "proposed_mechanism", ["mechanistic_study"],
             {"page": 4, "scheme": 5}, "Proposed mechanism for CuH-catalyzed cyclization of benzimidazoles with allenes",
             "proposed CuH hydrocupration/cyclization/beta-hydride-elimination cycle", [generic_ids["intramolecular"]])
r["participants"] = [
    participant("p_1c", "1c", "reactant", E[4]), participant("p_2c", "2c", "product", E[4], product_role="main_product"),
    participant("p_A", "LCuH_A", "catalyst", E[4]), participant("p_B", "allylCu_B", "unknown", E[4]),
    participant("p_C", "copper_nitride_C", "unknown", E[4]), participant("p_D", "silylated_D", "unknown", E[4]),
    participant("p_cu", "Cu(OAc)2", "catalyst", E[4]), participant("p_lig", "Xantphos", "ligand", E[4]),
    participant("p_silane", "(MeO)2MeSiH", "reagent", E[4]),
]
ops = [
    ("formation_A", ["p_cu", "p_lig", "p_silane", "p_A"], "proposed formation of LCuH A"),
    ("hydrocupration", ["p_A", "p_1c", "p_B"], "proposed hydrocupration of allene 1c to allylcopper B"),
    ("cyclization", ["p_B", "p_C"], "proposed intramolecular addition to C=N giving C"),
    ("path_a", ["p_C", "p_A", "p_2c"], "proposed beta-hydride elimination path a"),
    ("path_b_silylation", ["p_C", "p_silane", "p_D", "p_A"], "proposed silane interaction path b"),
    ("path_b_product", ["p_D", "p_2c"], "proposed net loss of silane giving 2c"),
]
r["steps"] = []
for sid, pids, op in ops:
    s = st.new_step(sid, E[4]); s["participant_ids"] = pids; s["operation"] = f(op, op, E[4]); r["steps"].append(s)
classification(r, E[4], True)
r["review"]["unresolved"].append("All intermediates and elementary arrows are author-proposed; no isolation or direct observation is claimed.")
reactions.append(r)

r = reaction("page4_proposed_inter_e_selectivity", E[4], "proposed_mechanism", ["mechanistic_study"],
             {"page": 4, "region": "left-column prose below Scheme 5"},
             "Proposed six-membered transition-state explanation of E/Z selectivity",
             "proposed stereochemical model for intermolecular C2 allylation", [generic_ids["intermolecular"]])
r["participants"] = [participant("p_3a", "3a", "reactant", E[4]), participant("p_4a", "4a", "reactant", E[4]),
                     participant("p_5aa", "5aa", "product", E[4], product_role="E-major product mixture")]
s = st.new_step("proposed_transition_state_comparison", E[4]); s["participant_ids"] = [p["participant_id"] for p in r["participants"]]
s["operation"] = f("TSA places cyclohexyl exocyclic axial; TSB places it endocyclic axial with steric repulsion",
                   "author-proposed transition-state comparison", E[4]); r["steps"] = [s]
classification(r, E[4], True)
r["review"]["unresolved"].append("Transition-state structures are described in prose and referred to the unavailable SI; no computed energies are reported in the supplied pages.")
reactions.append(r)

bundle["reactions"] = reactions
bundle["source_context"] = [
    {"context_id": f"{DOC}:article_metadata", "evidence_ids": E[1],
     "title": "Cu-Catalyzed C-H Allylation of Benzimidazoles with Allenes", "doi": "10.1021/acs.orglett.1c02346",
     "authors": ["Yaxi Dong", "Bernhard Breit"]},
    {"context_id": f"{DOC}:mechanistic_interpretation", "evidence_ids": E[3] + E[4],
     "type": "author_interpretation", "summary": "Deuterium scrambling is interpreted as participation of the benzimidazole C2 hydrogen in CuH formation; this interpretation is not stored as direct proof of the proposed cycle."},
    {"context_id": f"{DOC}:si_notice", "evidence_ids": E[4], "type": "supplement_notice",
     "summary": "The article advertises Supporting Information, but no SI file was supplied or explicitly associated in this workspace."},
    {"context_id": f"{DOC}:page5_references", "evidence_ids": E[5], "type": "bibliography",
     "summary": "Page 5 contains references only; cited works were not converted into current-paper experiments beyond the explicitly depicted Scheme 1 literature examples."},
    {"context_id": f"{DOC}:characterization_screen", "evidence_ids": E[1] + E[2] + E[3] + E[4] + E[5],
     "type": "coverage_note", "summary": "No standalone mass-spectrum peak table, molecular-weight table, or fragment relative-abundance table occurs on pages 1-5."},
]

# Raw transcription and source-mention hints are saved separately from the formal bundle.
transcription = {
    "document_id": DOC,
    "title": bundle["source_context"][0]["title"],
    "table1": [{"entry": e, "substrate": s, "silane_mol_percent": x, "ligand": l, "product": p, "yield": y, "byproduct": b}
               for e, s, x, l, p, y, b in table1],
    "scheme2": scheme2,
    "table2": [{"entry": e, "solvent": s, "yield": y, "E/Z": z} for e, s, y, z in table2],
    "scheme3": {k: {"yield": v[0], "E/Z": v[1]} for k, v in scheme3.items()},
    "scheme4": {
        "a": {"input": "1c-d, 0.4 mmol, 83% D", "product": "2c-d, 87% yield, 68% D"},
        "b": {"inputs": ["1c-d, 0.2 mmol, 83% D", "1o, 0.2 mmol"], "products": ["2c-d, 84% yield, 30% D", "2o-d, 58% yield, 28% D"]},
        "c": {"inputs": ["3a-d, 0.4 mmol, 90% D", "4a, 0.6 mmol"], "product": "5aa-d, printed 83 (unit absent), 64% D"},
        "d": {"inputs": ["3a-d, 0.2 mmol, 90% D", "3f, 0.2 mmol", "4a, 0.6 mmol"], "products": ["5aa-d, printed 72 (unit absent), 26% D", "5fa-d, printed 47 (unit absent), 28% D"]},
    },
    "duplicate_experiment_decisions": {
        "Table 1 entry 7 / Scheme 2 2c": table1_ids[7],
        "Table 2 entry 4 / Scheme 3 5aa": table2_ids[4],
    },
}
st.write(RUN / "raw_model/transcription.json", transcription)
st.write(RUN / "raw_model/stereo_diagnostics.json", diagnostics)
mentions = st.source_mentions(json.loads((RUN / "dataset.json").read_text()),
                              ["1c", "2c", "2b", "3a", "4a", "5aa", "1c-d", "3a-d"])
st.write(RUN / "raw_model/source_mentions.json", mentions)

# Exact five-page coverage required for this session.
def evidence_pages(record):
    return {int(eid.split(":page")[1].split(":")[0])
            for eid in record["evidence_ids"] if ":page" in eid}


p1_ids = [r["reaction_id"] for r in reactions if 1 in evidence_pages(r)]
p2_ids = [r["reaction_id"] for r in reactions if 2 in evidence_pages(r)]
p3_ids = [r["reaction_id"] for r in reactions if 3 in evidence_pages(r)]
p4_ids = [r["reaction_id"] for r in reactions if 4 in evidence_pages(r)]
coverage = {
    "document_id": DOC, "semantic_accuracy_verified": False,
    "pages": [
        {"page": 1, "state": "attempted", "regions": [
            {"region": "title/abstract and graphical abstract", "finding": "Article metadata and duplicate overview of the two current-work transformations; not separate experiments.", "evidence_ids": E[1]},
            {"region": "Scheme 1a", "finding": "One cited Buchwald literature example.", "evidence_ids": E[1]},
            {"region": "Scheme 1b", "finding": "Two cited Xu literature arrows (acyclic and tethered examples).", "evidence_ids": E[1]},
            {"region": "Scheme 1c", "finding": "Two generic current-work arrows (intra- and intermolecular); no specific yield on this page.", "evidence_ids": E[1]},
        ], "record_ids": p1_ids, "unresolved": ["Scheme 1 literature-example yields, temperatures, and durations are not reported here."]},
        {"page": 2, "state": "attempted", "regions": [
            {"region": "Table 1 including reaction arrow and footnotes a-d", "finding": "Ten optimization experiments; entry 4 also reports byproduct 2b′ at 16%; entry 10 replaces the hydrosilane.", "evidence_ids": E[2]},
            {"region": "Scheme 2 including footnotes", "finding": "Eighteen scope outcomes 2c-2t; 2c duplicates Table 1 entry 7 and is not double-counted.", "evidence_ids": E[2]},
            {"region": "Table 2 including reaction arrow and footnotes a-c", "finding": "Six solvent-optimization experiments with two E/Z ratios; entry 4 is repeated as Scheme 3 product 5aa.", "evidence_ids": E[2]},
        ], "record_ids": p2_ids, "unresolved": ["No SI procedures/workups were supplied.", "Table 1 prints THF 0.5 M and also 0.4 mmol substrate in 2 mL THF; concentration basis is not stated."]},
        {"page": 3, "state": "attempted", "regions": [
            {"region": "Scheme 3 including footnotes a-c", "finding": "Twenty-one intermolecular scope outcomes; 5aa duplicates Table 2 entry 4; E/Z ratios use crude 1H NMR.", "evidence_ids": E[3]},
            {"region": "Scheme 4 equations a-d", "finding": "Four real deuterium-labeling/competition experiments, kept separate from proposed mechanisms.", "evidence_ids": E[3]},
            {"region": "mechanistic discussion below Scheme 4", "finding": "Author interpretation of deuterium scrambling and lead-in to Scheme 5; not an additional experiment.", "evidence_ids": E[3]},
        ], "record_ids": p3_ids, "unresolved": ["Scheme 4c/d numerals beneath products omit units/yield labels and remain unnormalized.", "Scheme 4 yield measurement bases are not stated."]},
        {"page": 4, "state": "attempted", "regions": [
            {"region": "Scheme 5 arrows and intermediates A-D", "finding": "Author-proposed catalytic cycle; one proposed-mechanism record with six mapped arrow steps, not performed experiments.", "evidence_ids": E[4]},
            {"region": "left-column prose below Scheme 5", "finding": "Proposed TSA/TSB explanation of E-selectivity; SI details unavailable.", "evidence_ids": E[4]},
            {"region": "Associated Content", "finding": "SI is advertised but not supplied or explicitly associated.", "evidence_ids": E[4]},
        ], "record_ids": p4_ids, "unresolved": ["Proposed intermediates/transition states are not isolated or fully encoded.", "Referenced SI mechanism details are unavailable."]},
        {"page": 5, "state": "attempted", "regions": [
            {"region": "references 8-12 and remainder of page", "finding": "Bibliography only; no current-paper experimental arrow, table, footnote, or characterization peak table.", "evidence_ids": E[5]},
        ], "record_ids": [], "unresolved": []},
    ]
}
st.write(RUN / "coverage/pages_0001_0005.json", coverage)

tasks_path = RUN / "tasks.json"
tasks = json.loads(tasks_path.read_text())
for task in tasks:
    if task.get("document_id") == DOC and task.get("evidence_id", "").endswith(":text"):
        task["state"] = "attempted"
        task["reason"] = "Original page image and extracted text inspected; coverage decision saved in pages_0001_0005.json."
st.write(tasks_path, tasks)
st.write(RUN / "submitted.json", bundle)
print(json.dumps({"molecules": len(molecules), "reactions": len(reactions),
                  "stereo_diagnostics": len(diagnostics),
                  "coverage_pages": len(coverage["pages"])}, ensure_ascii=False, indent=2))
