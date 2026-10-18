#!/usr/bin/env python3
"""Build the page 1-6 v3.1 extraction bundle from a manually inspected inventory.

This is a deterministic transcription aid.  It does not call a model, infer SI
association, or claim that the machine structures have been source-verified.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
SKILL = ROOT / "skills" / "chemical-reaction-extraction"
sys.path.insert(0, str(SKILL / "scripts"))

from standard_records import (  # noqa: E402
    fact,
    new_molecule,
    new_outcome,
    new_participant,
    new_reaction,
    new_step,
)


RUN = ROOT / "runs" / "current"
DOC = "doc_832eb8cea3728513"


def ev(page: int, image: bool = True) -> list[str]:
    ids = [f"{DOC}:page{page}:text"]
    if image and page <= 6:
        ids.append(f"{DOC}:page{page}:image")
    return ids


def missing(evidence, status="unreported"):
    return fact(evidence=evidence, status=status, verification="missing")


def extracted(raw, normalized=None, evidence=()):
    return fact(raw=raw, normalized=raw if normalized is None else normalized,
                evidence=evidence, method="original_extraction")


def inferred(raw, normalized, evidence, derivation):
    return fact(raw=raw, normalized=normalized, evidence=evidence,
                method="inference", verification="inferred", derivation=derivation)


def not_applicable(evidence):
    return fact(evidence=evidence, method="not_applicable", status="not_applicable",
                verification="missing")


MOLECULES: dict[str, dict] = {}


def add_molecule(key, label, page, *, stereo="unresolved", geometry="unresolved",
                 representation="not_encoded", assignment=None, mixture=None,
                 identity_note=None):
    evidence = ev(page)
    mid = f"{DOC}:mol:{key}"
    m = new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = extracted(label, label, evidence)
    m["names"]["original"] = extracted(label, label, evidence)
    m["names"]["english"] = extracted(label, label, evidence)
    m["names"]["abbreviations"] = missing(evidence)
    m["names"]["aliases"] = missing(evidence)
    m["cas_numbers"] = missing(evidence)
    m["representation_type"] = "unresolved"
    m["smiles"] = missing(evidence)
    m["generic_template"] = None
    st = m["stereochemistry"]
    st["atom_centers"] = missing(evidence)
    st["double_bonds"] = missing(evidence)
    st["axial"] = missing(evidence)
    st["mixture"] = extracted(mixture, mixture, evidence) if mixture else missing(evidence)
    st["assignment"] = extracted(assignment, assignment, evidence) if assignment else missing(evidence)
    status_basis = "Status assigned from the inspected source drawing/label; no candidate SMILES was asserted."
    st["chirality_status"] = inferred(label, stereo, evidence, status_basis)
    st["geometry_status"] = inferred(label, geometry, evidence, status_basis)
    st["representation_status"] = inferred(label, representation, evidence, status_basis)
    st["isomeric_smiles"] = missing(evidence)
    st["identity_key"] = missing(evidence)
    note = identity_note or (
        "The source label/drawing was registered independently. Candidate SMILES is left empty "
        "because complete atom-by-atom transcription and stereochemical encoding were not established."
    )
    m["identity_note"] = inferred(label, note, evidence,
                                  "Conservative structure-gap decision after visual inspection.")
    m["external_properties"] = []
    MOLECULES[key] = m
    return mid


def add_agent(key, label, page):
    return add_molecule(key, label, page, stereo="unresolved", geometry="unresolved",
                        representation="not_encoded",
                        identity_note="Named participant is explicit in the source; structure and CAS were not imported or guessed.")


def mid(key):
    return MOLECULES[key]["molecule_id"]


# Source compounds.  Relative cis/trans and dr information is retained even
# though no machine SMILES is promoted to a verified structure.
for key, label, page, stereo, geometry, assignment, mixture in [
    ("11", "compound 11 (alkynoate)", 3, "achiral", "not_applicable", None, None),
    ("9", "compound 9 (allylboronate, Z/E mixture)", 2, "achiral", "mixture", "Z/E geometry reported", "Z/E mixture"),
    ("12", "compound 12, 3-phenyl-2-propynal", 3, "achiral", "not_applicable", None, None),
    ("13", "compound 13, 3-triisopropylsilyl-2-propynal", 3, "achiral", "not_applicable", None, None),
    ("13_tms", "3-trimethylsilyl-2-propynal (failed substrate)", 2, "achiral", "not_applicable", None, None),
    ("S4", "compound S4, methyl (Z,E)-5-(2-methyl-1,3-dioxolan-2-yl)pent-2-enoate", 3, "achiral", "mixture", "Z/E mixture reported", "Z/E mixture"),
    ("14a", "compound 14a (trans lactam)", 3, "partially_specified", "specified", "trans relative lactam configuration; absolute configuration not assigned", None),
    ("14b", "compound 14b (cis lactam)", 3, "partially_specified", "specified", "cis relative lactam configuration; absolute configuration not assigned", None),
    ("15a", "compound 15a (trans lactam)", 3, "partially_specified", "specified", "trans relative lactam configuration; absolute configuration not assigned", None),
    ("15b", "compound 15b (cis lactam)", 3, "partially_specified", "specified", "cis relative lactam configuration; absolute configuration not assigned", None),
    ("18a", "compound 18a (trans keto lactam)", 3, "partially_specified", "specified", "trans relative lactam configuration", None),
    ("18b", "compound 18b (cis keto lactam)", 3, "partially_specified", "specified", "cis relative lactam configuration", None),
    ("19a", "compound 19a (trans TIPS keto lactam)", 3, "partially_specified", "specified", "trans relative lactam configuration", None),
    ("19b", "compound 19b (cis TIPS keto lactam)", 3, "partially_specified", "specified", "cis relative lactam configuration", None),
    ("20a", "compound 20a (trans propargyl alcohol)", 3, "mixture", "specified", "trans lactam relationship; new alcohol-bearing center unassigned", "approximately 1:1 diastereomeric mixture at the newly created stereocenter"),
    ("20b", "compound 20b (cis propargyl alcohol)", 3, "mixture", "specified", "cis lactam relationship; new alcohol-bearing center unassigned", "approximately 1:1 diastereomeric mixture at the newly created stereocenter"),
    ("21a", "compound 21a (trans TIPS propargyl alcohol)", 3, "mixture", "specified", "trans lactam relationship; new alcohol-bearing center unassigned", "approximately 1:1 diastereomeric mixture at the newly created stereocenter"),
    ("21b", "compound 21b (cis TIPS propargyl alcohol)", 3, "mixture", "specified", "cis lactam relationship; new alcohol-bearing center unassigned", "approximately 1:1 diastereomeric mixture at the newly created stereocenter"),
    ("22a", "compound 22a (trans propargyl pivalate)", 3, "mixture", "specified", "trans lactam relationship", "1:1 diastereomeric mixture at C10"),
    ("22b", "compound 22b (cis propargyl pivalate)", 3, "mixture", "specified", "cis lactam relationship", "1:1 diastereomeric mixture at C10"),
    ("23a", "compound 23a (trans TIPS propargyl pivalate)", 3, "mixture", "specified", "trans lactam relationship", "1:1 diastereomeric mixture at C10"),
    ("23b", "compound 23b (cis TIPS propargyl pivalate)", 3, "mixture", "specified", "cis lactam relationship", "1:1 diastereomeric mixture at C10"),
    ("24a", "compound 24a (trans phenyl allene-yne)", 3, "partially_specified", "specified", "trans relative lactam configuration; terminal CH2 allene is not axially stereogenic", None),
    ("24b", "compound 24b (cis phenyl allene-yne)", 5, "partially_specified", "specified", "cis relative lactam configuration; terminal CH2 allene is not axially stereogenic", None),
    ("25a", "compound 25a (trans TIPS allene-yne)", 3, "partially_specified", "specified", "trans relative lactam configuration; terminal CH2 allene is not axially stereogenic", None),
    ("25b", "compound 25b (cis TIPS allene-yne)", 5, "partially_specified", "specified", "cis relative lactam configuration; terminal CH2 allene is not axially stereogenic", None),
    ("26", "compound 26 (terminal-alkyne allene-yne)", 3, "partially_specified", "specified", "trans relative lactam configuration", None),
    ("S8", "compound S8 (N-methyl ketal intermediate)", 4, "partially_specified", "specified", "predominantly trans relative lactam configuration", None),
    ("27a", "compound 27a (N-methyl keto lactam)", 3, "partially_specified", "specified", "trans relative lactam configuration", None),
    ("28a", "compound 28a (N-methyl propargyl alcohol)", 3, "mixture", "specified", "7:1 trans/cis lactam ratio; new center unassigned", "1:1 diastereomeric mixture at newly formed center"),
    ("29a", "compound 29a (N-methyl propargyl acetate)", 3, "mixture", "specified", "predominantly trans lactam relationship", "1:1 diastereomeric mixture"),
    ("30a", "compound 30a (N-methyl allene-yne)", 3, "partially_specified", "specified", "trans relative lactam configuration", None),
    ("31a", "compound 31a (phenyl reduced methyl-lactam byproduct)", 4, "partially_specified", "specified", "relative configuration depicted; structure reported as not fully confirmed", None),
    ("32a", "compound 32a (TIPS reduced methyl-lactam byproduct)", 4, "mixture", "specified", "relative configuration depicted", "9:1 dr reported in Eq 1 drawing"),
    ("cis_alkene", "unlabelled cis-alkene reduction product from 22a", 4, "unresolved", "specified", "cis alkene stated; full product identity not assigned in Table 1", None),
    ("33", "compound 33, iodobenzene", 5, "achiral", "not_applicable", None, None),
    ("34", "compound 34, 4-iodobenzonitrile", 5, "achiral", "not_applicable", None, None),
    ("35", "compound 35, 4-iodobenzotrifluoride", 5, "achiral", "not_applicable", None, None),
    ("36", "compound 36, 4-iodoanisole", 5, "achiral", "not_applicable", None, None),
    ("37", "compound 37 (N-phenyl allene-yne)", 5, "partially_specified", "specified", "relative ring stereochemistry depicted", None),
    ("38", "compound 38 (N-4-cyanophenyl allene-yne)", 5, "partially_specified", "specified", "relative ring stereochemistry depicted", None),
    ("39", "compound 39 (N-4-trifluoromethylphenyl allene-yne)", 5, "partially_specified", "specified", "relative ring stereochemistry depicted", None),
    ("40", "compound 40 (N-4-methoxyphenyl allene-yne)", 5, "partially_specified", "specified", "relative ring stereochemistry depicted", None),
    ("41", "compound 41, 2-iodothiophene", 5, "achiral", "not_applicable", None, None),
    ("42", "compound 42 (N-2-thiophenyl allene-yne)", 5, "partially_specified", "specified", "relative ring stereochemistry depicted", None),
    ("43", "compound 43 (N-tosyl allene-yne)", 5, "partially_specified", "specified", "relative ring stereochemistry depicted", None),
    ("44", "compound 44 (N-Boc allene-yne)", 5, "partially_specified", "specified", "relative ring stereochemistry depicted", None),
    ("45", "compound 45 (APKR product)", 6, "partially_specified", "specified", "relative fused-ring stereochemistry depicted and X-ray structure shown; absolute configuration not assigned", None),
    ("46", "compound 46 (TIPS APKR product)", 6, "partially_specified", "specified", "relative fused-ring stereochemistry depicted", None),
    ("47", "compound 47 (terminal-alkyne APKR product)", 6, "partially_specified", "specified", "relative fused-ring stereochemistry depicted", None),
    ("48", "compound 48 (N-methyl APKR product)", 6, "partially_specified", "specified", "relative fused-ring stereochemistry depicted", None),
    ("49", "compound 49 (N-phenyl APKR product)", 6, "partially_specified", "specified", "relative fused-ring stereochemistry depicted", None),
    ("50", "compound 50 (N-4-cyanophenyl APKR product)", 6, "partially_specified", "specified", "relative fused-ring stereochemistry depicted", None),
    ("51", "compound 51 (N-4-trifluoromethylphenyl APKR product)", 6, "partially_specified", "specified", "relative fused-ring stereochemistry depicted", None),
    ("52", "compound 52 (N-4-methoxyphenyl APKR product)", 6, "partially_specified", "specified", "relative fused-ring stereochemistry depicted", None),
    ("53", "compound 53 (N-2-thiophenyl APKR product)", 6, "partially_specified", "specified", "relative fused-ring stereochemistry depicted", None),
    ("54", "compound 54 (N-tosyl APKR product)", 6, "partially_specified", "specified", "relative fused-ring stereochemistry depicted", None),
    ("55", "compound 55 (N-Boc APKR product)", 6, "partially_specified", "specified", "relative fused-ring stereochemistry depicted", None),
    ("56", "compound 56 (trans,trans alpha-methyl lactam allene-yne)", 6, "mixture", "specified", "trans,trans relative assignment; C11 not absolutely assigned", "8:1 dr at C11"),
    ("57", "compound 57 (cis,trans alpha-methyl lactam allene-yne)", 6, "partially_specified", "specified", "cis,trans relative assignment", None),
    ("58", "compound 58 (alpha-methyl APKR product)", 6, "mixture", "specified", "relative assignment; C11 stereochemistry not fully assigned", "7:1 dr in Scheme 4 (text also reports 8:1)"),
    ("59", "compound 59 (phenyl alpha-methyl APKR product)", 6, "partially_specified", "specified", "single diastereomer; C11 assignment based on calculated/experimental coupling comparison", None),
    ("60", "compound 60 (N-acetyl APKR product)", 6, "partially_specified", "specified", "relative fused-ring stereochemistry follows compound 47; exact encoding not established", None),
    ("S17", "compound S17 (alpha-methylene-gamma-lactone analogue)", 6, "partially_specified", "specified", "relative fused-ring stereochemistry depicted in source context outside the covered image range", None),
]:
    add_molecule(key, label, page, stereo=stereo, geometry=geometry,
                 assignment=assignment, mixture=mixture)


# Exact, fully legible small structures were transcribed and checked with
# `reaction_tools.py standard stereo`.  The checker explicitly does not confer
# source verification; these remain source-identified machine records.
STEREO_CHECKS = {
    "12": ("O=CC#Cc1ccccc1", "1a0e4b1cfea89862371f042d4764f832715e53d2412f8b2778eedc76959cc9b4",
           [], [{"bond_index": 0, "stereo": "STEREONONE"}]),
    "33": ("Ic1ccccc1", "87e1da21bf5ad715fb197e6ecfec7d903bfa25b0e37e8604ea6be7b787b05cf9", [], []),
    "34": ("N#Cc1ccc(I)cc1", "f1c8ca11f1b6daca64b22e409032e16d8f9e5b034ad0c519f8db476537c4da8d", [], []),
    "35": ("FC(F)(F)c1ccc(I)cc1", "260d950ba165aa61116957f1dfba1f1400fa2ffd50ee58013d4a02c9e18a234e", [], []),
    "36": ("COc1ccc(I)cc1", "669aacc5d1162b2ae1164677e465c7d8a71a2f2e2efaf3fde1054a7bf0f0eb3d", [], []),
    "41": ("Ic1cccs1", "84f031d7e220a5ac7c6dea03edd10f660af93c436dbf70d5910d3decfa770ca8", [], []),
}
for key, (smiles, graph_key, centers, double_bonds) in STEREO_CHECKS.items():
    m = MOLECULES[key]
    evidence = m["occurrences"][0]["evidence_ids"]
    derivation = "SMILES transcribed from the fully legible source structure and graph-checked with standard stereo; source_verified=false."
    m["representation_type"] = "smiles"
    m["smiles"] = fact(raw=m["names"]["original"]["raw_value"], normalized=smiles,
                       evidence=evidence, method="structure_derivation", verification="inferred", derivation=derivation)
    st = m["stereochemistry"]
    st["atom_centers"] = inferred("standard stereo atom_centers", centers, evidence, derivation)
    st["double_bonds"] = inferred("standard stereo double_bonds", double_bonds, evidence, derivation)
    st["representation_status"] = inferred("standard stereo", "not_applicable", evidence,
                                            "No stereogenic representation is applicable to this graph.")
    st["isomeric_smiles"] = fact(raw=smiles, normalized=smiles, evidence=evidence,
                                  method="structure_derivation", verification="inferred", derivation=derivation)
    st["identity_key"] = inferred("standard stereo graph key", graph_key, evidence, derivation)
    m["identity_note"] = inferred(m["names"]["original"]["raw_value"],
                                   "Complete small-molecule graph transcribed; no stereogenic element detected by the graph checker.",
                                   evidence, derivation)


# Explicit allene-axis review.  These displayed allenes terminate in =CH2, so
# the two H substituents at that end are equivalent and the allene axis is not
# stereogenic.  This statement does not resolve the separate ring centers.
for key in ("24a", "24b", "25a", "25b", "26", "30a", "31a", "32a",
            "37", "38", "39", "40", "42", "43", "44", "56", "57"):
    m = MOLECULES[key]
    evidence = m["occurrences"][0]["evidence_ids"]
    m["stereochemistry"]["axial"] = inferred(
        "terminal =CH2 allene shown in source",
        {"status": "not_stereogenic", "basis": "one allene terminus bears two equivalent hydrogen substituents"},
        evidence,
        "Axial-stereogenicity review of the depicted cumulene substituents; ring stereocenters are assessed separately."
    )


for key, label, page in [
    ("CuI", "copper iodide (CuI)", 3), ("MeLi", "methyl lithium (MeLi)", 3),
    ("Et2O", "diethyl ether", 3), ("THF", "tetrahydrofuran (THF)", 3),
    ("toluene", "toluene", 3), ("HMPA", "hexamethylphosphoramide (HMPA)", 3),
    ("DIBALH", "diisobutylaluminum hydride (DIBALH)", 3), ("hexanes", "hexanes", 3),
    ("PinBCH2Cl", "2-(chloromethyl)-4,4,5,5-tetramethyl-1,3-dioxaborolane", 3),
    ("NH4OH", "ammonium hydroxide", 3), ("ethanol", "ethanol", 3),
    ("PPTS", "pyridinium para-toluenesulfonate (PPTS)", 3), ("acetone", "acetone", 3),
    ("water", "water", 3), ("ethynylMgBr", "ethynyl magnesium bromide", 3),
    ("ScOTf3", "scandium triflate, Sc(OTf)3", 3), ("piv2O", "pivalic anhydride", 3),
    ("MeCN", "acetonitrile", 3), ("Stryker", "triphenylphosphine copper hydride hexamer (Stryker's reagent)", 3),
    ("TBAF", "tetra-n-butylammonium fluoride (TBAF)", 3), ("NaH", "sodium hydride", 3),
    ("MeI", "iodomethane", 3), ("DMF", "dimethylformamide (DMF)", 3),
    ("Ac2O", "acetic anhydride", 3), ("Et3N", "triethylamine", 3),
    ("DMAP", "4-(N,N-dimethylamino)pyridine (DMAP)", 3), ("DCM", "dichloromethane", 3),
    ("DMEDA", "N,N'-dimethylethylenediamine", 5), ("Cs2CO3", "cesium carbonate", 5),
    ("pTsCl", "para-toluenesulfonyl chloride", 5), ("Boc2O", "di-tert-butyl dicarbonate", 5),
    ("Rhcat", "rhodium biscarbonyl chloride dimer, [Rh(CO)2Cl]2", 5), ("CO", "carbon monoxide", 5),
    ("cysteamine", "cysteamine (HSCH2CH2NH2)", 6), ("CDCl3", "deuterated chloroform (CDCl3)", 6),
]:
    add_agent(key, label, page)


REACTIONS = []


def q(number, unit, **extra):
    out = {"number": number, "unit": unit}
    out.update(extra)
    return out


def participant(key, role, evidence, *, equiv=None, molpct=None, stage="reaction"):
    pid = f"p{len(_PARTICIPANTS) + 1}"
    p = new_participant(pid, mid(key), role, evidence, stage=stage)
    p["product_role"] = extracted("product", "product", evidence) if role == "product" else not_applicable(evidence)
    p["charges"]["mass"] = missing(evidence)
    p["charges"]["amount"] = missing(evidence)
    p["charges"]["volume"] = missing(evidence)
    p["charges"]["equivalents"] = extracted(f"{equiv} equiv", q(equiv, "equiv"), evidence) if equiv is not None else missing(evidence)
    p["charges"]["mol_percent"] = extracted(f"{molpct} mol %", q(molpct, "mol%"), evidence) if molpct is not None else missing(evidence)
    if equiv is not None or molpct is not None:
        # Quantities in the displayed schemes/tables are stated relative to the
        # substrate.  Store the registry ID, as required by the v3.1 contract.
        reference_mid = _PARTICIPANTS[0]["molecule_id"] if _PARTICIPANTS else mid(key)
        p["equivalence_reference"] = extracted("relative to substrate", reference_mid, evidence)
    else:
        p["equivalence_reference"] = missing(evidence)
    _PARTICIPANTS.append(p)
    return pid


def step(step_id, evidence, operation, temperature=None, time=None, *, participants=(), atmosphere=None, pressure=None):
    s = new_step(step_id, evidence)
    s["participant_ids"] = list(participants)
    s["operation"] = extracted(operation, operation, evidence)
    if isinstance(temperature, str) and any(ch.isdigit() for ch in temperature):
        # A mixed numeric/qualitative ramp such as "0 °C to rt" cannot be
        # normalized to a numeric range without inventing the rt endpoint.
        s["temperature"] = fact(raw=temperature, normalized=None, evidence=evidence)
    else:
        s["temperature"] = extracted(str(temperature), temperature, evidence) if temperature is not None else missing(evidence)
    s["time"] = extracted(str(time), time, evidence) if time is not None else missing(evidence)
    s["atmosphere"] = extracted(atmosphere, atmosphere, evidence) if atmosphere is not None else missing(evidence)
    s["pressure"] = extracted(str(pressure), pressure, evidence) if pressure is not None else missing(evidence)
    return s


LABELS = {
    "C-H": 1, "C-C": 2, "C-O": 3, "C-S": 4, "C-N": 5, "C-B": 11,
    "DB-C-C": 13, "DB-C-O": 14, "N-S": 31,
}


def classification(r, labels, evidence, rationale):
    c = r["classification"]
    c["mapped_reaction_smiles"] = missing(evidence)
    atom_pairs = {
        "C-H": ("C", "H", 1), "C-C": ("C", "C", 1), "C-O": ("C", "O", 1),
        "C-S": ("C", "S", 1), "C-N": ("C", "N", 1), "C-B": ("C", "B", 1),
        "DB-C-C": ("C", "C", 2), "DB-C-O": ("C", "O", 2), "N-S": ("N", "S", 1),
    }
    change_values = []
    label_values = []
    for index, code in enumerate(labels, 1):
        a, b, order = atom_pairs[code]
        change_id = f"b{index}"
        change = {
            "change_id": change_id,
            "atoms": [{"atom_map": None, "element": a}, {"atom_map": None, "element": b}],
            "before_order": 0,
            "after_order": order,
            "change_type": "formed",
        }
        change_values.append(inferred(
            rationale, change, evidence,
            "Net formed-bond type was inferred from the inspected source drawing; atom identities are intentionally unmapped."
        ))
        label_values.append({"category_number": LABELS[code], "code": code, "supporting_change_ids": [change_id]})
    c["bond_changes"] = change_values
    c["labels"] = inferred(rationale, label_values, evidence,
                           "Formation labels inferred from source reactant/product drawings under bond-formation-20261016.")
    c["coverage"] = inferred(rationale, "partial", evidence,
                             "No atom-mapped reaction SMILES; only visually supported net formed-bond classes are asserted.")
    c["limitations"] = ["No atom-mapped reaction SMILES.", "Hydrogen source and complete bond-change inventory remain unresolved."]
    c["model_votes"] = []


def outcome(outcome_id, kind, targets, value, evidence, *, role="product", yield_type=None,
            basis=None, scope=None, status=None):
    o = new_outcome(outcome_id, kind, [mid(t) for t in targets], evidence)
    o["product_role"] = extracted(role, role, evidence) if role is not None else not_applicable(evidence)
    if status:
        o["value"] = fact(evidence=evidence, status=status, verification="missing")
    else:
        raw = value if isinstance(value, str) else json.dumps(value, ensure_ascii=False)
        o["value"] = extracted(raw, value, evidence)
    o["yield_type"] = extracted(yield_type, yield_type, evidence) if yield_type else not_applicable(evidence)
    o["measurement_basis"] = extracted(basis, basis, evidence) if basis else missing(evidence)
    o["scope"] = extracted(scope, scope, evidence) if scope else missing(evidence)
    return o


def conflict_fact(alternatives, evidence):
    return fact(raw=None, normalized=None, evidence=evidence, status="conflict", verification="unreviewed",
                alternatives=[extracted(raw, normalized, evs) for raw, normalized, evs in alternatives])


def reaction(rid, page, location, reactants, products, agents, solvents, steps, outcomes, labels,
             *, category="substrate_scope", name="", transformation="", scale=None, notes=(),
             extra_evidence=()):
    global _PARTICIPANTS
    evidence = list(dict.fromkeys(ev(page) + list(extra_evidence)))
    rid = f"{DOC}:rxn:{rid}"
    r = new_reaction(rid, DOC, evidence)
    r["record_type"] = "specific_experiment"
    r["experiment_categories"] = extracted(category, [category], evidence)
    r["location"] = extracted(location, {"page": page, "region": location}, evidence)
    r["general_scheme_ids"] = []
    r["reaction_name"] = extracted(name, name, evidence) if name else missing(evidence)
    r["transformation_type"] = extracted(transformation, transformation, evidence) if transformation else missing(evidence)
    _PARTICIPANTS = []
    for key in reactants:
        participant(key, "reactant", evidence)
    for key in products:
        participant(key, "product", evidence)
    for spec in agents:
        if isinstance(spec, str):
            participant(spec, "reagent", evidence)
        else:
            key, role, eqv, molpct = spec
            participant(key, role, evidence, equiv=eqv, molpct=molpct)
    for key in solvents:
        participant(key, "solvent", evidence)
    r["participants"] = _PARTICIPANTS
    if scale:
        scale = dict(scale)
        scale["reference_molecule_id"] = mid(reactants[0])
        r["reaction_scale"] = extracted(str(scale), scale, evidence)
    else:
        r["reaction_scale"] = missing(evidence)
    r["solvent_conditions"]["composition"] = extracted(
        ", ".join(MOLECULES[k]["names"]["original"]["normalized_value"] for k in solvents), solvents, evidence
    ) if solvents else missing(evidence)
    r["solvent_conditions"]["mixture_ratio"] = missing(evidence)
    r["solvent_conditions"]["total_volume"] = missing(evidence)
    r["solvent_conditions"]["concentration"] = missing(evidence)
    r["solvent_conditions"]["concentration_basis"] = missing(evidence)
    r["steps"] = steps(evidence, r) if callable(steps) else steps
    for k in r["special_conditions"]:
        r["special_conditions"][k] = not_applicable(evidence)
    r["procedure"]["original_text"] = extracted(location, location, evidence)
    r["procedure"]["general_procedure_reference"] = missing(evidence)
    r["procedure"]["modifications"] = missing(evidence)
    r["procedure"]["addition_order"] = missing(evidence)
    r["procedure"]["pretreatment"] = missing(evidence)
    r["outcomes"] = outcomes(evidence, r) if callable(outcomes) else outcomes
    classification(r, labels, evidence, transformation or name or location)
    r["notes"] = [extracted(note, note, evidence) for note in notes]
    r["review"] = {
        "status": "needs_review", "method": "model", "reviewer": None,
        "checked_fields": [],
        "unresolved": [
            "Machine record remains needs_review; schema validation is not chemistry verification.",
            "Candidate SMILES/atom mapping and complete stereochemical encoding were not established from pages 1-6.",
            "Unreported quantities, workup and purification details were not supplied from an unassociated SI."
        ],
    }
    r["completeness"] = "partial"
    REACTIONS.append(r)
    return r


def simple_steps(operation, temp, tm, *, atmosphere=None, pressure=None):
    return lambda e, r: [step("s1", e, operation, temp, tm, atmosphere=atmosphere, pressure=pressure)]


def yield_out(targets, pct, *, ytype="reported yield", basis=None, scope=None):
    return lambda e, r: [outcome("o1", "yield", targets, q(pct, "%"), e, yield_type=ytype, basis=basis, scope=scope)]


# Page 2: explicit failed substrate test.  The shared condition is inherited
# from Scheme 2(b) on page 3 rather than invented from convention.
r = reaction(
    "rxn_p2_tms_alkynal_failure", 2, "right-column final paragraph; failed TMS alkynal lactamization",
    ["9", "13_tms"], [], ["NH4OH"], ["ethanol"],
    lambda e, r: [step("s1", e, "subjected to the lactamization reaction conditions (Scheme 2b)", "r.t.", q(16, "h"))],
    lambda e, r: [outcome("o1", "decomposition", ["13_tms"], "unstable; presumed loss of trimethylsilyl group", e,
                               role=None, scope="failed substrate observation")],
    [], category="other", name="allylboration/lactamization substrate test",
    transformation="failed three-component lactamization",
    notes=("No product was shown; no product participant was invented.",), extra_evidence=ev(3)
)


# Scheme 2A reaction 11 -> 9.  The figure label, footnote and prose are kept as
# explicit alternatives because their yields and Z/E ratios differ.
def steps_11_9(e, r):
    return [
        step("s1", e, "combine CuI and MeLi in Et2O/THF", q(-30, "°C"), q(30, "min")),
        step("s2", e, "add toluene, HMPA and DIBALH in hexanes; hold", q(-30, "°C"), q(2, "h")),
        step("s3", e, "add compound 11; hold", q(-20, "°C"), q(5, "h")),
        step("s4", e, "add PinBCH2Cl at -20 °C", q(-20, "°C"), None),
        step("s5", e, "warm from -20 °C to room temperature; ramp time unreported, then hold", "-20 °C to r.t.", q(16, "h")),
    ]


def outs_11_9(e, r):
    oy = outcome("o1", "yield", ["9"], q(81, "%"), e, yield_type="average yield", scope="compound 9")
    oy["value"] = conflict_fact([
        ("81% (product label/prose average)", q(81, "%"), ev(2) + ev(3)),
        ("76% (Scheme 2 footnote a)", q(76, "%"), ev(3)),
    ], e)
    os = outcome("o2", "selectivity", ["9"], "2:1", e, role="product mixture", scope="Z/E ratio, Z:E")
    os["value"] = conflict_fact([
        ("2:1 Z/E (product label/prose)", "2:1", ev(2) + ev(3)),
        ("3:1 Z/E (Scheme 2 footnote a)", "3:1", ev(3)),
    ], e)
    ob = outcome("o3", "yield", ["S4"], q(10, "%"), e, role="byproduct", yield_type="reported byproduct amount",
                 scope="S4 Z/E mixture")
    ob["value"] = conflict_fact([
        ("10% (Scheme 2 footnote a)", q(10, "%"), ev(3)),
        ("5-10% (page 2 prose, repeated runs)", {"minimum": 5, "maximum": 10, "unit": "%", "raw_text": "5-10%"}, ev(2)),
    ], e)
    return [oy, os, ob]


r = reaction(
    "rxn_p3_s2a_11_to_9", 3, "Scheme 2A arrow (a), compound 11 to compound 9",
    ["11"], ["9", "S4"], ["CuI", "MeLi", "HMPA", "DIBALH", "PinBCH2Cl"],
    ["Et2O", "THF", "toluene", "hexanes"], steps_11_9, outs_11_9,
    ["C-C", "C-B", "C-H"], name="allylboronate formation", transformation="alkynoate to allylboronate",
    notes=("Figure/product label and Scheme 2 footnote disagree on yield and Z/E ratio; neither value was selected.",)
)
r["review"]["unresolved"].append("Resolve whether 81%/2:1 and 76%/3:1 are distinct runs or a source-level inconsistency.")


def add_lactamization(rid, aldehyde, prods, yld, ratio, page=3):
    return reaction(
        rid, page, f"Scheme 2A arrow (b), compound 9 + {aldehyde} to {','.join(prods)}",
        ["9", aldehyde], prods, ["NH4OH"], ["ethanol"], simple_steps("allylboration/lactamization", "r.t.", q(16, "h")),
        lambda e, r: [
            outcome("o1", "yield", prods, q(yld, "%"), e, role="combined product mixture", yield_type="reported yield"),
            outcome("o2", "selectivity", prods, ratio, e, role="product mixture", scope=f"trans/cis ratio in order {prods[0]}/{prods[1]}; range retained"),
        ], ["C-C", "C-N"], name="three-component allylboration/lactamization",
        transformation="2-alkynal imine allylboration followed by lactamization"
    )


add_lactamization("rxn_p3_s2a_9_12_to_14ab", "12", ["14a", "14b"], 58, "5-2:1")
add_lactamization("rxn_p3_s2a_9_13_to_15ab", "13", ["15a", "15b"], 71, "3.5-2:1")


def add_deprotection(rid, inputs, outputs, yld):
    return reaction(
        rid, 3, f"Scheme 2A arrow (c), {','.join(inputs)} to {','.join(outputs)}",
        inputs, outputs, ["PPTS"], ["acetone", "water"],
        simple_steps("acid-catalyzed ketal hydrolysis in acetone/water (15:1)", "reflux", q(16, "h")),
        lambda e, r: [
            outcome("o1", "yield", outputs, q(yld, "%"), e, role="combined product mixture", yield_type="reported yield"),
            outcome("o2", "selectivity", outputs, "4:1", e, role="product mixture", scope=f"trans/cis ratio in order {outputs[0]}/{outputs[1]}")
        ], ["DB-C-O"], name="ketal deprotection", transformation="ketal hydrolysis to methyl ketone",
        extra_evidence=ev(4)
    )


add_deprotection("rxn_p3_s2a_14ab_to_18ab", ["14a", "14b"], ["18a", "18b"], 83)
add_deprotection("rxn_p3_s2a_15ab_to_19ab", ["15a", "15b"], ["19a", "19b"], 73)


def add_grignard(rid, inputs, outputs, yld):
    return reaction(
        rid, 3, f"Scheme 2A arrow (d), {','.join(inputs)} to {','.join(outputs)}",
        inputs, outputs, ["ethynylMgBr"], ["THF"],
        simple_steps("addition of ethynyl magnesium bromide", q(0, "°C"), q(3, "h")),
        lambda e, r: [
            outcome("o1", "yield", outputs, q(yld, "%"), e, role="combined product mixture", yield_type="reported yield"),
            outcome("o2", "selectivity", outputs, "approximately 1:1", e, role="product mixture", scope="diastereomeric ratio at newly created center; no trans/cis ratio assigned to this step")
        ], ["C-C"], name="ethynylmagnesium bromide addition", transformation="methyl ketone to tertiary propargyl alcohol",
        extra_evidence=ev(4)
    )


add_grignard("rxn_p3_s2a_18ab_to_20ab", ["18a", "18b"], ["20a", "20b"], 75)
add_grignard("rxn_p3_s2a_19ab_to_21ab", ["19a", "19b"], ["21a", "21b"], 82)


def add_pivalation(rid, inputs, outputs, yields):
    return reaction(
        rid, 3, f"Scheme 2A arrow (e), {','.join(inputs)} to {','.join(outputs)}",
        inputs, outputs, ["piv2O", ("ScOTf3", "reagent", 0.4, None)], ["MeCN"],
        simple_steps("chemoselective pivaloylation", "r.t.", q(16, "h")),
        lambda e, r: [
            *[outcome(f"o{i+1}", "yield", [out], q(y, "%"), e, role="isolated stereoisomer", yield_type="reported yield") for i, (out, y) in enumerate(zip(outputs, yields))],
            *[outcome(f"o{len(outputs)+i+1}", "selectivity", [out], "1:1", e, role="product mixture", scope="diastereomeric ratio at C10") for i, out in enumerate(outputs)]
        ], ["C-O"], name="pivaloylation", transformation="propargyl alcohol to propargyl pivalate",
        extra_evidence=ev(4)
    )


add_pivalation("rxn_p3_s2a_20ab_to_22ab", ["20a", "20b"], ["22a", "22b"], [61, 16])
add_pivalation("rxn_p3_s2a_21ab_to_23ab", ["21a", "21b"], ["23a", "23b"], [39, 16])


def add_allene(rid, inp, out, yld):
    return reaction(
        rid, 3, f"Scheme 2A arrow (f), {inp} to {out}", [inp], [out], ["Stryker"], ["toluene"],
        simple_steps("Stryker-reagent conversion of propargyl pivalate to 3,3-disubstituted allene", q(-10, "°C"), q(2, "h")),
        yield_out([out], yld), ["DB-C-C", "C-H"], name="allene formation with Stryker's reagent",
        transformation="propargyl pivalate to allene", extra_evidence=ev(4)
    )


add_allene("rxn_p3_s2a_22a_to_24a", "22a", "24a", 33)
add_allene("rxn_p3_s2a_23a_to_25a", "23a", "25a", 46)
reaction(
    "rxn_p3_s2a_25a_to_26", 3, "Scheme 2A arrow (g), compound 25a to 26", ["25a"], ["26"],
    ["TBAF"], ["THF"], simple_steps("TIPS deprotection with TBAF", q(0, "°C"), q(45, "min")),
    yield_out(["26"], 79), ["C-H"], name="TIPS deprotection", transformation="silyl alkyne to terminal alkyne", extra_evidence=ev(4)
)


def zero_to_rt_steps(operation, total):
    return lambda e, r: [
        step("s1", e, operation + "; additions at 0 °C", q(0, "°C"), None),
        step("s2", e, "warm to room temperature; ramp time not reported", "0 °C to r.t.", total),
    ]


reaction(
    "rxn_p3_s2b_14a_to_S8", 3, "Scheme 2B arrow (h), N-methylation of 14a", ["14a"], ["S8"],
    ["NaH", "MeI"], ["DMF"], zero_to_rt_steps("add sodium hydride then iodomethane", q(15, "min")),
    yield_out(["S8"], 71), ["C-N"], name="N-methylation", transformation="secondary lactam to N-methyl lactam",
    notes=("Page 4 prose states the 14a material contained a small amount of cis 14b.",), extra_evidence=ev(4)
)
reaction(
    "rxn_p3_s2b_S8_to_27a", 3, "Scheme 2B arrow (c), ketal hydrolysis to 27a", ["S8"], ["27a"],
    ["PPTS"], ["acetone", "water"], simple_steps("ketal hydrolysis in acetone/water (15:1)", "reflux", q(16, "h")),
    yield_out(["27a"], 73), ["DB-C-O"], name="ketal deprotection", transformation="ketal to methyl ketone", extra_evidence=ev(4)
)
reaction(
    "rxn_p3_s2b_27a_to_28a", 3, "Scheme 2B arrow (d), 27a to 28a", ["27a"], ["28a"],
    ["ethynylMgBr"], ["THF"], simple_steps("ethynyl magnesium bromide addition", q(0, "°C"), q(3, "h")),
    lambda e, r: [
        outcome("o1", "yield", ["28a"], q(62, "%"), e, yield_type="reported yield"),
        outcome("o2", "selectivity", ["28a"], "7:1", e, scope="trans/cis ratio"),
        outcome("o3", "selectivity", ["28a"], "1:1", e, scope="diastereomeric ratio at newly formed center")
    ], ["C-C"], name="ethynylmagnesium bromide addition", transformation="methyl ketone to propargyl alcohol", extra_evidence=ev(4)
)
reaction(
    "rxn_p3_s2b_28a_to_29a", 3, "Scheme 2B arrow (i), 28a to 29a", ["28a"], ["29a"],
    ["Ac2O", "Et3N", "DMAP"], ["DCM"], zero_to_rt_steps("acetylation", q(3, "h")),
    lambda e, r: [outcome("o1", "yield", ["29a"], q(54, "%"), e, yield_type="reported yield"),
                  outcome("o2", "selectivity", ["29a"], "1:1", e, scope="diastereomeric ratio")],
    ["C-O"], name="acetylation", transformation="propargyl alcohol to propargyl acetate", extra_evidence=ev(4)
)
reaction(
    "rxn_p3_s2b_29a_to_30a", 3, "Scheme 2B arrow (f), 29a to 30a", ["29a"], ["30a"],
    ["Stryker"], ["toluene"], simple_steps("Stryker-reagent allene formation", q(-10, "°C"), q(2, "h")),
    yield_out(["30a"], 53), ["DB-C-C", "C-H"], name="allene formation with Stryker's reagent",
    transformation="propargyl acetate to allene",
    notes=("Page 4 prose reports no evidence of reduced lactam methylene.",), extra_evidence=ev(4)
)


# Table 1: every row is a distinct optimization experiment.
TABLE1 = [
    (1, "23a", 1.8, 0, 1.5, 64, "100:trace:0", 110, ["25a", "32a"], "aged Stryker's reagent"),
    (2, "23a", 1.0, 0, 1.0, 73, "60:40:0", 110, ["25a", "32a"], "newly opened Stryker's reagent"),
    (3, "23a", 0.8, -10, 2.0, 66, "63:25:12", 110, ["25a", "32a"], "newly opened Stryker's reagent"),
    (4, "23a", 0.9, -20, 2.0, 58, "80:20:0", 110, ["25a", "32a"], "newly opened Stryker's reagent"),
    (5, "23a", 0.9, -10, 2.0, 78, "77:23:0", 220, ["25a", "32a"], "newly opened Stryker's reagent"),
    (6, "23a", 0.9, -10, 2.0, 49, "72:19:9", 440, ["25a", "32a"], "newly opened Stryker's reagent"),
    (7, "22a", 0.9, -10, 2.0, 70, "70:16:0", 110, ["24a", "31a", "cis_alkene"], "newly opened Stryker's reagent"),
    (8, "22a", 0.9, -10, 1.5, 56, "90:10:0", 122, ["24a", "31a"], "newly opened Stryker's reagent"),
]


for entry, sm, eqv, temp, tm, total, ratio, scale_mg, prods, bottle in TABLE1:
    def make_steps(e, r, temp=temp, tm=tm):
        r["solvent_conditions"]["concentration"] = extracted("0.04 M", q(0.04, "M"), e)
        r["solvent_conditions"]["concentration_basis"] = extracted("propargyl carboxylate", "starting propargyl pivalate", e)
        return [step("s1", e, "Stryker-reagent allene formation in degassed toluene with water", q(temp, "°C"), q(tm, "h"))]
    def make_outs(e, r, total=total, ratio=ratio, prods=prods, entry=entry, sm=sm):
        targets = prods
        outs = [outcome("o1", "yield", targets, q(total, "%"), e, role="total isolated mixture",
                        yield_type="total yield", scope="combined products reported in Table 1")]
        ratio_targets = prods + ([sm] if ratio.split(":")[-1] != "0" else [])
        outs.append(outcome("o2", "selectivity", ratio_targets, ratio, e, role="product/recovered-starting-material mixture",
                            scope=f"Table order (24a/25a)/(31a/32a)/SM; starting material={sm}"))
        if entry == 1:
            outs.append(outcome("o3", "trace", ["32a"], "trace", e, role="byproduct", scope="reduced methyl lactam"))
        if entry == 7:
            outs.append(outcome("o3", "selectivity", ["cis_alkene"], {"number": 14, "unit": "%", "qualifier": "approx", "raw_text": "~14% of the product mixture"}, e,
                                role="byproduct", scope="fraction of total product mixture from alkyne reduction"))
        return outs
    rr = reaction(
        f"rxn_p4_t1_e{entry}", 4, f"Table 1 entry {entry}", [sm], prods,
        [("Stryker", "reagent", eqv, None), ("water", "additive", 2, None)], ["toluene"],
        make_steps, make_outs, ["DB-C-C", "C-H"], category="optimization",
        name="optimization of allene formation with Stryker's reagent",
        transformation="propargyl pivalate to allene with competing methylene/alkyne reduction",
        scale=q(scale_mg, "mg"), notes=(bottle, "The table reports total yield and a composition ratio; no unsupported individual product yields were calculated.")
    )
    if entry in (3, 6):
        rr["review"]["unresolved"].append("Starting material is present in the composition ratio, but the source does not report an independent recovery percentage or measurement basis.")
    if entry == 7:
        rr["review"]["unresolved"].append("The cis-alkene byproduct is not assigned a compound number or full structure on pages 1-6.")


# Page 5, reactions of cis lactams only reported as an unassigned group range.
for inp, out in [("22b", "24b"), ("23b", "25b")]:
    rr = reaction(
        f"rxn_p5_cis_{inp}_to_{out}", 5, f"left-column continuation: cis lactam {inp} with Stryker's reagent",
        [inp], [out], ["Stryker"], [], simple_steps("reaction with Stryker's reagent", None, None),
        lambda e, r, out=out: [outcome("o1", "yield", [out], None, e, yield_type="individual yield not assigned", status="unreported",
                                             scope="source gives only a shared 30-35% range for 24b/25b experiments")],
        ["DB-C-C", "C-H"], category="substrate_scope", name="cis-lactam allene formation",
        transformation="propargyl pivalate to allene",
        notes=("The 30-35% range is reported collectively for the two cis-lactam experiments and was not assigned to either product.",)
    )
    rr["review"]["unresolved"].append("Individual yield, exact conditions, byproduct amount and measurement basis require the explicitly cited SI, which was not provided/associated.")


# Scheme 3: four arylations, one heteroarylation, tosylation and Boc protection.
for reagent_key, out, yld in [("33", "37", 44), ("34", "38", 42), ("35", "39", 47), ("36", "40", 63), ("41", "42", 59)]:
    reaction(
        f"rxn_p5_s3_24a_{reagent_key}_to_{out}", 5,
        f"Scheme 3 {'Eq 3' if reagent_key == '41' else 'Eq 2'}: 24a + {reagent_key} to {out}",
        ["24a", reagent_key], [out], [("CuI", "reagent", None, 20), ("DMEDA", "ligand", None, 40), "Cs2CO3"], ["toluene"],
        simple_steps("Cu-catalyzed N-arylation/N-heteroarylation", q(80, "°C"), q(20, "h")),
        yield_out([out], yld), ["C-N"], name="N-arylation of lactam", transformation="lactam N-C(aryl) bond formation"
    )


reaction(
    "rxn_p5_s3_24a_to_43", 5, "Scheme 3 Eq 4: N-tosylation of 24a", ["24a"], ["43"],
    ["NaH", "pTsCl"], ["DMF"], simple_steps("N-tosylation using NaH and p-TsCl", q(0, "°C"), q(2, "h")),
    yield_out(["43"], 29), ["N-S"], name="N-tosylation", transformation="lactam N-S bond formation"
)
reaction(
    "rxn_p5_s3_24a_to_44", 5, "Scheme 3 Eq 5: N-Boc protection of 24a", ["24a"], ["44"],
    ["Boc2O", "DMAP", "Et3N"], ["DCM"], zero_to_rt_steps("N-Boc protection", q(2, "h")),
    yield_out(["44"], 47, ytype="average yield"), ["C-N"], name="N-Boc protection", transformation="lactam N-C(carbamate) bond formation"
)


# Table 2 APKR scope.  Entry 2 has a 3 h post-addition hold; all other rows use
# the common 30 min hold.  Footnote c is retained as an alternative operation,
# not fabricated into a second row with an invented yield.
TABLE2 = [
    (1, "24a", "45", 72, "average of two runs", True, False),
    (2, "25a", "46", 46, "average of two runs", False, True),
    (3, "26", "47", 53, "average of two runs", False, False),
    (4, "30a", "48", 75, "one run", False, False),
    (5, "37", "49", 67, "average of two runs", True, False),
    (6, "38", "50", 72, "one run", False, False),
    (7, "39", "51", 70, "average of two runs", False, False),
    (8, "40", "52", 77, "average of two runs", False, False),
    (9, "42", "53", 67, "average of two runs", False, False),
    (10, "43", "54", 71, "average of two runs", True, False),
    (11, "44", "55", 40, "average of two runs", True, False),
]


for entry, sm, prod, yld, basis, single_similar, long_hold in TABLE2:
    def apkr_steps(e, r, long_hold=long_hold):
        return [
            step("s1", e, "add allene-yne dropwise to Rh(I) catalyst over 1 h", q(110, "°C"), q(1, "h"), atmosphere="CO", pressure="1 atm"),
            step("s2", e, "heat after completion of allene-yne addition", q(110, "°C"), q(3, "h") if long_hold else q(30, "min"), atmosphere="CO", pressure="1 atm"),
        ]
    notes = []
    if single_similar:
        notes.append("Table 2 footnote c states that similar yields were obtained with single-portion addition, but gives no separate exact yield.")
    reaction(
        f"rxn_p6_t2_e{entry}", 6, f"Table 2 entry {entry}: {sm} to {prod}", [sm], [prod],
        [("Rhcat", "catalyst", None, 10), "CO"], ["toluene"], apkr_steps,
        yield_out([prod], yld, ytype="APKR yield", basis=basis), ["C-C", "DB-C-O"],
        category="substrate_scope", name="allenic Pauson-Khand reaction (APKR)",
        transformation="intramolecular allene-yne-carbon monoxide [2+2+1] cycloaddition", notes=tuple(notes), extra_evidence=ev(5)
    )


# Scheme 4 methyl-lactam APKR experiments.
reaction(
    "rxn_p6_s4_eq7_56_to_58", 6, "Scheme 4 Eq 7: compound 56 to 58", ["56"], ["58"],
    [("Rhcat", "catalyst", None, 10), "CO"], ["toluene"],
    simple_steps("APKR under CO (1 atm)", q(110, "°C"), q(1.5, "h"), atmosphere="CO", pressure="1 atm"),
    lambda e, r: [
        outcome("o1", "yield", ["58"], q(78, "%"), e, yield_type="reported yield"),
        outcome("o2", "selectivity", ["58"], "7:1", e, scope="product dr shown in Scheme 4"),
    ], ["C-C", "DB-C-O"], name="allenic Pauson-Khand reaction (APKR)",
    transformation="intramolecular allene-yne-carbon monoxide [2+2+1] cycloaddition",
    notes=("Page 6 text describes the product as an 8:1 mixture, while Scheme 4 labels it dr 7:1; the discrepancy is retained.",)
)
REACTIONS[-1]["outcomes"][1]["value"] = conflict_fact([
    ("7:1 dr (Scheme 4)", "7:1", ev(6)), ("8:1 mixture (page 6 prose)", "8:1", ev(6))
], ev(6))
REACTIONS[-1]["review"]["unresolved"].append("Resolve source discrepancy between 7:1 (scheme) and 8:1 (prose) for compound 58.")

reaction(
    "rxn_p6_s4_eq8_57_to_59", 6, "Scheme 4 Eq 8: compound 57 to 59", ["57"], ["59"],
    [("Rhcat", "catalyst", None, 10), "CO"], ["toluene"],
    simple_steps("single-portion addition followed by APKR under CO (1 atm)", q(110, "°C"), q(1, "h"), atmosphere="CO", pressure="1 atm"),
    lambda e, r: [outcome("o1", "yield", ["59"], q(77, "%"), e, yield_type="reported yield"),
                  outcome("o2", "selectivity", ["59"], "single diastereomer", e, scope="product diastereoselectivity; C11 assigned by coupling-constant comparison")],
    ["C-C", "DB-C-O"], name="allenic Pauson-Khand reaction (APKR)",
    transformation="intramolecular allene-yne-carbon monoxide [2+2+1] cycloaddition"
)


# Page 6 cysteamine reactivity experiments.  Half-life is an outcome directed
# at the tested electrophile.  No isolated thiol-adduct structure is invented.
KINETICS = [
    ("45", q(2.1, "d"), 6), ("49", q(2.1, "d"), 6), ("51", q(9.8, "h"), 6),
    ("52", q(8.5, "d"), 6), ("54", "complete in less than 10 min; half-life not determined", 6),
    ("60", q(5, "min"), 6), ("S17", q(7, "min"), 7),
]


for comp, result, result_page in KINETICS:
    extra = ev(result_page, image=False) if result_page != 6 else []
    def kinetic_outs(e, r, comp=comp, result=result):
        if comp == "54":
            return [outcome("o1", "other", [comp], result, e, role="tested electrophile",
                            basis="1H NMR disappearance of exocyclic methylene resonances",
                            scope="time to complete reaction; half-life explicitly not determined")]
        return [outcome("o1", "other", [comp], result, e, role="tested electrophile",
                        basis="pseudo-first-order reaction progress monitored by 1H NMR",
                        scope="reaction half-life with cysteamine")]
    rr = reaction(
        f"rxn_p6_cysteamine_{comp}", 6, f"page 6 right column: cysteamine reactivity of compound {comp}",
        [comp], [], [("cysteamine", "reagent", 15, None)], ["CDCl3"],
        simple_steps("mix electrophile and cysteamine; acquire 1H NMR spectra at regular intervals", q(22, "°C"), None),
        kinetic_outs, ["C-S", "C-H"], category="mechanistic_study",
        name="cysteamine pseudo-first-order thiol-reactivity experiment",
        transformation="hetero-Michael addition to alpha-methylene-gamma-lactam/lactone",
        notes=("The thiol adduct was not isolated or structurally registered from pages 1-6.",), extra_evidence=extra
    )
    # The source lists lactam plus cysteamine; avoid treating the second mention
    # as two independent chemical additions in the participant registry.
    seen = set()
    dedup = []
    for p in rr["participants"]:
        key = (p["molecule_id"], p["role"], p["stage"])
        if key in seen:
            continue
        seen.add(key)
        dedup.append(p)
    rr["participants"] = dedup
    if comp == "S17":
        rr["review"]["unresolved"].append("The 7 min half-life is taken from extracted page-7 continuation text; page 7 image was outside the requested visual-coverage range.")


def source_context(context_id, page, region, classification_name, text, *, extra=()):
    evidence = list(dict.fromkeys(ev(page) + list(extra)))
    return {
        "context_id": context_id,
        "document_id": DOC,
        "location": {"page": page, "region": region},
        "classification": classification_name,
        "summary": text,
        "evidence_ids": evidence,
        "review_status": "needs_review",
    }


contexts = [
    source_context("ctx_p1_bibliography", 1, "article header", "bibliographic_metadata",
                   "Title: Synthesis of Guaianolide Analogues with a Tunable α-Methylene−γ-lactam Electrophile and Correlating Bioactivity with Thiol Reactivity. DOI: 10.1021/acs.jmedchem.0c01464."),
    source_context("ctx_p1_graphical_abstract", 1, "graphical abstract beneath access bar", "generic_scheme",
                   "Generic APKR followed by thiol addition; summary graphic, not a separately reported experiment."),
    source_context("ctx_p1_intro", 1, "abstract and introduction columns", "background",
                   "Title, abstract and conceptual motivation; no page-1 specific synthetic experiment counted."),
    source_context("ctx_p2_fig1", 2, "Figure 1", "proposed_mechanism",
                   "Conceptual noncovalent binding followed by covalent thiol addition; not an experimental reaction record."),
    source_context("ctx_p2_fig2", 2, "Figure 2A-B", "literature_example",
                   "Compound 1 and retrosynthetic/literature comparison 4->2 versus proposed 5->3; not counted as current experiments."),
    source_context("ctx_p2_scheme1", 2, "Scheme 1", "generic_scheme",
                   "Retrosynthetic APKR plan 6<-7<-8<-9+10; no yield or specific run."),
    source_context("ctx_p2_precursor_history", 2, "right-column prose", "unavailable_si_context",
                   "Compound 11 is said to originate from S1 in three steps in SI; SI association is unestablished, so those steps were not invented."),
    source_context("ctx_p3_ts16_17", 3, "Scheme 2 structures 16 and 17 plus discussion", "proposed_mechanism",
                   "Zimmerman-Traxler/alternative transition-state rationale for trans/cis ratios; mechanistic proposal, not a separate experiment."),
    source_context("ctx_p4_acetylation_observation", 4, "left-column pivaloylation discussion", "unresolved_observation",
                   "Acetylation conditions are stated to give an N-acetylated product, but no product label, amount or yield is provided on pages 1-6."),
    source_context("ctx_p6_fig3", 6, "Figure 3", "characterization",
                   "Crystal structure of APKR product 45; characterization evidence, not a synthetic experiment."),
    source_context("ctx_p6_cbeta", 6, "right-column electrophilic-reactivity discussion", "characterization",
                   "Reported 13C-beta chemical shifts for compounds 45, 49, 51, 52, 54 and 60 are characterization/property context, not independent syntheses."),
]


template = json.loads((RUN / "submitted.json").read_text())
# The document/evidence block is immutable after `standard init`.  Re-read it
# from the prepared dataset so repeated deterministic builds cannot drift.
template["documents"] = json.loads((RUN / "dataset.json").read_text())["documents"]
template["molecules"] = list(MOLECULES.values())
template["reactions"] = REACTIONS
template["source_context"] = contexts
(RUN / "submitted.json").write_text(json.dumps(template, ensure_ascii=False, indent=2) + "\n")


page_record_ids = {i: [] for i in range(1, 7)}
for r in REACTIONS:
    for eid in r["evidence_ids"]:
        if eid.startswith(f"{DOC}:page"):
            try:
                p = int(eid.split(":page", 1)[1].split(":", 1)[0])
            except ValueError:
                continue
            if p in page_record_ids and r["reaction_id"] not in page_record_ids[p]:
                page_record_ids[p].append(r["reaction_id"])


coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {
            "page": 1, "state": "attempted",
            "regions": [
                {"region": "graphical abstract", "finding": "generic APKR and thiol-addition summary; not a specific experiment", "evidence_ids": ev(1)},
                {"region": "abstract/introduction", "finding": "background and hypotheses; no independent page-1 experiment", "evidence_ids": ev(1)},
            ],
            "record_ids": page_record_ids[1],
            "unresolved": ["Graphical structures are generic R-group depictions and were not promoted to concrete molecule SMILES."],
        },
        {
            "page": 2, "state": "attempted",
            "regions": [
                {"region": "Figure 1", "finding": "conceptual covalent-inhibitor mechanism, not a specific experiment", "evidence_ids": ev(2)},
                {"region": "Figure 2A-B", "finding": "compound example and cited/prior APKR comparison; literature context", "evidence_ids": ev(2)},
                {"region": "Scheme 1", "finding": "retrosynthetic generic plan, not an executed run", "evidence_ids": ev(2)},
                {"region": "right-column final paragraph", "finding": "actual failed 3-trimethylsilyl-2-propynal lactamization substrate test", "evidence_ids": ev(2)},
            ],
            "record_ids": page_record_ids[2],
            "unresolved": ["Three-step preparation of 11 is only referred to SI; SI was not provided/associated.", "Failed TMS-aldehyde test has no shown product and no independently reported conversion."],
        },
        {
            "page": 3, "state": "attempted",
            "regions": [
                {"region": "Scheme 2A", "finding": "each synthetic arrow transcribed; compound 9 yield/Z:E conflicts preserved", "evidence_ids": ev(3)},
                {"region": "Scheme 2B", "finding": "multistep sequence split at isolated intermediates/yields", "evidence_ids": ev(3)},
                {"region": "structures 16/17 and lower discussion", "finding": "proposed transition states and stereochemical rationale, not independent experiments", "evidence_ids": ev(3)},
                {"region": "Scheme 2 footnote", "finding": "ordered reagent additions, temperatures, holds and ramp statements assigned to reaction steps", "evidence_ids": ev(3)},
            ],
            "record_ids": page_record_ids[3],
            "unresolved": ["Full atom-mapped/isomeric SMILES were not established.", "The meaning of the printed 5-2:1 and 3.5-2:1 ranges remains as source text.", "No SI workup/purification details were inferred."],
        },
        {
            "page": 4, "state": "attempted",
            "regions": [
                {"region": "Table 1 including Eq 1 and footnotes", "finding": "all eight optimization rows separately recorded; aged/new reagent and cis-alkene footnotes retained", "evidence_ids": ev(4)},
                {"region": "left-column Scheme 2 continuation", "finding": "cross-page yields/selectivities linked to existing Scheme 2 records", "evidence_ids": ev(4)},
                {"region": "right-column optimization discussion", "finding": "0.04 M, 2 equiv water and composition details associated with table records", "evidence_ids": ev(4)},
            ],
            "record_ids": page_record_ids[4],
            "unresolved": ["Table 1 does not state the measurement basis for total yield.", "Entries 3 and 6 include SM in a ratio but do not report an independent recovery percentage.", "Entry 7 cis-alkene byproduct lacks a compound number/full structure."],
        },
        {
            "page": 5, "state": "attempted",
            "regions": [
                {"region": "Scheme 3 Eq 2", "finding": "four N-arylation experiments separately recorded", "evidence_ids": ev(5)},
                {"region": "Scheme 3 Eq 3-5 and footnote", "finding": "heteroarylation, N-tosylation and N-Boc protection separately recorded with footnote reagents", "evidence_ids": ev(5)},
                {"region": "left-column continuation above heading", "finding": "two cis-lactam Stryker tests retained without assigning shared 30-35% range to either product", "evidence_ids": ev(5)},
                {"region": "right-column APKR discussion", "finding": "conditions and outcomes linked forward to Table 2 records on page 6", "evidence_ids": ev(5)},
            ],
            "record_ids": page_record_ids[5],
            "unresolved": ["Individual yields for 24b and 25b cannot be assigned from the shared 30-35% statement.", "No SI details were imported."],
        },
        {
            "page": 6, "state": "attempted",
            "regions": [
                {"region": "Table 2 and Eq 6 including footnotes", "finding": "all eleven APKR scope rows separately recorded; addition/hold times and run basis preserved", "evidence_ids": ev(6)},
                {"region": "Scheme 4 Eq 7-8", "finding": "two methyl-lactam APKR experiments recorded; 58 dr conflict retained", "evidence_ids": ev(6)},
                {"region": "Figure 3", "finding": "X-ray characterization of 45, not a separate synthesis", "evidence_ids": ev(6)},
                {"region": "electrophilic-reactivity paragraph", "finding": "seven cysteamine kinetic experiments represented; S17 numerical continuation uses page-7 text only", "evidence_ids": ev(6)},
            ],
            "record_ids": page_record_ids[6],
            "unresolved": ["S17 half-life continues onto page 7, whose image was outside this visual-coverage batch.", "Thiol adduct structures were not isolated on pages 1-6.", "13C shifts were kept as characterization context, not reaction records."],
        },
    ],
}
(RUN / "coverage").mkdir(parents=True, exist_ok=True)
(RUN / "coverage" / "pages_0001_0006.json").write_text(json.dumps(coverage, ensure_ascii=False, indent=2) + "\n")


inventory = {
    "batch_id": "allene-stereo-v31-200-20261017",
    "document_id": DOC,
    "visually_inspected_pages": [1, 2, 3, 4, 5, 6],
    "reaction_candidate_count": len(REACTIONS),
    "molecule_candidate_count": len(MOLECULES),
    "reaction_ids": [r["reaction_id"] for r in REACTIONS],
    "source_context_ids": [c["context_id"] for c in contexts],
    "limitations": [
        "Machine records remain needs_review.",
        "No SI was associated.",
        "Page 7 text was used only to complete the S17 sentence crossing the requested page boundary; page 7 was not visually covered.",
        "No three-model vote was performed.",
    ],
}
(RUN / "raw_model" / "page_inventory.json").write_text(json.dumps(inventory, ensure_ascii=False, indent=2) + "\n")
(RUN / "raw_model" / "standard_stereo_checks.json").write_text(json.dumps({
    key: {
        "input_smiles": data[0],
        "canonical_isomeric_smiles": data[0],
        "stereo_graph_key": data[1],
        "atom_centers": data[2],
        "double_bonds": data[3],
        "source_verified": False,
        "limitations": [
            "Graph diagnostics cannot infer unreported configuration, ee, dr or racemic composition.",
            "Axial/allene, relative and ring cis/trans assignments need source review; toolkit support is limited.",
        ],
    }
    for key, data in STEREO_CHECKS.items()
}, ensure_ascii=False, indent=2) + "\n")
print(json.dumps({"molecules": len(MOLECULES), "reactions": len(REACTIONS), "source_context": len(contexts)}, ensure_ascii=False))
