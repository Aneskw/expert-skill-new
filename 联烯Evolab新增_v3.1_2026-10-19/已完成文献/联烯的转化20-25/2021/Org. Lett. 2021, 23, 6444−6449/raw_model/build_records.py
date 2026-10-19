#!/usr/bin/env python3
"""Build source-derived v3.1 candidate records for doc_0d74203a7b8cf397.

This is a paper-specific extraction artifact. It does not certify chemical accuracy.
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
DOC = "doc_0d74203a7b8cf397"
EV = {n: f"{DOC}:page{n}:image" for n in range(1, 7)}
TXT = {n: f"{DOC}:page{n}:text" for n in range(1, 7)}


def F(raw=None, normalized=None, evidence=(), method="original_extraction", derivation=None,
      status=None):
    return st.fact(raw, normalized, evidence, method=method, derivation=derivation, status=status)


def q(number, unit, raw_text=None, qualifier=None):
    value = {"number": number, "unit": unit}
    if raw_text is not None:
        value["raw_text"] = raw_text
    if qualifier is not None:
        value["qualifier"] = qualifier
    return value


bundle = json.loads((RUN / "submitted.json").read_text())
molecules = {}
stereo_diagnostics = {}


def add_molecule(key, label, name=None, evidence=(EV[1],), smiles=None,
                 chirality="unresolved", geometry="unresolved",
                 representation=None, assignment=None, atom_centers=None,
                 double_bonds=None, mixture=None, identity_note=None,
                 generic=None):
    mid = f"{DOC}:mol:{key}"
    m = st.new_molecule(mid, DOC, evidence)
    m["occurrences"][0]["compound_label"] = F(label, label, evidence)
    m["names"]["original"] = F(label, label, evidence)
    if name:
        m["names"]["english"] = F(name, name, evidence)
    m["cas_numbers"] = F(evidence=evidence)
    if generic is not None:
        m["representation_type"] = "generic_template"
        m["smiles"] = F(evidence=evidence)
        variables = []
        for var in generic.get("variables", []):
            variables.append({
                "label": var,
                "attachment_points": F(f"{var} attachment shown in source drawing", [var], evidence),
                "definition": F(f"Source variable {var}", var, evidence),
                "allowed_values": F(evidence=evidence),
                "constraints": F(evidence=evidence),
            })
        m["generic_template"] = {
            "fixed_structure": F(generic["raw"], generic["fixed"], evidence),
            "variables": variables,
            "instance_mapping": [],
            "reason": F("Generic R-group drawing", generic["reason"], evidence),
        }
        representation = representation or "not_encoded"
    elif smiles:
        m["representation_type"] = "smiles"
        deriv = "Candidate connectivity/isomeric SMILES transcribed from the labelled source structure; RDKit parsing is only a graph diagnostic."
        m["smiles"] = F(f"source drawing/label {label}", smiles, evidence,
                         method="structure_derivation", derivation=deriv)
        diag = st.stereo_inspect(smiles)
        stereo_diagnostics[mid] = diag
        m["stereochemistry"]["isomeric_smiles"] = F(
            f"source drawing/label {label}", diag["canonical_isomeric_smiles"], evidence,
            method="structure_derivation", derivation=deriv)
        m["stereochemistry"]["identity_key"] = F(
            None, diag["stereo_graph_key"], evidence, method="structure_derivation",
            derivation="Computed from the candidate isomeric molecular graph; not a sample identity or source verification.")
        representation = representation or ("partial" if any(
            p["specified"] == "Unspecified" for p in diag["potential_stereo"]) else "encoded")
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = F(evidence=evidence)
        representation = representation or "unresolved"
    for field, value in (("chirality_status", chirality), ("geometry_status", geometry),
                         ("representation_status", representation)):
        m["stereochemistry"][field] = F(
            None, value, evidence, method="inference",
            derivation="Status assigned from the source drawing/name and the limits of the retained molecular representation.")
    if assignment is not None:
        m["stereochemistry"]["assignment"] = F(assignment, assignment, evidence)
    if atom_centers is not None:
        m["stereochemistry"]["atom_centers"] = F(atom_centers, atom_centers, evidence)
    if double_bonds is not None:
        m["stereochemistry"]["double_bonds"] = F(double_bonds, double_bonds, evidence)
    if mixture is not None:
        m["stereochemistry"]["mixture"] = F(mixture, mixture, evidence)
    m["identity_note"] = F(
        identity_note or "Identity retained at the precision supported by the supplied main-text drawing; no CAS or external lookup used.",
        identity_note or "Identity retained at the precision supported by the supplied main-text drawing; no CAS or external lookup used.", evidence)
    molecules[key] = m
    return mid


def mid(key):
    return molecules[key]["molecule_id"]


# Generic concepts and the four arrows in Figure 1.
add_molecule("generic_1", "1", "generic amino-substituted allylmetal reagent", evidence=(EV[1],),
             generic={"raw": "Figure 1A: L_nM-allyl-NR2", "fixed": {"class": "alpha,gamma-aminoanion allylmetal equivalent"},
                      "variables": ["M", "L_n", "R"], "reason": "Figure 1A is a variable conceptual reagent."})
add_molecule("generic_2", "2", "generic carbonyl electrophile", evidence=(EV[1],),
             generic={"raw": "Figure 1A: R1-C(=O)-R2", "fixed": {"class": "carbonyl"},
                      "variables": ["R1", "R2"], "reason": "Figure 1A contains undefined carbonyl substituents."})
add_molecule("generic_3", "3", "generic aminoallylation product", evidence=(EV[1],),
             chirality="unspecified", geometry="not_applicable",
             generic={"raw": "Figure 1A product 3", "fixed": {"class": "alkenyl 1,2-aminoalcohol synthon"},
                      "variables": ["R", "R1", "R2"], "reason": "Figure 1A is a generic product drawing."})
add_molecule("generic_4", "4", "generic 1,2-aminoalcohol", evidence=(EV[1],),
             chirality="unspecified", geometry="not_applicable",
             generic={"raw": "Figure 1A product 4", "fixed": {"class": "vicinal aminoalcohol"},
                      "variables": ["R", "R1", "R2"], "reason": "Figure 1A is a retrosynthetic motif."})
add_molecule("generic_8", "8", "generic ketone 8", evidence=(EV[1], EV[3]),
             generic={"raw": "Figure 1C/D and Scheme 1: R1-C(=O)-R2", "fixed": {"class": "ketone"},
                      "variables": ["R1", "R2"], "reason": "Current-work scheme uses undefined ketone substituents."})
add_molecule("generic_12", "12", "generic protected aminoallylation product 12", evidence=(EV[1], EV[3]),
             chirality="specified", geometry="not_applicable", representation="not_encoded",
             assignment="Two stereogenic relationships are depicted, but R1/R2 are variable and absolute descriptors are not printed.",
             atom_centers="Two depicted stereogenic centers in the generic product.",
             generic={"raw": "Figure 1D and Scheme 1 generic product 12", "fixed": {"class": "protected alkenyl 1,2-aminoalcohol"},
                      "variables": ["R1", "R2"], "reason": "Current-work product is a generic R-group drawing."})
for key, label, name in [
    ("lit_5", "5", "primary alcohol"), ("lit_6", "6", "phthalimido-allene"),
    ("lit_7", "7", "Ir-catalyzed aminoallylation product"),
    ("lit_9", "9", "Evans-auxiliary-derived chiral allenamide"),
    ("lit_10", "10", "auxiliary-controlled aminoallylation product")]:
    add_molecule(key, label, name, evidence=(EV[1],), chirality="specified" if key in {"lit_9", "lit_10"} else "unresolved",
                 geometry="not_applicable", representation="not_encoded" if key in {"lit_9", "lit_10"} else "unresolved",
                 assignment="Source drawing specifies stereochemistry but no complete descriptor was transcribed." if key in {"lit_9", "lit_10"} else None)
add_molecule("ir_cat", "(R)-H8-BINAP-Ir", "(R)-H8-BINAP iridium catalyst", evidence=(EV[1],),
             chirality="specified", representation="not_encoded", assignment="(R) ligand designation in Figure 1B")
add_molecule("simes", "SIMES", "SIMES ligand", evidence=(EV[1],))

# Current reaction reagents, solvents, ligands and workup agents.
add_molecule("11", "11", "achiral allenamide 11", evidence=(EV[1], EV[2], EV[3], EV[4]),
             smiles="C=C=CN1CCOC1=O", chirality="achiral", geometry="not_applicable",
             identity_note="Candidate N-allenyl oxazolidinone connectivity; terminal allene is not axially chiral. Exact graph remains needs-review against the drawing.")
reagent_evidence = {
    "cuoac2": (EV[1], EV[2], EV[3], EV[4]), "silane": (EV[1], EV[2], EV[3], EV[4]),
    "toluene": (EV[1], EV[2], EV[3]), "phcf3": (EV[2], EV[3], EV[4]),
    "tbuoh": (EV[2], EV[3], EV[4]), "nh4f": (EV[1], EV[2], EV[3], EV[4]),
    "meoh": (EV[3],), "nah": (EV[3], EV[4]), "thf": (EV[4],), "cui": (EV[4],),
    "kotbu": (EV[4],), "tscl": (EV[4],), "dbu": (EV[4],), "nai": (EV[4],),
    "h2so4": (EV[4],), "hexane": (EV[4],), "etoac": (EV[4],),
}
for key, label, name, smi in [
    ("cuoac2", "Cu(OAc)2", "copper(II) acetate", None),
    ("silane", "Me(MeO)2SiH", "methyldimethoxysilane", "CO[SiH](C)OC"),
    ("toluene", "toluene", "toluene", "Cc1ccccc1"),
    ("phcf3", "PhCF3", "benzotrifluoride", "FC(F)(F)c1ccccc1"),
    ("tbuoh", "t-BuOH", "tert-butanol", "CC(C)(C)O"),
    ("nh4f", "NH4F", "ammonium fluoride", "[NH4+].[F-]"),
    ("meoh", "MeOH", "methanol", "CO"),
    ("nah", "NaH", "sodium hydride", "[Na+].[H-]"),
    ("thf", "THF", "tetrahydrofuran", "C1CCOC1"),
    ("cui", "CuI", "copper(I) iodide", "[Cu]I"),
    ("kotbu", "KOtBu", "potassium tert-butoxide", "CC(C)(C)[O-].[K+]"),
    ("tscl", "TsCl", "p-toluenesulfonyl chloride", "Cc1ccc(S(=O)(=O)Cl)cc1"),
    ("dbu", "DBU", "1,8-diazabicyclo[5.4.0]undec-7-ene", None),
    ("nai", "NaI", "sodium iodide", "[Na+].[I-]"),
    ("h2so4", "H2SO4", "sulfuric acid", "OS(=O)(=O)O"),
    ("hexane", "hex", "hexanes", None), ("etoac", "EtOAc", "ethyl acetate", "CCOC(C)=O")]:
    add_molecule(key, label, name, evidence=reagent_evidence[key], smiles=smi,
                 chirality="achiral" if smi else "unresolved", geometry="not_applicable",
                 representation="encoded" if smi else "unresolved")

ligands = {
    "ss_phbpe": ("(S,S)-Ph-BPE", "(S,S)-Ph-BPE"), "r_binap": ("(R)-BINAP", "(R)-BINAP"),
    "r_segphos": ("(R)-Segphos", "(R)-Segphos"), "w3": ("W3", "Walphos W3"),
    "w8": ("W8", "Walphos W8"), "j6": ("J6", "Josiphos J6"), "j11": ("J11", "Josiphos J11")}
for key, (label, name) in ligands.items():
    add_molecule(key, label, name, evidence=(EV[1], EV[2], EV[3], EV[4]) if key == "w8" else (EV[2],),
                 chirality="specified" if label.startswith("(") else "partially_specified",
                 geometry="not_applicable", representation="not_encoded",
                 assignment=(label + " source designation") if label.startswith("(") else "Chiral ligand drawing/usage is explicit; full absolute descriptor is not printed in the supplied main text.")

# Ketone scope. Concrete graphs are candidates only; dense/ambiguous 8r remains unresolved.
ketone_smiles = {
    "8a": "CC(=O)c1ccccc1", "8b": "CCC(=O)c1ccccc1", "8c": "O=C(Cc1ccccc1)c1ccccc1",
    "8d": "CC(=O)c1ccc(C(F)(F)F)cc1", "8e": "CC(=O)c1ccc(Cl)cc1", "8f": "COc1ccc(C(C)=O)cc1",
    "8g": "CC(=O)c1ccc(N(C)C)cc1", "8h": "CC(=O)c1ccc(C)cc1", "8i": "CC(=O)c1ccc2OCOc2c1",
    "8j": "CC(=O)c1ccc2ccccc2c1", "8k": "CC(=O)c1cccc(C(F)(F)F)c1", "8l": "COC(=O)c1cccc(C(C)=O)c1",
    "8m": "CC(=O)c1cccc(Br)c1", "8n": "COc1cccc(C(C)=O)c1", "8o": "CC(=O)c1cccc(N(C)C)c1",
    "8p": "CC(=O)c1ccco1", "8q": "CC(=O)c1ccc(C)s1", "8s": "O=C1CCCc2ccccc21",
    "8t": "CC(C)=O", "8u": "O=C(c1ccccc1)c1ccccc1"}
for letter in "abcdefghijklmnopqrstu":
    key = "8" + letter
    smi = ketone_smiles.get(key)
    add_molecule(key, key, f"ketone {key}", evidence=(EV[2], EV[3], EV[4]) if key in {"8a", "8b"} else (EV[3],),
                 smiles=smi, chirality="achiral" if smi else "unresolved", geometry="not_applicable",
                 representation="encoded" if smi else "unresolved",
                 identity_note=("Candidate connectivity derived from the Scheme 1 substituent drawing; positional assignment needs expert review."
                                if smi else "Dense N-tosyl heteroaryl ketone drawing retained unresolved rather than assigning unsafe connectivity."))

# Protected products, migration products and application/mechanistic products.
for letter in "abcdefghijklmnopqrst":
    add_molecule("12" + letter, "12" + letter, f"branched protected 1,2-aminoalcohol 12{letter}",
                 evidence=(EV[2], EV[3], EV[4]) if letter in "ab" else (EV[3],),
                 chirality="specified", geometry="not_applicable", representation="not_encoded",
                 assignment="Two stereogenic centers are depicted by solid/dashed stereobonds; absolute R/S descriptors are not printed in Scheme 1.",
                 atom_centers="Two depicted stereogenic centers; exact atom indices unavailable because complete SMILES is unresolved.")
add_molecule("l12u", "l-12u", "linear protected aminoallylation product l-12u", evidence=(EV[3],),
             chirality="unspecified", geometry="unspecified", representation="not_encoded",
             assignment="Scheme 1 does not print an absolute configuration or er for l-12u.")
add_molecule("l12s", "linear isomer of 12s", "linear isomer accompanying 12s", evidence=(EV[3],),
             chirality="unspecified", geometry="unspecified", representation="unresolved")
for key in ("13a", "13b", "13c"):
    add_molecule(key, key, f"carbonate-migration product {key}", evidence=(EV[2], EV[3], EV[4]),
                 chirality="specified", geometry="not_applicable", representation="not_encoded",
                 assignment="Source drawing/result identifies a single depicted diastereomer; absolute descriptors are not fully encoded.",
                 atom_centers="Stereocenters depicted in source; complete atom-index assignment unresolved.")
for key, label, name, chiral in [
    ("b17", "b-17", "branched allylation product b-17", True),
    ("l17", "l-17", "linear allylation product l-17", False),
    ("18", "18", "protonation product 18", False), ("19", "19", "protonation product 19", False),
    ("20", "20", "oxazolidinone 20", True), ("21", "21", "aminoalcohol 21", True),
    ("lcuotbu", "LCuOtBu", "ligated copper tert-butoxide catalyst", False)]:
    add_molecule(key, label, name, evidence=(EV[4],), chirality="specified" if chiral else "unresolved",
                 geometry="unspecified" if key in {"b17", "l17", "18", "19", "20", "21"} else "not_applicable",
                 representation="not_encoded" if chiral else "unresolved",
                 assignment="Source drawing/er specifies a chiral product but no complete absolute descriptor was transcribed." if chiral else None)


def participant(pid, key, role, evidence, stage="reaction", charges=None, ref=None, product_role=None):
    p = st.new_participant(pid, mid(key), role, evidence, stage=stage)
    if product_role:
        p["product_role"] = F(product_role, product_role, evidence)
    for field, value in (charges or {}).items():
        raw = value.get("raw_text") if isinstance(value, dict) and value.get("raw_text") else value
        p["charges"][field] = F(raw, value, evidence)
    if ref:
        p["equivalence_reference"] = F(ref, mid(ref), evidence)
    return p


LABEL_NUM = {row["code"]: row["number"] for row in json.loads(
    (ROOT / "skills/chemical-reaction-extraction/references/bond_taxonomy.json").read_text())["categories"]}


def classification(r, codes, evidence, limitations):
    if codes is None:
        r["classification"]["labels"] = F(None, [], evidence, method="inference",
                                                derivation="No defensible formed-bond assignment was made from the available main-text drawing.")
        r["classification"]["coverage"] = F(None, "unresolved", evidence, method="inference",
                                                  derivation="Structures or atom correspondence are insufficient for bond-change classification.")
        r["classification"]["limitations"] = limitations
        return
    changes, labels = [], []
    atoms_for = {"C-C": ("C", "C", "formed"), "C-H": ("C", "H", "hydrogen_gain"),
                 "O-H": ("O", "H", "hydrogen_gain"), "C-O": ("C", "O", "formed"),
                 "N-H": ("N", "H", "hydrogen_gain")}
    for i, code in enumerate(codes, 1):
        a, b, ctype = atoms_for[code]
        cid = f"b{i}"
        change = {"change_id": cid, "atoms": [{"atom_map": None, "element": a}, {"atom_map": None, "element": b}],
                  "before_order": 0, "after_order": 1, "change_type": ctype}
        if ctype == "hydrogen_gain":
            change["hydrogen_count"] = 1
        changes.append(F(None, change, evidence, method="structure_derivation",
                         derivation="Category-level net change inferred from the depicted reactant/product classes; atom-complete mapping is unavailable."))
        labels.append({"category_number": LABEL_NUM[code], "code": code, "supporting_change_ids": [cid]})
    r["classification"]["bond_changes"] = changes
    r["classification"]["labels"] = F(None, labels, evidence, method="structure_derivation",
                                            derivation="Labels follow all retained category-level formed-bond changes and the supplied 41-class taxonomy.")
    r["classification"]["coverage"] = F(None, "partial", evidence, method="inference",
                                              derivation="No atom-mapped complete structures; implicit-H gains are net changes only.")
    r["classification"]["limitations"] = limitations


def reaction(rid, record_type, categories, evidence, location, name, transformation,
             unresolved, codes=None, general=()):
    r = st.new_reaction(f"{DOC}:{rid}", DOC, evidence)
    r["record_type"] = record_type
    r["experiment_categories"] = F(categories, categories, evidence)
    r["location"] = F(location, location, evidence)
    r["general_scheme_ids"] = list(general)
    r["reaction_name"] = F(name, name, evidence)
    r["transformation_type"] = F(None, transformation, evidence, method="inference",
                                     derivation="Normalized transformation type inferred from the source caption/drawing.")
    r["review"]["unresolved"] = list(unresolved)
    classification(r, codes, evidence, [
        "Category-level inference only; complete source structures and atom mapping are unavailable.",
        "Implicit-H gains do not identify the elementary proton/hydride source."] if codes is not None else ["Classification left unresolved rather than inventing an atom mapping."])
    return r


def step(sid, evidence, participant_ids, operation, temperature=None, time=None, stage="reaction"):
    s = st.new_step(sid, evidence, stage=stage)
    s["participant_ids"] = list(participant_ids)
    s["operation"] = F(operation, operation, evidence)
    if temperature is not None:
        s["temperature"] = F(temperature if isinstance(temperature, str) else temperature.get("raw_text", temperature), temperature, evidence)
    if time is not None:
        s["time"] = F(time if isinstance(time, str) else time.get("raw_text", time), time, evidence)
    return s


def outcome(oid, kind, targets, raw, normalized, evidence, role=None, basis=None, scope=None, yield_type=None):
    o = st.new_outcome(oid, kind, [mid(t) for t in targets], evidence)
    o["value"] = F(raw, normalized, evidence)
    if role:
        o["product_role"] = F(role, role, evidence)
    if basis:
        o["measurement_basis"] = F(basis, basis, evidence)
    if scope:
        o["scope"] = F(scope, scope, evidence)
    if yield_type:
        o["yield_type"] = F(yield_type, yield_type, evidence)
    return o


reactions = []

# Figure 1 arrows: prior/conceptual arrows are never counted as current experiments.
for rid, inputs, prod, name in [
    ("fig1a_concept", ["generic_1", "generic_2"], "generic_3", "Aminoallylation approach to chiral vicinal aminoalcohols"),
    ("fig1b_krische", ["lit_5", "lit_6", "ir_cat"], "lit_7", "Ir-catalyzed enantioselective aldehyde aminoallylation"),
    ("fig1c_prior_auxiliary", ["generic_8", "lit_9", "cuoac2", "simes", "silane"], "lit_10", "Cu-catalyzed auxiliary-based ketone aminoallylation")]:
    r = reaction(rid, "literature_example", ["other"], (EV[1],), {"page": 1, "figure": "Figure 1"}, name,
                 "literature/contextual aminoallylation", ["Literature/conceptual arrow only; not an experiment performed in this paper."], None)
    r["participants"] = [participant(f"p{i+1}", x, "reactant" if x not in {"ir_cat", "cuoac2", "simes", "silane"} else ("catalyst" if x in {"ir_cat", "cuoac2"} else "reagent"), (EV[1],)) for i, x in enumerate(inputs)]
    r["participants"].append(participant("product", prod, "product", (EV[1],), product_role="depicted literature/concept product"))
    reactions.append(r)

GEN = f"{DOC}:fig1d_current_generic"
r = reaction("fig1d_current_generic", "generic_scheme", ["other"], (EV[1], EV[3]),
             {"pages": [1, 3], "regions": ["Figure 1D", "Scheme 1 header"]},
             "Cu-catalyzed asymmetric ketone aminoallylation", "enantioselective reductive coupling",
             ["Generic R groups; individual experiments are separate records."], ["C-C", "C-H", "O-H"])
r["participants"] = [participant("ketone", "generic_8", "reactant", (EV[1], EV[3])),
                     participant("allenamide", "11", "reactant", (EV[1], EV[3])),
                     participant("product", "generic_12", "product", (EV[1], EV[3]), product_role="generic protected 1,2-aminoalcohol")]
reactions.append(r)

# Table 1 optimization records.
table = [
    (1, "ss_phbpe", 64, "89:11", "20:80", {}), (2, "r_binap", 82, "83:17", "18:82", {}),
    (3, "r_segphos", 51, "86:14", "30:70", {}), (4, "w3", 64, "90:10", "57:43", {}),
    (5, "w8", 77, "81:19", "93:7", {"er13": "50:50"}), (6, "j6", 58, "91:9", "15:85", {}),
    (7, "j11", 60, "78:22", "28:72", {}), (8, "w8", 71, "87:13", "87:13", {"silane_eq": 10}),
    (9, "w8", 58, ">99:1", "88:12", {"tbuoh": 2}), (10, "w8", 61, "90:10", "91:9", {"tbuoh": 2, "solvent": "phcf3"}),
    (11, "w8", 45, "52:48", "85:15", {"ketone": "8b", "er13": "60:40"}),
    (12, "w8", 50, ">99:1", "97:3", {"ketone": "8b", "tbuoh": 2}),
    (13, "w8", 77, ">99:1", "96:4", {"ketone": "8b", "tbuoh": 2, "solvent": "phcf3"})]
for entry, lig, yld, ratio, er, mods in table:
    ket = mods.get("ketone", "8a")
    branched = "12b" if ket == "8b" else "12a"
    migration = "13b" if ket == "8b" else "13a"
    solv = mods.get("solvent", "toluene")
    sil_eq = mods.get("silane_eq", 2)
    r = reaction(f"table1_entry_{entry:02d}", "specific_experiment", ["optimization"], (EV[2],),
                 {"page": 2, "table": "Table 1", "entry": entry}, "Chiral ligand survey",
                 "Cu-catalyzed reductive aminoallylation", [
                     "Supporting Information with full procedure is not supplied.",
                     "Absolute configuration of the er order is not printed; ratio retained as reported.",
                     "No isolated workup/purification details are available in the supplied main text."],
                 ["C-C", "C-H", "O-H", "C-O"], (GEN,))
    r["participants"] = [
        participant("ketone", ket, "reactant", (EV[2],), charges={"amount": q(0.25, "mmol", "0.25 mmol")}),
        participant("allenamide", "11", "reactant", (EV[2],), charges={"amount": q(0.375, "mmol", "0.375 mmol"), "equivalents": q(1.5, "equiv", "1.5 equiv")}, ref=ket),
        participant("copper", "cuoac2", "catalyst", (EV[2],), charges={"mol_percent": q(5, "mol%", "5 mol %")}, ref=ket),
        participant("ligand", lig, "ligand", (EV[2],), charges={"mol_percent": q(6, "mol%", "6 mol %")}, ref=ket),
        participant("silane", "silane", "reagent", (EV[2],), charges={**({"amount": q(0.50, "mmol", "0.50 mmol")} if sil_eq == 2 else {}), "equivalents": q(sil_eq, "equiv", f"{sil_eq} equiv")}, ref=ket),
        participant("solvent", solv, "solvent", (EV[2],), charges={"volume": q(0.5, "mL", "0.5 mL")}),
        participant("workup", "nh4f", "workup_agent", (EV[2],), stage="workup"),
        participant("branched", branched, "product", (EV[2],), product_role="branched product"),
        participant("migration", migration, "product", (EV[2],), product_role="carbonate-migration product")]
    if mods.get("tbuoh"):
        r["participants"].append(participant("tbuoh", "tbuoh", "additive", (EV[2],), charges={"equivalents": q(2, "equiv", "2 equiv")}, ref=ket))
    r["reaction_scale"] = F("ketone 0.25 mmol", {"reference_molecule_id": mid(ket), "amount": q(0.25, "mmol")}, (EV[2],))
    r["solvent_conditions"]["composition"] = F("toluene" if solv == "toluene" else "PhCF3", [mid(solv)], (EV[2],))
    r["solvent_conditions"]["total_volume"] = F("0.5 mL", q(0.5, "mL"), (EV[2],))
    r["steps"] = [step("reaction", (EV[2],), ["ketone", "allenamide", "copper", "ligand", "silane", "solvent"] + (["tbuoh"] if mods.get("tbuoh") else []),
                       "Cu-catalyzed reductive coupling", "rt", q(24, "h", "24 h")),
                  step("nh4f_workup", (EV[2],), ["workup", "branched", "migration"], "then NH4F", stage="workup")]
    r["outcomes"] = [
        outcome("yield", "yield", [branched], f"{yld}%", q(yld, "%", f"{yld}%"), (EV[2],), "branched product", "1H NMR of the unpurified mixture using dimethyl fumarate as standard", "Table 1 % product", "NMR yield"),
        outcome("product_ratio", "selectivity", [branched, migration], ratio, {"metric": f"{branched}:{migration}", "raw_text": ratio}, (EV[2],), scope=f"{branched}:{migration}", basis="1H NMR of the unpurified mixture"),
        outcome("branched_er", "selectivity", [branched], f"{er} er", {"metric": "er", "major_minor": er}, (EV[2],), scope="reported enantiomer ratio; absolute major-enantiomer descriptor unassigned", basis="chiral HPLC")]
    if mods.get("er13"):
        r["outcomes"].append(outcome("migration_er", "selectivity", [migration], f"{mods['er13']} er", {"metric": "er", "major_minor": mods["er13"]}, (EV[2],), scope="carbonate-migration product er; absolute order unassigned", basis="chiral HPLC"))
    reactions.append(r)

# Scheme 1 substrate scope: 21 Method A experiments and 15 Method B experiments.
scope_a = {
    "a": (75, "93:7"), "b": (73, "67:33"), "c": (90, "95:5"), "d": (67, "72:28"),
    "e": (70, "85:15"), "f": (58, "82:18"), "g": (67, "77:23"), "h": (65, "66:34"),
    "i": (75, "80:20"), "j": (68, "76:24"), "k": (43, "90:10"), "l": (47, "97:3"),
    "m": (73, "92:8"), "n": (54, "96:4"), "o": (55, "92:8"), "p": (76, "92:8"),
    "q": (89, "90:10"), "r": (80, "76:24"), "s": (28, "72:28"), "t": (51, "77:23"),
    "u": (33, None)}
scope_b = {"a": (61, "91:9"), "b": (77, "96:4"), "c": (55, "97:3"), "d": (66, "88:12"),
           "e": (49, "86:14"), "f": (14, "90:10"), "g": (0, None), "h": (40, "83:17"),
           "i": (18, "78:22"), "j": (50, "85:15"), "k": (38, "83:17"), "l": (33, "85:15"),
           "m": (58, "95:5"), "n": (33, "78:22"), "o": (37, "85:15")}


def add_scope(letter, method, yld, er):
    ket = "8" + letter
    depicted = "l12u" if letter == "u" else "12" + letter
    converted = method == "A" and letter in {"b", "c"}
    target = "13" + letter if converted else depicted
    solv = "toluene" if method == "A" else "phcf3"
    ev = (EV[3],)
    r = reaction(f"scheme1_{letter}_method_{method.lower()}", "specific_experiment", ["substrate_scope"], ev,
                 {"page": 3, "scheme": "Scheme 1", "compound": depicted, "method": method},
                 "Enantioselective Cu-catalyzed reductive coupling", "ketone aminoallylation",
                 ["Full experimental procedure, workup and purification are in unprovided SI.",
                  "Yield measurement basis is not stated in the supplied main text." if not (letter == "s" and method == "A") else "Main-product yield basis is unreported; the 58% linear isomer is explicitly described as isolated.",
                  "The reported er does not identify the absolute configuration/order of the major enantiomer."],
                 ["C-C", "C-H", "O-H"] + (["C-O"] if converted else []), (GEN,))
    r["participants"] = [
        participant("ketone", ket, "reactant", ev, charges={"amount": q(0.25, "mmol", "0.25 mmol")}),
        participant("allenamide", "11", "reactant", ev),
        participant("copper", "cuoac2", "catalyst", ev, charges={"mol_percent": q(5, "mol%", "5 mol %")}, ref=ket),
        participant("ligand", "w8", "ligand", ev, charges={"mol_percent": q(6, "mol%", "6 mol %")}, ref=ket),
        participant("silane", "silane", "reagent", ev),
        participant("solvent", solv, "solvent", ev),
        participant("workup", "nh4f", "workup_agent", ev, stage="workup"),
        participant("product", target, "product", ev, product_role="reported scope product")]
    if converted:
        r["participants"].extend([participant("intermediate12", depicted, "product", ev, product_role="pre-conversion 12/13 mixture component"),
                                  participant("nah", "nah", "reagent", ev)])
    if method == "B":
        r["participants"].append(participant("tbuoh", "tbuoh", "additive", ev, charges={"equivalents": q(2, "equiv", "2 equiv")}, ref=ket))
    if letter == "s" and method == "A":
        r["participants"].append(participant("linear", "l12s", "product", ev, product_role="linear regioisomer"))
    r["reaction_scale"] = F("0.25 mmol ketone", {"reference_molecule_id": mid(ket), "amount": q(0.25, "mmol")}, ev)
    r["participants"][1]["charges"]["amount"] = F("Table 1: 11 (0.375 mmol)", q(0.375, "mmol", "0.375 mmol"), (EV[2],), method="inherited_general_condition")
    r["participants"][1]["charges"]["equivalents"] = F("Table 1: 11 (1.5 equiv)", q(1.5, "equiv", "1.5 equiv"), (EV[2],), method="inherited_general_condition")
    r["participants"][1]["equivalence_reference"] = F("Table 1 ketone denominator", mid(ket), (EV[2],), method="inherited_general_condition")
    r["participants"][4]["charges"]["amount"] = F("Table 1: 0.50 mmol Me(MeO)2SiH", q(0.50, "mmol", "0.50 mmol"), (EV[2],), method="inherited_general_condition")
    r["participants"][5]["charges"]["volume"] = F("Table 1: 0.5 mL solvent", q(0.5, "mL", "0.5 mL"), (EV[2],), method="inherited_general_condition")
    r["solvent_conditions"]["composition"] = F("toluene" if method == "A" else "PhCF3", [mid(solv)], ev)
    r["solvent_conditions"]["total_volume"] = F("0.5 mL inherited from optimized Table 1 conditions", q(0.5, "mL"), ev,
                                                       method="inherited_general_condition")
    r["steps"] = [step("reaction", ev, ["ketone", "allenamide", "copper", "ligand", "silane", "solvent"] + (["tbuoh"] if method == "B" else []),
                       f"Method {method} reductive coupling", "rt"),
                  step("nh4f_workup", ev, ["workup", "product"], "NH4F/MeOH", stage="workup")]
    r["steps"][0]["time"] = F("Table 1: 24 h", q(24, "h", "24 h"), (EV[2],), method="inherited_general_condition")
    if converted:
        r["steps"].append(step("convert_mixture_to_13", ev, ["intermediate12", "nah", "product"], "Convert 12/13 mixture to 13 with NaH"))
    r["outcomes"] = [outcome("reported_yield", "yield", [target], f"{yld}%", q(yld, "%", f"{yld}%"), ev,
                                      "reported scope product", None, "Scheme 1 reported yield", "reported yield")]
    if er:
        r["outcomes"].append(outcome("reported_er", "selectivity", [target], f"{er} er", {"metric": "er", "major_minor": er}, ev,
                                             scope="reported er; absolute enantiomer order unassigned"))
    if letter == "s" and method == "A":
        r["outcomes"].append(outcome("linear_isomer_yield", "yield", ["l12s"], "58% yield of linear isomer also isolated", q(58, "%", "58%"), ev,
                                             "linear regioisomer", "isolated", "linear-isomer yield", "isolated yield"))
    reactions.append(r)


for letter, (yld, er) in scope_a.items():
    add_scope(letter, "A", yld, er)
for letter, (yld, er) in scope_b.items():
    add_scope(letter, "B", yld, er)

# Scheme 2 is an author-proposed stereochemical/mechanistic model, not an experiment.
r = reaction("scheme2_proposed_model", "proposed_mechanism", ["mechanistic_study"], (EV[3], EV[4]),
             {"pages": [3, 4], "scheme": "Scheme 2"}, "Working mechanistic stereochemical model",
             "proposed CuH hydrocupration/ketone allylcupration model",
             ["Intermediates 14-16 and TS1 are proposed, not directly observed synthetic experiments.",
              "The (W8)Cu quadrant model is extrapolated from a cited (W1)PdCl2 structure."], None)
reactions.append(r)

# Prose-only failed mechanistic attempt: 12a under LCuOtBu gives 13a only.
r = reaction("prose_12a_lcuotbu_attempt", "specific_experiment", ["control", "mechanistic_study"], (EV[4],),
             {"page": 4, "region": "left-column paragraph before Scheme 3"},
             "Attempted regeneration of intermediate 16 from 12a", "carbonate migration under LCuOtBu treatment",
             ["SI is absent; catalyst loading, solvent, temperature, time and yield are unreported in the supplied main text.",
              "The text reports only that 13a was observed; detection method and detection limit are unreported."], ["C-O"])
r["participants"] = [participant("substrate", "12a", "reactant", (EV[4],)), participant("catalyst", "lcuotbu", "catalyst", (EV[4],)),
                     participant("product", "13a", "product", (EV[4],), product_role="only observed rearrangement product")]
r["steps"] = [step("treatment", (EV[4],), ["substrate", "catalyst", "product"], "treatment with an LCuOtBu catalyst")]
r["outcomes"] = [outcome("observed_product", "other", ["13a"], "only rearrangement product 13a was observed", "observed", (EV[4],),
                                 "rearrangement product", scope="qualitative observation; method/limit unreported")]
reactions.append(r)

# Scheme 3a retro-allylation experiment.
r = reaction("scheme3a_retro_allylation", "specific_experiment", ["mechanistic_study"], (EV[4],),
             {"page": 4, "scheme": "Scheme 3a"}, "Retro-allylation experiment", "Cu-catalyzed retro-allylation/protonation",
             ["Reaction time, atmosphere, workup, purification and yield measurement bases are unreported in the supplied main text.",
              "Outcome percentages are retained without assuming isolated/NMR bases.",
              "70:30 er is for recovered b-17; absolute enantiomer order is not assigned."], None)
r["participants"] = [participant("substrate", "b17", "reactant", (EV[4],)),
                     participant("cui", "cui", "catalyst", (EV[4],), charges={"mol_percent": q(10, "mol%", "10 mol %")}, ref="b17"),
                     participant("kotbu", "kotbu", "reagent", (EV[4],), charges={"mol_percent": q(10, "mol%", "10 mol %")}, ref="b17"),
                     participant("w8", "w8", "ligand", (EV[4],), charges={"mol_percent": q(10, "mol%", "10 mol %")}, ref="b17"),
                     participant("thf", "thf", "solvent", (EV[4],)),
                     participant("p18", "18", "product", (EV[4],), product_role="protonation product"),
                     participant("p19", "19", "product", (EV[4],), product_role="protonation product"),
                     participant("p8a", "8a", "product", (EV[4],), product_role="regenerated ketone"),
                     participant("recovered", "b17", "reactant", (EV[4],), product_role="recovered substrate"),
                     participant("linear", "l17", "product", (EV[4],), product_role="linear allylation product")]
r["solvent_conditions"]["composition"] = F("THF", [mid("thf")], (EV[4],))
r["steps"] = [step("retro_allylation", (EV[4],), ["substrate", "cui", "kotbu", "w8", "thf"], "retro-allylation experiment", "rt")]
r["outcomes"] = [
    outcome("yield18", "yield", ["18"], "23%", q(23, "%", "23%"), (EV[4],), "protonation product"),
    outcome("yield19", "yield", ["19"], "19%", q(19, "%", "19%"), (EV[4],), "protonation product"),
    outcome("yield8a", "yield", ["8a"], "28%", q(28, "%", "28%"), (EV[4],), "regenerated ketone"),
    outcome("recovery_b17", "recovery", ["b17"], "44%", q(44, "%", "44%"), (EV[4],), scope="recovered b-17"),
    outcome("recovered_b17_er", "selectivity", ["b17"], "70:30 er", {"metric": "er", "major_minor": "70:30"}, (EV[4],), scope="recovered b-17 er; starting material was 82:18 er"),
    outcome("yield_l17", "yield", ["l17"], "3%", q(3, "%", "3%"), (EV[4],), "linear allylation product")]
reactions.append(r)

# Scheme 3b scale-up.
r = reaction("scheme3b_scaleup", "specific_experiment", ["scale_up"], (EV[4],),
             {"page": 4, "scheme": "Scheme 3b"}, "1.0 mmol scale reaction", "Cu-catalyzed reductive aminoallylation",
             ["Reaction duration, atmosphere, workup/purification and yield measurement basis are unreported in the supplied main text.",
              "W8 recovery measurement basis is unreported."], ["C-C", "C-H", "O-H"], (GEN,))
r["participants"] = [participant("ketone", "8b", "reactant", (EV[4],), charges={"amount": q(1.0, "mmol", "1.0 mmol")}),
                     participant("allenamide", "11", "reactant", (EV[4],), charges={"equivalents": q(1.5, "equiv", "1.5 equiv")}, ref="8b"),
                     participant("copper", "cuoac2", "catalyst", (EV[4],), charges={"mol_percent": q(5, "mol%", "5 mol %")}, ref="8b"),
                     participant("w8", "w8", "ligand", (EV[4],), charges={"mol_percent": q(6, "mol%", "6 mol %")}, ref="8b"),
                     participant("silane", "silane", "reagent", (EV[4],)), participant("tbuoh", "tbuoh", "additive", (EV[4],)),
                     participant("solvent", "phcf3", "solvent", (EV[4],)), participant("workup", "nh4f", "workup_agent", (EV[4],), stage="workup"),
                     participant("product", "12b", "product", (EV[4],), product_role="branched protected aminoalcohol")]
r["reaction_scale"] = F("1.0 mmol scale", {"reference_molecule_id": mid("8b"), "amount": q(1.0, "mmol")}, (EV[4],))
r["solvent_conditions"]["composition"] = F("PhCF3", [mid("phcf3")], (EV[4],))
r["steps"] = [step("reaction", (EV[4],), ["ketone", "allenamide", "copper", "w8", "silane", "tbuoh", "solvent"], "1.0 mmol scale reductive coupling", "rt"),
              step("workup", (EV[4],), ["workup", "product"], "then NH4F", stage="workup")]
r["outcomes"] = [outcome("yield", "yield", ["12b"], "72%", q(72, "%", "72%"), (EV[4],), "branched product"),
                 outcome("er", "selectivity", ["12b"], "96:4 er", {"metric": "er", "major_minor": "96:4"}, (EV[4],), scope="product er; absolute order unassigned"),
                 outcome("w8_recovery", "recovery", ["w8"], "81% recovery", q(81, "%", "81%"), (EV[4],), scope="W8 recovery")]
reactions.append(r)

# Scheme 3c: one arrow/operation sequence per distinct transformation.
r = reaction("scheme3c_12a_to_13a", "specific_experiment", ["derivatization"], (EV[4],),
             {"page": 4, "scheme": "Scheme 3c", "arrow": "12a to 13a"}, "Carbonate migration and recrystallization",
             "base-promoted carbonate migration", ["Reaction time, temperature and yield measurement basis are not printed.",
             "73% is the yield after hexane/EtOAc crystallization; denominator is not explicitly defined."], ["C-O"])
r["participants"] = [participant("substrate", "12a", "reactant", (EV[4],)), participant("nah", "nah", "reagent", (EV[4],)),
                     participant("thf", "thf", "solvent", (EV[4],)), participant("product", "13a", "product", (EV[4],), product_role="carbonate-migration product"),
                     participant("hex", "hexane", "purification_agent", (EV[4],), stage="purification"), participant("etoac", "etoac", "purification_agent", (EV[4],), stage="purification")]
r["solvent_conditions"]["composition"] = F("THF", [mid("thf")], (EV[4],))
r["steps"] = [step("migration", (EV[4],), ["substrate", "nah", "thf", "product"], "NaH, THF"),
              step("crystallization", (EV[4],), ["product", "hex", "etoac"], "hex/EtOAc crystallization", stage="purification")]
r["outcomes"] = [outcome("migration_yield", "yield", ["13a"], "99%", q(99, "%", "99%"), (EV[4],), "carbonate-migration product", scope="pre-recrystallization reported yield"),
                 outcome("initial_er", "selectivity", ["13a"], "93:7 er", {"metric": "er", "major_minor": "93:7"}, (EV[4],), scope="before recrystallization"),
                 outcome("crystallized_yield", "yield", ["13a"], "hex/EtOAc Xs. 73%", q(73, "%", "73%"), (EV[4],), "recrystallized product", "crystallization/isolation", "yield after crystallization"),
                 outcome("crystallized_er", "selectivity", ["13a"], ">99:1 er", {"metric": "er", "major_minor": ">99:1", "qualifier": "gt"}, (EV[4],), scope="after recrystallization")]
reactions.append(r)

r = reaction("scheme3c_13a_to_20", "specific_experiment", ["derivatization"], (EV[4],),
             {"page": 4, "scheme": "Scheme 3c", "arrow": "13a to 20"}, "Three-step telescoped oxazolidinone removal sequence",
             "tosylation/elimination/enamide hydrolysis", ["Telescoped sequence was performed without intermediate purification.",
             "Solvent, reagent quantities, temperatures, times, workup/purification and yield measurement basis are not reported in the supplied main text."], None)
r["participants"] = [participant("substrate", "13a", "reactant", (EV[4],)), participant("tscl", "tscl", "reagent", (EV[4],)),
                     participant("dbu", "dbu", "reagent", (EV[4],)), participant("nai", "nai", "reagent", (EV[4],)),
                     participant("h2so4", "h2so4", "reagent", (EV[4],)), participant("product", "20", "product", (EV[4],), product_role="three-step telescoped product")]
r["steps"] = [step("tosylation", (EV[4],), ["substrate", "tscl"], "1. TsCl"), step("elimination", (EV[4],), ["dbu", "nai"], "2. DBU/NaI"),
              step("hydrolysis", (EV[4],), ["h2so4", "product"], "3. H2SO4")]
r["outcomes"] = [outcome("three_step_yield", "yield", ["20"], "78% (3 steps)", q(78, "%", "78%"), (EV[4],),
                                 "telescoped product", scope="overall three-step yield from 13a to 20", yield_type="three-step overall yield")]
reactions.append(r)

r = reaction("scheme3c_20_to_21", "specific_experiment", ["derivatization"], (EV[4],),
             {"page": 4, "scheme": "Scheme 3c", "arrow": "20 to 21"}, "Hydrolysis to aminoalcohol 21",
             "oxazolidinone hydrolysis", ["Arrow cites ref. 10a; conditions, time, temperature, workup, purification and yield basis are not reproduced in the supplied main text.",
             "The 86% value is retained for this arrow but the unavailable cited/SI procedure is not inferred."], None)
r["participants"] = [participant("substrate", "20", "reactant", (EV[4],)), participant("product", "21", "product", (EV[4],), product_role="aminoalcohol")]
r["steps"] = [step("hydrolysis", (EV[4],), ["substrate", "product"], "ref. 10a hydrolysis")]
r["procedure"]["general_procedure_reference"] = F("ref. 10a", "reference 10a", (EV[4],))
r["outcomes"] = [outcome("yield", "yield", ["21"], "86%", q(86, "%", "86%"), (EV[4],), "aminoalcohol")]
reactions.append(r)

bundle["molecules"] = list(molecules.values())
bundle["reactions"] = reactions
bundle["source_context"] = [
    {"context_id": f"{DOC}:table1_measurement", "document_id": DOC, "evidence_ids": [EV[2]],
     "kind": "measurement_context", "text": "Table 1 yields and 12:13 ratios were determined by 1H NMR on the unpurified mixture using dimethyl fumarate as standard; er by chiral HPLC. The standard is characterization-stage, not a reaction reagent."},
    {"context_id": f"{DOC}:single_diastereomer", "document_id": DOC, "evidence_ids": [EV[2]],
     "kind": "stereochemical_context", "text": "Table 1 footnote: in all cases, a single diastereomer of 12a and 13a was obtained by 1H NMR analysis; exact dr/detection limit not reported."},
    {"context_id": f"{DOC}:si_notice", "document_id": DOC, "evidence_ids": [EV[4]],
     "kind": "availability_context", "text": "The article states SI exists, but no SI file was supplied or explicitly associated in this workspace."},
    {"context_id": f"{DOC}:references_only", "document_id": DOC, "evidence_ids": [EV[5], EV[6]],
     "kind": "coverage_context", "text": "PDF pages 5 and 6 contain references only; no current-paper synthetic experiment or standalone mass-spectrum table was found."},
]

(RUN / "submitted.json").write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")
(RUN / "raw_model/stereo_diagnostics.json").write_text(json.dumps(stereo_diagnostics, ensure_ascii=False, indent=2) + "\n")
(RUN / "raw_model/extraction_summary.json").write_text(json.dumps({
    "document_id": DOC, "source": "manual visual extraction from original page images",
    "pages_viewed": [1, 2, 3, 4, 5, 6], "reaction_count": len(reactions),
    "molecule_count": len(molecules), "table1_entries": len(table),
    "scheme1_method_a": len(scope_a), "scheme1_method_b": len(scope_b),
    "chemical_accuracy_verified": False,
    "notes": ["Raw source-derived candidate data; standard merge will quarantine all formal records as needs_review.",
              "No SI was associated; missing procedure details remain missing."]}, ensure_ascii=False, indent=2) + "\n")

mention_terms = ["8a", "8b", "11", "12a", "12b", "13a", "13b", "b-17", "18", "19", "20", "21"]
mentions = st.source_mentions(json.loads((RUN / "dataset.json").read_text()), mention_terms, context=220)
(RUN / "raw_model/source_mentions.json").write_text(json.dumps(mentions, ensure_ascii=False, indent=2) + "\n")

# Mark tasks attempted; this is page-level workflow state, not a claim of complete recall.
tasks = json.loads((RUN / "tasks.json").read_text())
for task in tasks:
    if task["document_id"] == DOC:
        task["state"] = "attempted"
        task["reason"] = "Original page image and extracted text inspected; semantic accuracy and completeness remain unverified."
(RUN / "tasks.json").write_text(json.dumps(tasks, ensure_ascii=False, indent=2) + "\n")

page1 = [r["reaction_id"] for r in reactions if r["reaction_id"].startswith(f"{DOC}:fig1")]
page2 = [r["reaction_id"] for r in reactions if ":table1_entry_" in r["reaction_id"]]
page3 = [r["reaction_id"] for r in reactions if ":scheme1_" in r["reaction_id"] or r["reaction_id"].endswith("scheme2_proposed_model")]
page4 = [r["reaction_id"] for r in reactions if any(x in r["reaction_id"] for x in ("prose_12a", "scheme3"))]
coverage = {"document_id": DOC, "semantic_accuracy_verified": False, "pages": [
    {"page": 1, "state": "attempted", "regions": [
        {"region": "Figure 1A-D and graphical abstract", "finding": "Three conceptual/cited-work arrows and one current-work generic scheme; graphical abstract duplicates the current transformation.", "evidence_ids": [EV[1]]},
        {"region": "title, abstract and opening prose", "finding": "Article metadata and background; no additional current-paper concrete experiment.", "evidence_ids": [EV[1], TXT[1]]}],
     "record_ids": page1, "unresolved": ["Generic R-group structures and cited-work procedures are not expanded."]},
    {"page": 2, "state": "attempted", "regions": [
        {"region": "Table 1 entries 1-13 and footnotes a-k", "finding": "Thirteen distinct optimization experiments, including 8b substitutions, excess silane, t-BuOH and PhCF3 variants.", "evidence_ids": [EV[2], TXT[2]]},
        {"region": "body prose around Table 1", "finding": "Author discussion of 12/13 competition and reversible allylcupration; no separate concrete experiment beyond the table on this page.", "evidence_ids": [EV[2], TXT[2]]}],
     "record_ids": page2, "unresolved": ["Full SI procedures and isolated workups unavailable.", "Absolute major-enantiomer identities for reported er values unassigned."]},
    {"page": 3, "state": "attempted", "regions": [
        {"region": "Scheme 1 Method A/B scope and footnotes a-c", "finding": "Twenty-one Method A experiments and fifteen distinct Method B experiments; footnote-b conversions and 12s linear byproduct retained.", "evidence_ids": [EV[3], TXT[3]]},
        {"region": "Scheme 2 and surrounding prose", "finding": "Author-proposed stereochemical/mechanistic model, not a performed synthetic experiment.", "evidence_ids": [EV[3], TXT[3]]}],
     "record_ids": page3, "unresolved": ["Scope yield measurement bases and full procedures are not stated in the supplied main text.", "Dense product SMILES and absolute R/S assignments remain unresolved."]},
    {"page": 4, "state": "attempted", "regions": [
        {"region": "left-column prose before Scheme 3", "finding": "Prose-only attempted 12a/LCuOtBu experiment gave only 13a.", "evidence_ids": [EV[4], TXT[4]]},
        {"region": "Scheme 3a", "finding": "One retro-allylation mechanistic experiment with products 18, 19, 8a, recovered b-17 and l-17.", "evidence_ids": [EV[4], TXT[4]]},
        {"region": "Scheme 3b", "finding": "One 1.0 mmol scale aminoallylation plus W8 recovery outcome.", "evidence_ids": [EV[4], TXT[4]]},
        {"region": "Scheme 3c", "finding": "Three distinct application arrows: 12a to 13a/recrystallization, telescoped 13a to 20, and ref.-10a hydrolysis 20 to 21.", "evidence_ids": [EV[4], TXT[4]]},
        {"region": "Associated Content", "finding": "SI availability notice only; no SI file is associated with this run.", "evidence_ids": [EV[4], TXT[4]]}],
     "record_ids": page4, "unresolved": ["Most Scheme 3 timings, temperatures and measurement bases are absent from main text.", "The 20-to-21 arrow cites ref. 10a without reproducing conditions."]},
    {"page": 5, "state": "attempted", "regions": [
        {"region": "reference list, citations 1-20", "finding": "References only; cited papers were not converted into experiments of this paper.", "evidence_ids": [EV[5], TXT[5]]}],
     "record_ids": [], "unresolved": []},
    {"page": 6, "state": "attempted", "regions": [
        {"region": "continuation of references 20-21", "finding": "References only; no current-work reaction arrow, table or characterization peak list.", "evidence_ids": [EV[6], TXT[6]]}],
     "record_ids": [], "unresolved": []}
]}
(RUN / "coverage/pages_0001_0006.json").write_text(json.dumps(coverage, ensure_ascii=False, indent=2) + "\n")

print(json.dumps({"molecules": len(molecules), "reactions": len(reactions), "stereo_diagnostics": len(stereo_diagnostics)}, ensure_ascii=False))
